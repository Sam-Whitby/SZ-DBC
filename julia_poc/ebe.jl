"""
SZ-DBC Phase 3 in Julia — compact structure approach.

Usage: julia --project=<dir> ebe.jl <compact_in.json> <output.json>

Mathematica runs Phase 1 (condition extraction) and Phase 2 (non-strict BFS +
degenerate filter), exports a compact structure with:
  - 'feasible_sigmas': the genuine open chambers found by Phase 2
  - per-unique-weight compact Piecewise structure (sigma-substitution)

This script runs Phase 3 only: for each sigma pattern, pre-computes active
Boltzmann terms, groups them by integer exponent coefficient vector, and checks
that the DB residual is exactly zero for every communicating state pair.

Correctness:
  Grouping by integer coefficient vector (v_coeffs + energy_coeffs[state]) is
  algebraically exact: two terms cancel iff their combined exponent function
  dot(v, J) is identical for ALL J, which happens iff the integer vectors are
  equal.  This is equivalent to Mathematica's \$dbcIsExpZero grouping and does
  not require evaluating at any specific J*.

  The 'feasible_sigmas' from Mathematica's Phase 2 cover all genuine open
  chambers of the hyperplane arrangement.  Phase 2 uses a non-strict BFS that
  bridges disconnected octants through degenerate boundary patterns, then filters
  those patterns, leaving only chambers with positive measure.
"""

using JSON3, Printf

# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

struct Term
    c_num    :: Int64
    c_den    :: Int64
    v_coeffs :: Vector{Int64}
end

struct Case
    active_sigma :: Vector{Int64}   # 0/1 per active condition
    terms        :: Vector{Term}
end

struct UniqueWeight
    active_cond_idxs :: Vector{Int64}   # 0-based indices into conditions
    cases            :: Vector{Case}
end

# ---------------------------------------------------------------------------
# JSON loading
# ---------------------------------------------------------------------------

function load_compact(path::String)
    data = open(path, "r") do f; JSON3.read(f); end

    n_atoms  = Int(data.n_atoms)
    n_conds  = Int(data.n_conds)
    n_states = Int(data.n_states)
    n_pairs  = Int(data.n_pairs)
    n_leaves = Int(data.n_leaves)

    state_energy_coeffs = [
        Int64[Int64(data.state_energy_coeffs[s][j]) for j in 1:n_atoms]
        for s in 1:n_states]

    # Feasible sigma patterns from Mathematica Phase 2
    feasible_sigmas = [Int[Int(x) for x in s] for s in data.feasible_sigmas]

    function parse_terms(raw_terms)
        Term[Term(
            Int64(t.c_num),
            Int64(t.c_den),
            Int64[Int64(x) for x in t.v_coeffs])
            for t in raw_terms]
    end

    function parse_case(rc)
        Case(
            Int64[Int64(x) for x in rc.active_sigma],
            parse_terms(rc.terms))
    end

    unique_weights = UniqueWeight[
        UniqueWeight(
            Int64[Int64(x) for x in uw.active_cond_idxs],
            Case[parse_case(c) for c in uw.cases])
        for uw in data.unique_weights]

    leaf_weight_idx = Int64[Int64(x) + 1 for x in data.leaf_weight_idx]  # 1-based

    pair_states  = Tuple{Int64,Int64}[(Int64(p[1]), Int64(p[2])) for p in data.pair_states]
    pair_ij_srcs = [Int64[Int64(x)+1 for x in s] for s in data.pair_ij_srcs]  # 1-based
    pair_ji_srcs = [Int64[Int64(x)+1 for x in s] for s in data.pair_ji_srcs]

    return (n_atoms, n_conds, n_states, n_pairs, n_leaves,
            state_energy_coeffs, feasible_sigmas,
            unique_weights, leaf_weight_idx,
            pair_states, pair_ij_srcs, pair_ji_srcs)
end

# ---------------------------------------------------------------------------
# Phase 3: DB check using compact structure
# ---------------------------------------------------------------------------

"""
Look up the Boltzmann terms for a unique weight under a given sigma pattern.
Projects sigma to the weight's active conditions and finds the matching case.
Avoids allocations: uses element-wise comparison without building a new vector.
"""
function get_terms(uw::UniqueWeight, sigma01::Vector{Int}) :: Vector{Term}
    @inbounds for c in uw.cases
        match = true
        for (j, idx) in enumerate(uw.active_cond_idxs)
            if sigma01[idx + 1] != c.active_sigma[j]
                match = false; break
            end
        end
        match && return c.terms
    end
    return Term[]  # unreachable for valid export data
end

"""
Run Phase 3 DB check over all chambers.

Key optimisations:
  1. Pre-compute which Case applies to each unique weight once per chamber
     (n_chambers × n_weights lookups instead of one per pair-leaf).
  2. Group by integer coefficient vector (v_coeffs + energy_coeffs) rather
     than by the Rational dot-product value.  The grouping is algebraically
     exact: two terms belong to the same Boltzmann class iff their combined
     integer coefficient vector is identical for all J.  This eliminates all
     BigInt arithmetic from Phase 3.

Returns (pass::Bool, violations::Vector{Tuple{Int,Int,Int}}).
"""
function run_phase3(
    feasible_sigmas     :: Vector{Vector{Int}},
    unique_weights      :: Vector{UniqueWeight},
    leaf_weight_idx     :: Vector{Int64},
    pair_states         :: Vector{Tuple{Int64,Int64}},
    pair_ij_srcs        :: Vector{Vector{Int64}},
    pair_ji_srcs        :: Vector{Vector{Int64}},
    state_energy_coeffs :: Vector{Vector{Int64}},
    n_atoms             :: Int)

    n_pairs   = length(pair_states)
    n_weights = length(unique_weights)

    violations = Tuple{Int,Int,Int}[]
    groups     = Dict{Vector{Int64}, Rational{Int64}}()

    # Per-chamber buffer: active terms for each unique weight
    active_terms = Vector{Vector{Term}}(undef, n_weights)
    # Reusable key buffer (avoids allocation in hot loop)
    key_buf      = Vector{Int64}(undef, n_atoms)

    for (reg_idx, sigma) in enumerate(feasible_sigmas)
        # Pre-compute which terms are active for each unique weight this chamber.
        for wi in 1:n_weights
            active_terms[wi] = get_terms(unique_weights[wi], sigma)
        end

        for p in 1:n_pairs
            isempty(pair_ij_srcs[p]) && isempty(pair_ji_srcs[p]) && continue

            empty!(groups)
            i, j = pair_states[p]
            ei   = state_energy_coeffs[i]
            ej   = state_energy_coeffs[j]

            # T(i→j) × Exp(-β·E_i): add contributions
            for leaf_idx in pair_ij_srcs[p]
                for term in active_terms[leaf_weight_idx[leaf_idx]]
                    iszero(term.c_num) && continue
                    @inbounds @. key_buf = term.v_coeffs + ei
                    c = Rational{Int64}(term.c_num, term.c_den)
                    key = copy(key_buf)
                    groups[key] = get(groups, key, zero(Rational{Int64})) + c
                end
            end

            # T(j→i) × Exp(-β·E_j): subtract contributions
            for leaf_idx in pair_ji_srcs[p]
                for term in active_terms[leaf_weight_idx[leaf_idx]]
                    iszero(term.c_num) && continue
                    @inbounds @. key_buf = term.v_coeffs + ej
                    c = Rational{Int64}(term.c_num, term.c_den)
                    key = copy(key_buf)
                    groups[key] = get(groups, key, zero(Rational{Int64})) - c
                end
            end

            # Check all group coefficients are zero
            for coeff in values(groups)
                if !iszero(coeff)
                    push!(violations, (i, j, reg_idx))
                    break
                end
            end
        end
    end

    return isempty(violations), violations
end

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

function main()
    if length(ARGS) != 2
        println(stderr, "Usage: julia --project=<dir> ebe.jl <input.json> <output.json>")
        exit(1)
    end

    t_load = @elapsed begin
        (n_atoms, n_conds, n_states, n_pairs, n_leaves,
         state_energy_coeffs, feasible_sigmas,
         unique_weights, leaf_weight_idx,
         pair_states, pair_ij_srcs, pair_ji_srcs) = load_compact(ARGS[1])
    end
    @printf(stderr, "Load: %.2fs  (%d conds, %d states, %d pairs, %d leaves, %d unique weights, %d chambers)\n",
            t_load, n_conds, n_states, n_pairs, n_leaves, length(unique_weights), length(feasible_sigmas))

    # Phase 3: DB check
    t_p3 = @elapsed begin
        pass, violations = run_phase3(
            feasible_sigmas, unique_weights, leaf_weight_idx,
            pair_states, pair_ij_srcs, pair_ji_srcs,
            state_energy_coeffs, n_atoms)
    end
    @printf(stderr, "Phase 3: %.3fs  result=%s\n", t_p3, pass ? "PASS" : "FAIL")

    open(ARGS[2], "w") do f
        if pass
            write(f, "{\"pass\":true}")
        else
            vlist = [[v[1], v[2], v[3]]
                     for v in violations[1:min(length(violations), 10)]]
            JSON3.write(f, Dict("pass"=>false, "violations"=>vlist))
        end
    end
end

main()
