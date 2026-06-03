"""
SZ-DBC Phase 2+3 in Julia — HiGHS LP chamber enumeration + integer-key DB check.

Usage: julia --project=<dir> ebe.jl <compact_in.json> <output.json>

Mathematica runs Phase 1 (condition extraction) and exports a compact structure.
When feasible_sigmas is empty, this script runs:
  Phase 2 — BFS over the hyperplane arrangement using HiGHS LP feasibility checks.
            Non-strict negation (>= for False strict conditions) bridges all octants.
            Degenerate boundary patterns are filtered algebraically after BFS.
  Phase 3 — For each genuine chamber, groups DB terms by integer exponent coefficient
            vector and checks that each group sums to zero (exact Rational{Int64}).

When feasible_sigmas is non-empty (precomputed, e.g. reused from D4 EBE check),
Phase 2 is skipped and only Phase 3 runs.

Phase 2 correctness (LP encoding):
  sigma[i]=1 with strict condition:     condEffLhs[i]·J >= eps  (>0 in the interior)
  sigma[i]=1 with non-strict condition: condEffLhs[i]·J >= 0
  sigma[i]=0 with strict condition:    -condEffLhs[i]·J >= 0    (>=0, includes boundary)
  sigma[i]=0 with non-strict condition:-condEffLhs[i]·J >= eps  (>0, strict negation)
  Degenerate filter: strict pairs (both False) and non-strict pairs (both True) removed.

Phase 3 correctness:
  Grouping by integer coefficient vector (v_coeffs + energy_coeffs[state]) is
  algebraically exact: two terms cancel iff their combined exponent function
  dot(v, J) is identical for ALL J, which happens iff the integer vectors are equal.
"""

using HiGHS, JSON3, Printf

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

    # Condition data for Phase 2
    cond_eff_lhs = Vector{Int64}[
        Int64[Int64(x) for x in v] for v in data.cond_eff_lhs]
    cond_is_strict = Bool[Bool(x) for x in data.cond_is_strict]
    initial_sigma  = Int[Int(x) for x in data.initial_sigma]

    state_energy_coeffs = [
        Int64[Int64(data.state_energy_coeffs[s][j]) for j in 1:n_atoms]
        for s in 1:n_states]

    # Feasible sigma patterns (empty → Julia runs Phase 2)
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
            cond_eff_lhs, cond_is_strict, initial_sigma,
            state_energy_coeffs, feasible_sigmas,
            unique_weights, leaf_weight_idx,
            pair_states, pair_ij_srcs, pair_ji_srcs)
end

# ---------------------------------------------------------------------------
# Phase 2: BFS over hyperplane arrangement using HiGHS LP
# ---------------------------------------------------------------------------

"""
Test whether a sigma pattern is feasible using a HiGHS LP.

For each condition i:
  sigma[i]=1, strict:     condEffLhs[i]·J >= eps   (interior of positive halfspace)
  sigma[i]=1, non-strict: condEffLhs[i]·J >= 0
  sigma[i]=0, strict:    -condEffLhs[i]·J >= 0     (non-strict negation ≥ 0)
  sigma[i]=0, non-strict:-condEffLhs[i]·J >= eps   (strict negation > 0)

The non-strict representation for sigma=0 of a strict condition is what allows
the BFS to visit degenerate boundary patterns and bridge disconnected octants.
"""
function is_feasible(sigma       :: Vector{Int},
                     cond_eff_lhs :: Vector{Vector{Int64}},
                     cond_is_strict :: Vector{Bool},
                     n_atoms     :: Int,
                     eps         :: Float64) :: Bool
    isempty(sigma) && return true  # k=0: no conditions, always feasible
    h = HiGHS.Highs_create()
    HiGHS.Highs_setBoolOptionValue(h, "output_flag", false)
    HiGHS.Highs_setBoolOptionValue(h, "solver_output_flag", false)
    for _ in 1:n_atoms
        HiGHS.Highs_addVar(h, -1e10, 1e10)
    end
    inds = Int32[j-1 for j in 1:n_atoms]
    for i in eachindex(sigma)
        v    = cond_eff_lhs[i]
        sign = sigma[i] == 1 ? 1.0 : -1.0
        coeffs = Float64[sign * v[j] for j in 1:n_atoms]
        # lb = eps when the constraint is strict in the active direction:
        #   sigma=1 + strict condition  →  condEffLhs·J > 0  →  lb = eps
        #   sigma=1 + non-strict        →  condEffLhs·J ≥ 0  →  lb = 0
        #   sigma=0 + strict (non-strict negation ≥0) →  lb = 0
        #   sigma=0 + non-strict (strict negation >0) →  lb = eps
        lb = ((sigma[i] == 1) == cond_is_strict[i]) ? eps : 0.0
        HiGHS.Highs_addRow(h, lb, 1e30, n_atoms, inds, coeffs)
    end
    HiGHS.Highs_run(h)
    status = HiGHS.Highs_getModelStatus(h)
    HiGHS.Highs_destroy(h)
    return status == 7  # 7 = kOptimal (feasible LP)
end

"""
Remove degenerate sigma patterns (those lying on a measure-zero boundary hyperplane).
For each contradictory pair (i,j) with condEffLhs[i] == -condEffLhs[j]:
  - Strict pairs (both Less/Greater):           filter when sigma[i]==0 && sigma[j]==0
  - Non-strict pairs (both LessEqual/GreaterEqual): filter when sigma[i]==1 && sigma[j]==1
  - Mixed pairs: no degenerate boundary exists; skip.
"""
function filter_degenerate(feasible_sigmas :: Vector{Vector{Int}},
                            cond_eff_lhs   :: Vector{Vector{Int64}},
                            cond_is_strict :: Vector{Bool},
                            n_conds        :: Int) :: Vector{Vector{Int}}
    contra_ff = Tuple{Int,Int}[]
    contra_tt = Tuple{Int,Int}[]
    for i in 1:n_conds, j in (i+1):n_conds
        if cond_eff_lhs[i] == -cond_eff_lhs[j]
            if cond_is_strict[i] && cond_is_strict[j]
                push!(contra_ff, (i, j))
            elseif !cond_is_strict[i] && !cond_is_strict[j]
                push!(contra_tt, (i, j))
            end
        end
    end
    isempty(contra_ff) && isempty(contra_tt) && return feasible_sigmas

    n_before = length(feasible_sigmas)
    filtered = filter(feasible_sigmas) do sigma
        !any(sigma[i] == 0 && sigma[j] == 0 for (i, j) in contra_ff) &&
        !any(sigma[i] == 1 && sigma[j] == 1 for (i, j) in contra_tt)
    end
    n_removed = n_before - length(filtered)
    if n_removed > 0
        @printf(stderr, "  EBE: %d genuine open chamber(s) (removed %d degenerate boundary patterns)\n",
                length(filtered), n_removed)
    end
    filtered
end

"""
BFS over the hyperplane arrangement to find all genuine open chambers.
Starts from initial_sigma (an interior point of one chamber) and explores
neighbours by flipping one bit at a time, testing feasibility via HiGHS LP.
After BFS, applies the degenerate filter.
"""
function run_phase2(initial_sigma  :: Vector{Int},
                    cond_eff_lhs   :: Vector{Vector{Int64}},
                    cond_is_strict :: Vector{Bool},
                    n_atoms        :: Int,
                    n_conds        :: Int) :: Vector{Vector{Int}}
    if n_conds == 0
        @printf(stderr, "  EBE: 1 feasible region (no branch conditions)\n")
        return [Int[]]
    end

    eps     = 1e-6
    visited = Set{Vector{Int}}()
    push!(visited, copy(initial_sigma))
    feasible = [copy(initial_sigma)]
    queue    = [copy(initial_sigma)]

    while !isempty(queue)
        sigma = popfirst!(queue)
        for i in 1:n_conds
            sigma2    = copy(sigma)
            sigma2[i] = 1 - sigma2[i]
            sigma2 in visited && continue
            push!(visited, copy(sigma2))
            if is_feasible(sigma2, cond_eff_lhs, cond_is_strict, n_atoms, eps)
                push!(feasible, copy(sigma2))
                push!(queue, copy(sigma2))
            end
        end
    end

    @printf(stderr, "  EBE: %d feasible region(s)  (%d sign patterns visited via BFS)\n",
            length(feasible), length(visited))

    filter_degenerate(feasible, cond_eff_lhs, cond_is_strict, n_conds)
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
  1. Pre-compute which Case applies to each unique weight once per chamber.
  2. Group by integer coefficient vector (v_coeffs + energy_coeffs) rather
     than by the Rational dot-product value.  Algebraically exact; no BigInt.

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

    active_terms = Vector{Vector{Term}}(undef, n_weights)
    key_buf      = Vector{Int64}(undef, n_atoms)

    for (reg_idx, sigma) in enumerate(feasible_sigmas)
        for wi in 1:n_weights
            active_terms[wi] = get_terms(unique_weights[wi], sigma)
        end

        for p in 1:n_pairs
            isempty(pair_ij_srcs[p]) && isempty(pair_ji_srcs[p]) && continue

            empty!(groups)
            i, j = pair_states[p]
            ei   = state_energy_coeffs[i]
            ej   = state_energy_coeffs[j]

            for leaf_idx in pair_ij_srcs[p]
                for term in active_terms[leaf_weight_idx[leaf_idx]]
                    iszero(term.c_num) && continue
                    @inbounds @. key_buf = term.v_coeffs + ei
                    c = Rational{Int64}(term.c_num, term.c_den)
                    key = copy(key_buf)
                    groups[key] = get(groups, key, zero(Rational{Int64})) + c
                end
            end

            for leaf_idx in pair_ji_srcs[p]
                for term in active_terms[leaf_weight_idx[leaf_idx]]
                    iszero(term.c_num) && continue
                    @inbounds @. key_buf = term.v_coeffs + ej
                    c = Rational{Int64}(term.c_num, term.c_den)
                    key = copy(key_buf)
                    groups[key] = get(groups, key, zero(Rational{Int64})) - c
                end
            end

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
         cond_eff_lhs, cond_is_strict, initial_sigma,
         state_energy_coeffs, feasible_sigmas,
         unique_weights, leaf_weight_idx,
         pair_states, pair_ij_srcs, pair_ji_srcs) = load_compact(ARGS[1])
    end
    @printf(stderr, "Load: %.2fs  (%d conds, %d states, %d pairs, %d leaves, %d unique weights)\n",
            t_load, n_conds, n_states, n_pairs, n_leaves, length(unique_weights))

    # Phase 2: enumerate chambers via HiGHS LP BFS (if not precomputed)
    t_p2 = 0.0
    if isempty(feasible_sigmas)
        t_p2 = @elapsed begin
            feasible_sigmas = run_phase2(initial_sigma, cond_eff_lhs, cond_is_strict,
                                         n_atoms, n_conds)
        end
        @printf(stderr, "Phase 2: %.3fs  (%d genuine chambers)\n", t_p2, length(feasible_sigmas))
    else
        @printf(stderr, "Phase 2: skipped (using %d precomputed chambers)\n",
                length(feasible_sigmas))
    end

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
