"""
SZ-DBC Phase 2 + Phase 3 in Julia — compact structure approach.

Usage: julia --project=<dir> ebe.jl <compact_in.json> <output.json>

Mathematica exports a compact Phase-1 structure (Piecewise condition
indices + coefficient vectors) for each unique leaf weight.  This script:

  Phase 2 — Chamber enumeration via CDDLib BFS (exact rational LP)
             Finds all feasible coupling-constant chambers by BFS over
             sign patterns, using CDDLib for strict-interior-point queries.

  Phase 3 — Exact DB check using compact evaluation
             For each chamber and each communicating pair, evaluates leaf
             weights from the compact structure (dot products with jStar)
             and checks the Boltzmann-exponent grouped DB residual = 0.

Both phases are fully deterministic: all arithmetic is exact rational (BigInt).
"""

using CDDLib, Polyhedra, LinearAlgebra, JSON3, Printf

# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

struct Term
    c_num   :: Int64
    c_den   :: Int64
    v_coeffs :: Vector{Int64}
end

struct Clause
    true_cond_idxs :: Vector{Int64}   # 0-based indices into cond_eff_lhs
    terms          :: Vector{Term}
end

struct UniqueWeight
    clauses       :: Vector{Clause}
    default_terms :: Vector{Term}
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

    cond_eff_lhs = Rational{BigInt}[
        Rational{BigInt}(data.cond_eff_lhs[i][j])
        for i in 1:n_conds, j in 1:n_atoms]

    state_energy_coeffs = [
        Int64[Int64(data.state_energy_coeffs[s][j]) for j in 1:n_atoms]
        for s in 1:n_states]

    initial_sigma = Int[Int(x) for x in data.initial_sigma]
    initial_jstar = Rational{BigInt}[
        Rational{BigInt}(Int64(data.initial_jstar[j][1]),
                         Int64(data.initial_jstar[j][2]))
        for j in 1:n_atoms]

    function parse_terms(raw_terms)
        Term[Term(
            Int64(t.c_num),
            Int64(t.c_den),
            Int64[Int64(x) for x in t.v_coeffs])
            for t in raw_terms]
    end

    function parse_clause(rc)
        Clause(
            Int64[Int64(x) for x in rc.true_cond_idxs],
            parse_terms(rc.terms))
    end

    unique_weights = UniqueWeight[
        UniqueWeight(
            Clause[parse_clause(c) for c in uw.clauses],
            parse_terms(uw.default_terms))
        for uw in data.unique_weights]

    leaf_weight_idx = Int64[Int64(x) + 1 for x in data.leaf_weight_idx]  # 1-based

    pair_states  = Tuple{Int64,Int64}[(Int64(p[1]), Int64(p[2])) for p in data.pair_states]
    pair_ij_srcs = [Int64[Int64(x)+1 for x in s] for s in data.pair_ij_srcs]  # 1-based
    pair_ji_srcs = [Int64[Int64(x)+1 for x in s] for s in data.pair_ji_srcs]

    return (n_atoms, n_conds, n_states, n_pairs, n_leaves,
            cond_eff_lhs, state_energy_coeffs,
            initial_sigma, initial_jstar,
            unique_weights, leaf_weight_idx,
            pair_states, pair_ij_srcs, pair_ji_srcs)
end

# ---------------------------------------------------------------------------
# Phase 2: Chamber enumeration
# ---------------------------------------------------------------------------

const R = Rational{BigInt}

"""
Build the H-representation polytope for a given sigma pattern.
Mirrors Mathematica's FindInstance exactly:
  sigma[i]=1: condition is True  → eff_lhs[i]·J >= 0  (non-strict, b = 0)
  sigma[i]=0: condition is False → eff_lhs[i]·J <  0  (strict,    b = -1)
The strict False constraint excludes the origin, so the polytope is non-trivial.
"""
function make_chamber_poly(sigma01::Vector{Int}, eff_lhs::Matrix{R}, lib)
    k, n = size(eff_lhs)
    rows = Vector{R}[]; bs = R[]
    for i in 1:k
        if sigma01[i] == 1
            # True: eff_lhs * J >= 0  →  -eff_lhs * J <= 0
            push!(rows, -eff_lhs[i,:]); push!(bs, R(0))
        else
            # False: eff_lhs * J < 0  →  eff_lhs * J <= -1  (strict via -1)
            push!(rows,  eff_lhs[i,:]); push!(bs, R(-1))
        end
    end
    A = reduce(vcat, [reshape(r, 1, length(r)) for r in rows])
    polyhedron(hrep(A, bs), lib)
end

"""
Check strict feasibility: does the sigma pattern have a strictly interior point?
Returns (is_feasible, poly).
"""
function chamber_feasible(sigma01, eff_lhs, lib)
    poly = make_chamber_poly(sigma01, eff_lhs, lib)
    return !isempty(poly), poly
end

"""
Get a rational point for a feasible chamber using CDDLib vertex enumeration.
CDDLib's vertices are on the boundary of the feasible region, which is fine
for Mathematica compatibility: the DB check only uses J* for substitution,
not for strictness testing.
Returns nothing if polytope is empty (no feasible J).
"""
function get_interior_point(sigma01::Vector{Int}, eff_lhs::Matrix{R}, poly)
    vr  = vrep(poly)
    pts = collect(Polyhedra.points(vr))
    return isempty(pts) ? nothing : first(pts)
end

"""
BFS over sign patterns to enumerate all feasible chambers.
Starts from the Mathematica-provided (initial_sigma, initial_jstar) as the
first chamber (even if it has no strictly interior point — jStar on boundary
is still a valid evaluation point for the DB check).
Returns list of (sigma, jStar) pairs.
"""
function enumerate_chambers(
    eff_lhs::Matrix{R},
    initial_sigma::Vector{Int},
    initial_jstar::Vector{R},
    lib)

    k, n = size(eff_lhs)
    chambers = Tuple{Vector{Int}, Vector{R}}[]
    visited  = Dict{Vector{Int}, Bool}()

    # Always include the initial (Mathematica) chamber as the first entry
    push!(chambers, (initial_sigma, initial_jstar))
    visited[initial_sigma] = true
    queue = [initial_sigma]

    n_checked = 0
    while !isempty(queue)
        sigma = popfirst!(queue)
        for i in 1:k
            sigma2 = copy(sigma)
            sigma2[i] = 1 - sigma2[i]
            haskey(visited, sigma2) && continue
            n_checked += 1
            feas, poly = chamber_feasible(sigma2, eff_lhs, lib)
            if feas
                J2 = get_interior_point(sigma2, eff_lhs, poly)
                if J2 !== nothing
                    visited[sigma2] = true
                    push!(chambers, (sigma2, J2))
                    push!(queue, sigma2)
                else
                    visited[sigma2] = false
                end
            else
                visited[sigma2] = false
            end
        end
    end

    @printf(stderr, "Phase 2: %d chambers  (%d patterns checked)\n",
            length(chambers), n_checked)
    return chambers
end

# ---------------------------------------------------------------------------
# Phase 3: DB check using compact structure
# ---------------------------------------------------------------------------

"""
Find the matching clause for a given sigma pattern and unique weight.
Uses first-match Piecewise semantics: first clause whose true_cond_idxs
are all satisfied (sigma[idx+1]==1) is used; falls back to default_terms.
"""
function get_terms(uw::UniqueWeight, sigma01::Vector{Int}) :: Vector{Term}
    for clause in uw.clauses
        if all(sigma01[idx + 1] == 1 for idx in clause.true_cond_idxs)
            return clause.terms
        end
    end
    return uw.default_terms
end

"""
Evaluate Boltzmann exponent: v = dot(v_coeffs, jStar).
"""
@inline function eval_v(v_coeffs::Vector{Int64}, jStar::Vector{R}) :: R
    s = R(0)
    for j in eachindex(v_coeffs)
        iszero(v_coeffs[j]) && continue
        s += R(v_coeffs[j]) * jStar[j]
    end
    s
end

"""
Run Phase 3 DB check over all chambers.
Returns (pass::Bool, violations::Vector{Tuple{Int,Int,Int}}).
"""
function run_phase3(
    chambers       :: Vector{Tuple{Vector{Int}, Vector{R}}},
    unique_weights :: Vector{UniqueWeight},
    leaf_weight_idx :: Vector{Int64},
    pair_states    :: Vector{Tuple{Int64,Int64}},
    pair_ij_srcs   :: Vector{Vector{Int64}},
    pair_ji_srcs   :: Vector{Vector{Int64}},
    state_energy_coeffs :: Vector{Vector{Int64}},
    n_atoms        :: Int)

    n_pairs   = length(pair_states)
    n_regions = length(chambers)

    violations = Tuple{Int,Int,Int}[]
    groups     = Dict{R, R}()

    for (reg_idx, (sigma, jStar)) in enumerate(chambers)
        # Energy values per state (dot product with jStar)
        energy_jstar = [eval_v(state_energy_coeffs[s], jStar)
                        for s in eachindex(state_energy_coeffs)]

        for p in 1:n_pairs
            isempty(pair_ij_srcs[p]) && isempty(pair_ji_srcs[p]) && continue

            empty!(groups)
            i, j = pair_states[p]

            # T(i→j) × Exp(-β·E_i): add contributions
            for leaf_idx in pair_ij_srcs[p]
                uw = unique_weights[leaf_weight_idx[leaf_idx]]
                for term in get_terms(uw, sigma)
                    iszero(term.c_num) && continue
                    v   = eval_v(term.v_coeffs, jStar)
                    key = v + energy_jstar[i]
                    c   = R(term.c_num, term.c_den)
                    groups[key] = get(groups, key, R(0)) + c
                end
            end

            # T(j→i) × Exp(-β·E_j): subtract contributions
            for leaf_idx in pair_ji_srcs[p]
                uw = unique_weights[leaf_weight_idx[leaf_idx]]
                for term in get_terms(uw, sigma)
                    iszero(term.c_num) && continue
                    v   = eval_v(term.v_coeffs, jStar)
                    key = v + energy_jstar[j]
                    c   = R(term.c_num, term.c_den)
                    groups[key] = get(groups, key, R(0)) - c
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
         cond_eff_lhs, state_energy_coeffs,
         initial_sigma, initial_jstar,
         unique_weights, leaf_weight_idx,
         pair_states, pair_ij_srcs, pair_ji_srcs) = load_compact(ARGS[1])
    end
    @printf(stderr, "Load: %.2fs  (%d conds, %d states, %d pairs, %d leaves, %d unique weights)\n",
            t_load, n_conds, n_states, n_pairs, n_leaves, length(unique_weights))

    lib = CDDLib.Library(:exact)

    # Phase 2: chamber enumeration
    t_p2 = @elapsed begin
        chambers = if n_conds == 0 || n_atoms == 0
            # Trivial: single chamber with empty jStar
            [(Int[], R[])]
        else
            enumerate_chambers(cond_eff_lhs, initial_sigma, initial_jstar, lib)
        end
    end
    @printf(stderr, "Phase 2: %.2fs\n", t_p2)

    # Phase 3: DB check
    t_p3 = @elapsed begin
        pass, violations = run_phase3(
            chambers, unique_weights, leaf_weight_idx,
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
