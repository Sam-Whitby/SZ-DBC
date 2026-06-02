"""
SZ-DBC Phase 3 Step C in Julia.

Usage: julia --project=<dir> phase3.jl <input.json> <output.json>

Mathematica does Step A (substituting J* into leaf weights for every feasible
region) and exports rational (c, v) pairs where weight = c * exp(-beta * v).
Julia performs only Step C: for each region and each communicating state pair,
sum the leaf-weight contributions, form T(i->j)*exp(-b*Ei) - T(j->i)*exp(-b*Ej),
group terms by Boltzmann exponent, and verify all rational coefficient sums = 0.

This is identical to Mathematica's dbcIsExpZero but compiled (~30-100x faster).
The BFS and FindInstance phases remain in Mathematica (Mathematica-only features).
"""

using JSON3
using Printf

# One leaf weight term: (c_num, c_den, v_num, v_den)  →  c * exp(-beta * v)
const Term = NTuple{4, Int64}

# Step C: DB check for all regions given pre-evaluated leaf weights
function run_step_c(
    # leaf_weights[region][leaf_idx] = Vector{Term}
    leaf_weights  :: Vector{Vector{Vector{Term}}},
    pair_states   :: Vector{NTuple{2,Int}},
    pair_ij_srcs  :: Vector{Vector{Int}},   # flat leaf indices for T(i→j)
    pair_ji_srcs  :: Vector{Vector{Int}},   # flat leaf indices for T(j→i)
    # state_energies[region][state_idx] = (e_num, e_den)
    state_energies :: Vector{Vector{NTuple{2,Int64}}}
)
    n_pairs   = length(pair_states)
    n_regions = length(leaf_weights)

    violations = Tuple{Int,Int,Int}[]
    groups = Dict{Rational{Int64}, Rational{Int64}}()

    for reg in 1:n_regions
        lw = leaf_weights[reg]
        en = state_energies[reg]

        for p in 1:n_pairs
            i, j = pair_states[p]
            ei = Rational{Int64}(en[i][1], en[i][2])
            ej = Rational{Int64}(en[j][1], en[j][2])

            isempty(pair_ij_srcs[p]) && isempty(pair_ji_srcs[p]) && continue

            empty!(groups)

            # T(i→j) * exp(-beta * E_i)
            @inbounds for leaf_idx in pair_ij_srcs[p]
                for (cn, cd, vn, vd) in lw[leaf_idx]
                    iszero(cn) && continue
                    key = Rational{Int64}(vn, vd) + ei
                    groups[key] = get(groups, key, zero(Rational{Int64})) +
                                  Rational{Int64}(cn, cd)
                end
            end

            # T(j→i) * exp(-beta * E_j)  [subtract]
            @inbounds for leaf_idx in pair_ji_srcs[p]
                for (cn, cd, vn, vd) in lw[leaf_idx]
                    iszero(cn) && continue
                    key = Rational{Int64}(vn, vd) + ej
                    groups[key] = get(groups, key, zero(Rational{Int64})) -
                                  Rational{Int64}(cn, cd)
                end
            end

            for coeff in values(groups)
                if !iszero(coeff)
                    push!(violations, (i, j, reg))
                    break
                end
            end
        end
    end

    return isempty(violations), violations
end

function load_data(path::String)
    data = open(path, "r") do f
        JSON3.read(f)
    end

    n_regions = Int(data.n_regions)
    n_leaves  = Int(data.n_leaves)
    n_states  = Int(data.n_states)
    n_pairs   = Int(data.n_pairs)

    # leaf_weights[region][leaf_idx] = Vector{Term}
    leaf_weights = Vector{Vector{Vector{Term}}}(undef, n_regions)
    for (rg, region_data) in enumerate(data.leaf_weights)
        region_lw = Vector{Vector{Term}}(undef, n_leaves)
        for (lk, leaf_data) in enumerate(region_data)
            region_lw[lk] = Term[(Int64(t[1]), Int64(t[2]), Int64(t[3]), Int64(t[4]))
                                  for t in leaf_data]
        end
        leaf_weights[rg] = region_lw
    end

    pair_states  = NTuple{2,Int}[(Int(p[1]), Int(p[2])) for p in data.pair_states]
    pair_ij_srcs = [Int[Int(x) for x in s] for s in data.pair_ij_srcs]
    pair_ji_srcs = [Int[Int(x) for x in s] for s in data.pair_ji_srcs]

    # state_energies[region][state] = (e_num, e_den)
    state_energies = Vector{Vector{NTuple{2,Int64}}}(undef, n_regions)
    for (rg, reg_en) in enumerate(data.state_energies)
        state_energies[rg] = NTuple{2,Int64}[(Int64(e[1]), Int64(e[2])) for e in reg_en]
    end

    return leaf_weights, pair_states, pair_ij_srcs, pair_ji_srcs, state_energies
end

function main()
    if length(ARGS) != 2
        println(stderr, "Usage: julia --project=<dir> phase3.jl <input.json> <output.json>")
        exit(1)
    end

    t_load = @elapsed begin
        lw, pst, pij, pji, se = load_data(ARGS[1])
    end

    # JIT warmup on a small subset
    run_step_c(lw[1:1], pst[1:min(20,end)], pij[1:min(20,end)],
               pji[1:min(20,end)], se[1:1])

    t_p3 = @elapsed begin
        pass, violations = run_step_c(lw, pst, pij, pji, se)
    end

    open(ARGS[2], "w") do f
        if pass
            write(f, "{\"pass\":true}")
        else
            vlist = [[v[1], v[2], v[3]] for v in violations[1:min(length(violations),10)]]
            JSON3.write(f, Dict("pass"=>false, "violations"=>vlist))
        end
    end

    @printf(stderr,
        "Julia Phase 3 (Step C only): load=%.2fs  step_c=%.3fs  regions=%d  pairs=%d  result=%s\n",
        t_load, t_p3, length(lw), length(pst), pass ? "PASS" : "FAIL")
end

main()
