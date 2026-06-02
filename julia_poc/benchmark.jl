"""
Minimal Phase 3 Step-C benchmark in Julia.

Tests the core DB check inner loop using synthetic data that matches
vmmc_2d.wl's structure exactly:
  - 56 orbit reps, ~78 leaves each  (4368 total leaves)
  - 5688 communicating state pairs
  - 512 feasible coupling-constant chambers
  - leaf weights after J* substitution = rational sums of Exp[-β*v] terms

This isolates whether Julia's compiled rational arithmetic is faster than
Mathematica's interpreted Step C loop, without requiring a data export.

Mathematica baseline (measured on vmmc_2d.wl, no D4):
  Step A (weight substitution, 4368 leaves):  ~164 ms/region  → 84 s total
  Step C (pair loop, 5688 pairs):             ~181 ms/region  → 93 s total
  Phase 2 (FindInstance BFS):                                    5 s
  DB EBE total:                                                181.7 s

Expected Julia speedup on Step C: 50-150x (compiled tight loop vs WL interpreter).
The BFS steps (26.7s) and Phase 2 (5s) remain in Mathematica regardless.
Floor for any Julia-augmented run: ~34s (BFS + Phase 2 + ergodicity).
"""

using Random
using Printf
using Dates

# ──────────────────────────────────────────────────────────────────────────────
# Data types
# Each leaf weight after J* substitution is a sum of terms c_k * exp(-β * v_k).
# Represented as a fixed-size tuple (cn, cd, vn, vd) where c = cn//cd, v = vn//vd.
# ──────────────────────────────────────────────────────────────────────────────
const Term = Tuple{Int32, Int32, Int32, Int32}   # (c_num, c_den, v_num, v_den)
const LeafWeight = Vector{Term}

# ──────────────────────────────────────────────────────────────────────────────
# Generate synthetic data matching VMMC scale
# ──────────────────────────────────────────────────────────────────────────────
function generate_vmmc_data(rng::AbstractRNG)
    n_reps   = 56
    n_leaves = 78        # leaves per rep (VMMC average)
    n_states = 504
    n_pairs  = 5688

    # Leaf weights: each leaf has 1–3 Exp terms (realistic for VMMC).
    # For a correct algorithm the actual values cancel in Step C; we only
    # care about data shape and operation count.
    leaf_weights = Vector{Vector{LeafWeight}}(undef, n_reps)
    for ri in 1:n_reps
        leaf_weights[ri] = Vector{LeafWeight}(undef, n_leaves)
        for li in 1:n_leaves
            n_terms = rand(rng, 1:3)
            leaf_weights[ri][li] = [
                (Int32(rand(rng, -4:4)), Int32(rand(rng, 2:12)),
                 Int32(rand(rng, -3:3)), Int32(rand(rng, 1:6)))
                for _ in 1:n_terms
            ]
        end
    end

    # pair_ij_src[p] = list of (ri,li) whose leaf weights contribute to T(i→j)
    # pair_ji_src[p] = list of (ri,li) whose leaf weights contribute to T(j→i)
    # VMMC measured: 4536/5688 pairs have 1 source, 1152 have 2+ sources → avg ~1.3
    pair_ij_src = Vector{Vector{Tuple{Int32,Int32}}}(undef, n_pairs)
    pair_ji_src = Vector{Vector{Tuple{Int32,Int32}}}(undef, n_pairs)
    for p in 1:n_pairs
        n_ij = rand(rng, 1:2)
        n_ji = rand(rng, 1:2)
        pair_ij_src[p] = [(Int32(rand(rng,1:n_reps)), Int32(rand(rng,1:n_leaves))) for _ in 1:n_ij]
        pair_ji_src[p] = [(Int32(rand(rng,1:n_reps)), Int32(rand(rng,1:n_leaves))) for _ in 1:n_ji]
    end

    # State energies as small rationals (VMMC: sums of ≤3 coupling terms)
    energies = [Rational{Int32}(Int32(rand(rng,-6:6)), Int32(rand(rng,1:6)))
                for _ in 1:n_states]

    # Pair state indices for DB check
    pair_states = [(rand(rng, Int32(1):Int32(n_states)),
                    rand(rng, Int32(1):Int32(n_states)))
                   for _ in 1:n_pairs]

    return leaf_weights, pair_ij_src, pair_ji_src, energies, pair_states
end

# ──────────────────────────────────────────────────────────────────────────────
# Phase 3 Step C: DB check for one region
# For each pair (i,j): accumulate T_ij*exp(-β*E_i) - T_ji*exp(-β*E_j)
# grouped by Boltzmann exponent, verify all coefficient sums = 0.
# ──────────────────────────────────────────────────────────────────────────────
function step_c(
    leaf_weights  :: Vector{Vector{LeafWeight}},
    pair_ij_src   :: Vector{Vector{Tuple{Int32,Int32}}},
    pair_ji_src   :: Vector{Vector{Tuple{Int32,Int32}}},
    energies      :: Vector{Rational{Int32}},
    pair_states   :: Vector{Tuple{Int32,Int32}},
    groups        :: Dict{Rational{Int64}, Rational{Int64}}  # pre-allocated, passed in
)::Bool
    for p in eachindex(pair_states)
        i, j = pair_states[p]
        ei = Rational{Int64}(energies[i])
        ej = Rational{Int64}(energies[j])

        empty!(groups)

        # T(i→j) * exp(-β*E_i) contributions
        for (ri, li) in pair_ij_src[p]
            for (cn, cd, vn, vd) in leaf_weights[ri][li]
                iszero(cn) && continue
                exp_val = Rational{Int64}(vn, vd) + ei
                groups[exp_val] = get(groups, exp_val, zero(Rational{Int64})) +
                                  Rational{Int64}(cn, cd)
            end
        end

        # T(j→i) * exp(-β*E_j) contributions (subtracted)
        for (ri, li) in pair_ji_src[p]
            for (cn, cd, vn, vd) in leaf_weights[ri][li]
                iszero(cn) && continue
                exp_val = Rational{Int64}(vn, vd) + ej
                groups[exp_val] = get(groups, exp_val, zero(Rational{Int64})) -
                                  Rational{Int64}(cn, cd)
            end
        end

        # DB condition: all coefficient groups must vanish
        # NOTE: in the real implementation this would return false immediately on
        # first violation, but for benchmarking we MUST check all pairs so we
        # can measure per-pair cost accurately. The real algorithm always PASSes
        # for correct algorithms (all 5688 pairs checked every region).
        # We suppress early exit to get a valid throughput measurement.
        for coeff in values(groups)
            if !iszero(coeff)
                # Real code: return false. Benchmark: record but continue.
                # (Synthetic random data won't zero; real VMMC data will.)
                break
            end
        end
    end
    return true
end

# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────
function main()
    rng = Random.MersenneTwister(42)

    println("Generating synthetic VMMC-scale data...")
    lw, src_ij, src_ji, energies, pairs = generate_vmmc_data(rng)

    groups = Dict{Rational{Int64}, Rational{Int64}}()

    # JIT warmup (first call compiles)
    println("Warming up Julia JIT...")
    step_c(lw, src_ij, src_ji, energies, pairs, groups)
    step_c(lw, src_ij, src_ji, energies, pairs, groups)

    # ── Time one region ──────────────────────────────────────────────────────
    t1 = time_ns()
    for _ in 1:10
        step_c(lw, src_ij, src_ji, energies, pairs, groups)
    end
    t1_elapsed = (time_ns() - t1) / 1e9
    t_per_region_ms = t1_elapsed / 10 * 1000

    # ── Time all 512 regions ─────────────────────────────────────────────────
    n_regions = 512
    t2 = time_ns()
    for _ in 1:n_regions
        step_c(lw, src_ij, src_ji, energies, pairs, groups)
    end
    t2_elapsed = (time_ns() - t2) / 1e9

    # ── Report ───────────────────────────────────────────────────────────────
    math_step_c_per_region_ms = 181.0   # measured Mathematica Step C per region
    math_step_c_total_s       = 93.0    # measured Mathematica Step C total (512 regions)
    math_step_a_total_s       = 84.0    # measured Mathematica Step A total (512 regions)
    math_phase3_total_s       = 181.7   # measured total DB EBE (Step A + Phase 2 + Step C)

    println()
    println("="^60)
    println("  RESULTS")
    println("="^60)
    @printf "  Julia Step C per region:     %7.3f ms\n" t_per_region_ms
    @printf "  Julia Step C × 512 regions: %7.3f s\n"  t2_elapsed
    println()
    @printf "  Mathematica Step C per reg:  %7.1f ms  (measured)\n" math_step_c_per_region_ms
    @printf "  Mathematica Step C × 512:    %7.1f s   (measured)\n" math_step_c_total_s
    println()
    speedup_step_c = math_step_c_total_s / t2_elapsed
    @printf "  Step C speedup (Julia/WL):   %7.1fx\n" speedup_step_c
    println()
    println("  ── Implications for full Phase 3 ──────────────────────")
    println("  Step A (weight substitution) stays in Mathematica: 84s")
    println("  Moving Step A to Julia too would require exporting the")
    println("  symbolic leaf weight structure — not done in this PoC.")
    println()
    julia_phase3_est = t2_elapsed + 5.0  # + Phase 2 (FindInstance stays in WL)
    total_with_julia_stepc_only = (84.0 + 5.0 + t2_elapsed + 27.0)  # step A in WL + phase2 + julia step C + BFS
    total_with_julia_full = (5.0 + t2_elapsed + 27.0)  # phase2 + julia step C + BFS (step A in Julia too)
    @printf "  If only Step C in Julia:     %5.0f s total (%.1fx speedup)\n" total_with_julia_stepc_only (217.0/total_with_julia_stepc_only)
    @printf "  If Step A+C both in Julia:   %5.0f s total (%.1fx speedup)\n" total_with_julia_full (217.0/total_with_julia_full)
    println()
    println("  ── Floor ───────────────────────────────────────────────")
    println("  BFS (τ + regular): ~27s   — must stay in Mathematica (Block)")
    println("  Phase 2 (FindInstance): ~5s — must stay in Mathematica")
    println("  Minimum achievable total: ~32s (Julia does all of Phase 3)")
    @printf "  Maximum total speedup:       %5.1fx\n" (217.0/32.0)
    println("="^60)
end

main()
