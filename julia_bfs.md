# Julia τ-BFS: Feasibility Analysis and Proof of Concept

**Date:** 2026-06-03  
**Context:** SZ-DBC currently delegates EBE Phase 2+3 to Julia (HiGHS LP + integer-key DB check). This document investigates whether the remaining Mathematica phase — the τ-BFS (translational invariance check + random-path enumeration) — can also move to Julia, eliminating the Mathematica dependency entirely.

**PoC code:** [`julia_poc/tau_bfs_poc.jl`](julia_poc/tau_bfs_poc.jl)

---

## Background: What the τ-BFS Does

Mathematica's `$dbcCheckTranslational` and `$dbcBuildStateLeaves` perform three things simultaneously in a single BFS pass:

1. **τ-check** — verifies translational invariance by augmenting all particle positions with symbolic offsets `(r + τr, c + τc)` and checking that no path weight depends on τr or τc.
2. **Path enumeration** — exhaustively explores the algorithm's random-number tree by replacing `RandomReal`, `RandomInteger`, `RandomChoice`, etc. with interceptors that read from a fixed bit sequence and track the exact rational path probability.
3. **Leaf weight production** — produces symbolic leaf weights of the form `c × exp(−β × ΔE)` (where `c` is a rational selection probability and `ΔE` is a linear combination of coupling constants) that feed into EBE Phases 1–3.

The mechanism relies on Mathematica's `Block[{RandomReal = ..., RandomInteger = ..., ...}, alg[state]]`, which overrides random functions globally within any nested call depth for the duration of the expression. This is the single feature that is hardest to replicate in Julia.

---

## What Was Tested

Six scenarios were implemented and run in [`julia_poc/tau_bfs_poc.jl`](julia_poc/tau_bfs_poc.jl):

### Scenario A — TauPos: algebraic τ-cancellation

A `TauPos` struct represents a particle position as `(vr + tr·τr, vc + tc·τc)` where `vr`, `vc`, `tr`, `tc` are all plain integers. Initial τ-augmented positions have `tr = tc = 1`. All arithmetic is defined on this struct.

**Key property:** When computing pairwise distances, subtraction yields `(vr_a − vr_b, 0, 0)` — the τ coefficients cancel exactly because both particles carry the same offset. This is an algebraic, not numerical, check: if τ ever leaked into a distance computation, `pbc_d2` would throw an exception immediately.

**Result:** PASS. ΔE computed on τ-augmented states equals ΔE on non-augmented states. τ cancels provably in all pairwise distances for any configuration.

---

### Scenario B — BitSeqRNG: exact rational path tracking

A `BitSeqRNG <: AbstractRNG` reads from a fixed bit sequence and tracks the exact `Rational{Int64}` path weight. Two distinct read operations are implemented, matching Mathematica's behaviour:

- `read_bit!` — used for `RandomInteger` calls: increments position AND multiplies weight by `1//2`.
- `accept_test!` — used for `RandomReal[] < p` acceptance tests: increments position directly (no ×½), then multiplies weight by `condP` (accept) or `1 − condP` (reject).

**Critical bug found during testing:** calling `read_bit!` from inside `accept_test!` applied an extra ×½ factor, causing total leaf weights to sum to 1/2 instead of 1. The fix — incrementing `pos` directly in `accept_test!` — exactly matches Mathematica's `acceptTestI`, which calls `pos++` (not `readBit[]`) before multiplying by the conditional probability.

**Result:** PASS. 6 leaves for a 3-way choice + binary acceptance; weights sum to exactly `1//1`.

---

### Scenario C — Single Metropolis τ-BFS

`single_metropolis.wl` was translated to Julia (~60 lines). The Julia version:
- Uses `rand_choice!(rng, n)` for uniform selection (reads ⌈log2(n)⌉ bits, matches `RandomChoice`).
- Uses `new_real!(rng)` + `accept_test!(rng, j, p_rat)` for the Metropolis acceptance step.
- Computes ΔE as a sparse `Dict{NTuple{3,Int}, Int}` — a DeltaE vector over coupling atom indices.
- Substitutes `τ → 0` on all leaf next-states (equivalent to Mathematica's `./{τr→0, τc→0}`).

**Results (nGrid=3, 3 particles, packed seed state `(1,1)(1,2)(1,3)`):**

| Metric | Value |
|--------|-------|
| Leaves | 24 |
| Weights sum | `1//1` ✓ |
| Next-states τ-free after substitution | PASS ✓ |
| Unique ΔE branch conditions extracted | 6 |
| τ-invariance (pbc_d2 would throw on leak) | PASS ✓ |
| Leaf count matches non-τ BFS | PASS ✓ |
| Time per representative state | **~3 ms** |

Sample conditions extracted (the `allConds` equivalent):
```
cond 1: -1×J(1,2,1) + 1×J(1,2,2)
cond 2: -1×J(1,2,1) + 1×J(1,2,2) - 1×J(1,3,1) + 1×J(1,3,2)
cond 3: -1×J(1,3,1) + 1×J(1,3,2)
```
These are the correct linear-combination branch conditions (ΔE vectors) that Mathematica would extract from the Piecewise clauses.

**Timing comparison:**

| Phase | Mathematica | Julia |
|-------|-------------|-------|
| Per-state BFS | ~200 ms | ~3 ms |
| All ~12 reps (single_metropolis) | ~2 s | ~0.04 s |
| Speedup | — | **~50–100×** |

---

### Scenario D — Broken algorithm: τ-non-invariance detection

An algorithm that selects a particle with probability proportional to its absolute row index (a τ-non-invariant selection) was tested. Detection relies on explicitly checking `p.tr != 0` whenever a position is used as a selection weight.

**Result:** PASS. Violation correctly detected in the τ-augmented run; no false positive in the τ=0 run.

---

### Scenario F — Concrete vs. always-read tree structure

This test exposes the key difference between concrete-coupling and symbolic-coupling BFS trees.

**Packed seed (all downhill):** With all couplings equal (J=0.3), or with distinct but all-repulsive couplings from a packed state, all valid moves disperse particles → ΔE ≤ 0 always → p=1. The concrete shortcut (p ≥ 1 → accept, no bit) and the always-read policy both produce **24 leaves**.

**Spread seed (uphill moves exist):** From a spread-out initial state `(1,1)(2,2)(3,3)`, some moves bring particles closer → ΔE > 0 → p < 1 → accept/reject branching. Results:

| Policy | Leaves |
|--------|--------|
| Concrete shortcut (downhill → no bit) | 24 |
| Always-read (Mathematica-equivalent) | **42** |

Both sum to weight 1. The 18 extra leaves are the reject branches for uphill moves. **The always-read policy matches Mathematica's symbolic tree structure.** Without it, downhill moves collapse to a single accept path, producing fewer leaves and a different branch-condition structure.

---

## Findings

### What Julia can do

| Capability | Status | Notes |
|-----------|--------|-------|
| τ-cancellation check in pairwise distances | ✅ Algebraically exact | TauPos throws on any τ leak |
| Exact rational path-weight tracking | ✅ Correct | `Rational{Int64}`, weights sum to 1 |
| BFS path enumeration | ✅ Correct | Matches Mathematica leaf counts |
| τ-non-invariance detection (selection weights) | ✅ Working | Explicit `tr != 0` flag |
| Branch condition extraction (ΔE vectors) | ✅ Working | Sparse coupling-atom Dict |
| Symbolic tree structure (always-read policy) | ✅ Implementable | 42 vs 24 leaves on spread seed |
| Speed vs Mathematica | ✅ ~50–100× faster | 3 ms vs ~200 ms per rep |
| No Symbolics.jl required for τ-check | ✅ Confirmed | Pure integer arithmetic |

### What Julia cannot do (or requires significant work)

**1. Automatic algorithm translation — the dominant practical barrier.**

Mathematica's `Block[{RandomReal = ...}, alg[state]]` is a dynamic binding that intercepts all random calls anywhere in the call stack during evaluation. Julia has no language-level equivalent. There is no way to run a `.wl` algorithm file in Julia without manually rewriting it as a Julia function that accepts a `BitSeqRNG` and uses `rand_choice!`, `new_real!`, `accept_test!` in place of Mathematica's random functions. Every algorithm file needs a manual translation.

Effort estimate per algorithm:
- `single_metropolis`, `kawasaki`: 1–2 hours each (simple logic, few random calls)
- `vmmc_2d`: 2–4 days (cluster formation, multi-step acceptance, complex state transformations)
- `quadratic_field`, `broken_*`: 1–3 hours each

**2. Always-read policy is required but non-trivial for complex algorithms.**

To match Mathematica's symbolic tree structure, acceptance tests must always read a bit even for downhill moves (p=1). The always-read policy discards the "impossible" reject path (weight=0) via `OutOfRangeException`. This works for single_metropolis but requires careful handling in algorithms where the Piecewise structure is more complex (e.g., VMMC's partial unpairing step where multiple conditions compound).

**3. Exact symbolic leaf weights for DB check require DeltaE-vector storage.**

To produce leaf weights in the compact format needed for `ebe.jl` Phase 2+3, each leaf must store:
- A rational coefficient `c_num/c_den` (product of selection probabilities — already in `rng.weight`)
- A ΔE vector (already extracted as the branch condition)

Currently the PoC uses a concrete float `p = min(1, exp(−β × ΔE_concrete))` for the acceptance test. For a complete DB-check integration, the acceptance weight must be stored as `{ΔE_vector, accepted}` rather than a concrete float. This is achievable — the DeltaE struct is already in place — but requires restructuring how leaf weights are accumulated.

**4. Exotic τ-violation coverage.**

The TauPos approach guarantees detection of τ-violations that flow through pairwise distance computations and through explicit position-dependent selection weights. It would NOT automatically catch:
- τ-violations in sorting/ordering logic (e.g., picking the particle with smallest absolute row)
- τ-violations through external lookup tables indexed by absolute position
- τ-violations in algorithms where the random decision depends on a position modulo something that doesn't cancel

In Mathematica, every expression is genuinely symbolic: τr and τc propagate through any function call and their presence in any path weight is detected. Julia's coverage is proportional to how thoroughly the algorithm translation instruments τ-tracking. For well-written physical algorithms (energy-based, PBC-consistent), this risk is low. For arbitrary algorithms, it is a real concern.

**5. Ergodicity check, D4 symmetry, orbit computation not implemented.**

The τ-BFS feeds into orbit expansion, ergodicity checking, and D4 verification. Replicating those phases would require translating additional Mathematica infrastructure: `$dbcComputeOrbits`, `$dbcCheckErgodicityFromLeaves`, `$dbcApplyGElem`, etc. This is straightforward algorithmic translation but adds further scope.

---

## Path to a Fully Julia-Native Checker

The PoC demonstrates that a complete Julia implementation is architecturally sound. The components are:

```
Algorithm (Julia)
  ↓  uses BitSeqRNG
τ-BFS (Julia)
  → TauPos + pbc_d2 (τ-check)
  → always-read BitSeqRNG (path enumeration)
  → DeltaE extraction (allConds equivalent)
  → tau_sub_zero on next-states
  ↓
ebe.jl (already in Julia)
  Phase 2: HiGHS LP BFS over hyperplane arrangement
  Phase 3: integer-key DB check
```

The only element not yet implemented is the DeltaE-to-compact-JSON export step, which would replace `$dbcExportCompact` in Mathematica. Given that `ebe.jl` already reads the compact JSON format, this would require:
1. Producing the `cond_eff_lhs`, `cond_is_strict`, `initial_sigma` fields from BFS leaves
2. Producing `unique_weights`, `leaf_weight_idx`, `pair_states`, `pair_ij_srcs`, `pair_ji_srcs` from the leaf weight structure
3. Emitting JSON

This is non-trivial (200–400 lines) but well-defined. The existing `$dbcExportCompact` implementation in `dbc_core.wl` is the specification.

**Expected total time for a fully Julia single_metropolis end-to-end:** ~50 ms (vs current ~14 s with Mathematica + Julia). The 5 s Julia cold-start overhead would also vanish since the checker would be a single Julia process.

---

## Conclusion

**Julia can fully replicate the τ-BFS phase.** The algebraic infrastructure — `TauPos`, `BitSeqRNG`, `DeltaE` vectors — is proven correct in the PoC, is faster by 50–100×, requires no symbolic algebra library, and produces the correct branch conditions.

**The bottleneck is not the technology but the algorithm translation effort.** Each `.wl` algorithm file must be manually rewritten as a Julia function using the `BitSeqRNG` API. For the current set of examples this is a bounded but non-trivial task (roughly 1–4 days total depending on algorithm complexity). For future algorithms, users would write Julia directly rather than Mathematica, which is the natural long-term position for an open-source tool.

**A Julia-native algorithm API is the recommended path forward if full open-source independence is a priority.** If the goal is purely performance and Mathematica remains acceptable, the current architecture (Mathematica for τ-BFS + Julia for Phase 2+3) already achieves the major speedups and the incremental gain from moving τ-BFS to Julia is modest (~2 s saved on vmmc_2d). The architectural gain — eliminating the Mathematica dependency entirely — is the stronger argument.
