# SZ-DBC: Technical Analysis

Timing measurements, bottleneck analysis, and evaluation of potential optimisations.
All timings are measured on an Apple M-series Mac (8 logical cores, Julia 1.12.6,
Mathematica 14.3) running `vmmc_2d.wl` with `$particleTypes={1,2,3}`, `$nGrid=3`
unless otherwise noted.

---

## 1. Workflow and Timing Breakdown

### Pipeline overview

The checker runs eight sequential phases for every algorithm file.

**Step 1 — State enumeration** (`$dbcEnumerateNParticleStates`): generates all
`N!/(n₁!n₂!...)` placements on the `nGrid²` torus. For `{1,2,3}` on 3×3: 504 states.
~0.1s for all examples.

**Step 2 — Orbit computation** (`$dbcComputeOrbits`): groups states into translation
orbits, reducing 504 states to 56 orbit representatives for vmmc_2d. ~0.1s.

**Step 4 — τ-BFS** (`$dbcCheckTranslational`): runs the symbolic BFS engine
(`RunWithBitsAT`) on each orbit rep with positions augmented by `{τr, τc}`. Every
`RandomReal/RandomInteger/RandomChoice` call is intercepted and replaced with a
symbolic bit-reader; each leaf accumulates a path weight as a Piecewise expression in
the coupling constants. After the BFS, checks that all leaf weights are τ-free.
Converts τ-BFS leaves to concrete leaves by substituting τ→0, avoiding a second BFS
pass (saving ~13s on vmmc_2d).

| Example | Orbit reps | τ-BFS time |
|---|---|---|
| kawasaki, broken\_\* | 3–7 | ~0.5–2s |
| single\_metropolis | ~9 | ~4s |
| vmmc\_2d | 56 | ~16s |

The dominant cost is the depth of the decision tree per rep. VMMC has a cluster-building
inner loop generating trees of depth ~12 bits (3 `RandomReal` calls per candidate
neighbour, up to 4 candidates per cluster step). Each of 56 reps takes ~0.3s.

**Step 5 — BFS leaves**: reused from τ-BFS at zero cost (τ→0 substitution in
τ-BFS leaves gives exactly the same output as a separate BFS from the same rep).

**Step 7 — Ergodicity**: graph BFS from the seed state over the orbit-expanded
transition graph. ~2s for vmmc_2d, ~0.1s for simple examples.

**Step 8 Phase 1 — Condition extraction**: scans all leaf weights and collects
distinct Piecewise conditions via `Cases[lf[[3]], HoldPattern[Piecewise[cl_,_]] :> ...]`.
For vmmc_2d: k=12 conditions from 4368 leaves across 56 reps. ~0.1s.

**Step 8 Phase 2 — Hyperplane arrangement BFS** (`FindInstance` calls): BFS over
2^k sigma patterns, calling `FindInstance[constraints, allSymParams, Rationals, 1]` to
test feasibility of each. For vmmc_2d: 2048 patterns visited, 512 feasible, filtered to
216 genuine chambers after degenerate-boundary removal. Measured: 6 ms per feasible
call, 2 ms per infeasible, ~10s total.

**Step 8 Phase 3 — DB check per chamber**:

- *Mathematica path*: substitutes rational J\* into all leaf weights, forms
  `T(i→j)·Exp[-β·E(i)] - T(j→i)·Exp[-β·E(j)]`, groups by `Exp[-β·rᵢ]` exponent class
  via `$dbcIsExpZero`. 216 chambers × 5688 pairs = 1.23 M pair evaluations. ~195s.
- *Julia path (`-julia`)*: reads compact structure, groups by integer coefficient vector
  `v_coeffs + energy_coeffs[state]`, checks rational coefficient sums.
  0.7s for 216 chambers.

### Full timing table

| Phase | kawasaki / broken\_\* | single\_metropolis | vmmc\_2d | vmmc\_2d -julia |
|---|---|---|---|---|
| State enum + orbits | 0.1s | 0.2s | 0.2s | 0.2s |
| τ-BFS | 0.5–2s | ~4s | ~16s | ~16s |
| Ergodicity | 0.1s | 1s | 2s | 2s |
| EBE Phase 2 | <0.1s | ~2s | ~10s | ~10s |
| EBE Phase 3 | <0.5s | ~15s | ~195s | ~0.7s |
| Export + Julia startup | — | — | — | ~2.5s |
| **Total** | **~3–5s** | **~35s** | **~226s** | **~33s** |

---

## 2. Current Bottlenecks

**Without `-julia` (vmmc_2d)**:
1. Phase 3 (195s, 86%) — `$dbcIsExpZero` does symbolic work on 1.23 M pair evaluations
   with β as a free symbol; dominant cost.
2. τ-BFS (16s, 7%) — symbolic BFS of 56 complex VMMC trees.
3. Phase 2 BFS (10s, 4%) — 2048 `FindInstance` calls at 2–6 ms each.

**With `-julia` (vmmc_2d)**:
1. τ-BFS (16s, 48%) — now the hard floor.
2. Phase 2 BFS (10s, 30%) — becomes the second-largest cost.
3. Julia startup + JSON export (2.5s, 8%).
4. Julia Phase 3 (0.7s, 2%) — essentially free.

For simpler examples, τ-BFS dominates everything; Phase 3 is negligible.

---

## 3. Effect of Setting β = 1 Everywhere

The short answer: no speedup for the Julia path, and it would break exact
verification on the Mathematica path.

**Why β matters in Mathematica Phase 3**: `$dbcIsExpZero` relies on the fact that
`{Exp[-β·r₁], Exp[-β·r₂], ...}` are linearly independent functions of β for distinct
rationals rᵢ. It checks that the coefficient of each distinct `Exp[-β·rᵢ]` sums to
zero — an algebraic proof that DB holds for *all* β simultaneously. With β=1, exponents
become concrete floats and grouping degenerates to a floating-point sum: verification
at one β value only, breaking exactness.

**Julia Phase 3 is already β-free**: `$dbcCompactTerms` extracts `vCoeffs` by
computing `Coefficient[arg /. β→1, symParamList[[j]]]`. β is substituted out at export
time; Julia never sees β at runtime. It groups by integer v\_coeffs + energy\_coeffs and
checks rational sums — correct for all β and all J simultaneously. Setting β=1 in the
algorithm file changes nothing in Julia's runtime.

**In the τ-BFS**: β appears in VMMC accept weights like `1 - Exp[β*(E_init - E_fwd)]`.
These stay symbolic during the BFS because the BFS explores branching structure, not
weight values. The Piecewise conditions that determine chamber structure are comparisons
between coupling-constant expressions (`eInit < eFwd`), which are β-free. Setting β=1
would produce slightly simpler leaf weight expressions but would not change the branching
tree, the conditions extracted, or the chamber count.

**Summary**: β=1 produces no speedup for either path and breaks exact verification for
the Mathematica Phase 3 path. The compact export already removes β implicitly at the
coefficient-extraction step.

---

## 4. Running Everything in Julia

### What is and is not feasible

| Phase | Feasible in Julia? | Notes |
|---|---|---|
| State enumeration | Yes, trivially | Not a bottleneck |
| Orbit computation | Yes | Not a bottleneck |
| τ-BFS | **No (fundamental)** | Requires symbolic evaluation of arbitrary Mathematica code |
| Phase 2 BFS (FindInstance) | **Yes, with HiGHS** | Measured 19× per-call speedup |
| Phase 3 DB check | Yes — already implemented | 280× measured speedup |
| Ergodicity | Yes, trivially | Not a bottleneck |

### τ-BFS cannot move to Julia

The τ-BFS runs the user's `Algorithm[state_]` function symbolically by intercepting
`RandomReal`, `RandomInteger`, and `RandomChoice` via Mathematica's `Block` mechanism.
It accumulates Piecewise leaf weights as symbolic expressions in the coupling constants.
This requires Mathematica's symbolic evaluation engine. Moving it to Julia would require
either a Mathematica → Julia compiler for arbitrary algorithm code, or redesigning the
algorithm description format as a Julia-embeddable DSL — neither is practical.

### Phase 2 in Julia: measured numbers

| | Per call | 2048 calls | Cold start | Total |
|---|---|---|---|---|
| Mathematica `FindInstance` | 6 ms feasible / 2 ms infeasible | ~10s | (included) | ~10s |
| Julia HiGHS LP (warm) | **0.32 ms** | **0.65s** | 2.22s (cached) | **2.9s** |
| Julia HiGHS LP (sysimage) | 0.32 ms | 0.65s | ~0.3s | **~1s** |

HiGHS is 19× faster per LP call. With a PackageCompiler.jl sysimage, cold-start
drops from 2.22s to ~0.3s (comparable to JSON3-only startup).

Since Phase 2 and Phase 3 would run in the same Julia subprocess (Phase 1 exports
conditions JSON → Julia does Phase 2 → exports feasible sigmas JSON → Julia does Phase 3),
the 2.22s cold start is paid once. Combined EBE in one Julia process:

```
Mathematica Phase 1 export:   ~0.1s  (conditions + compact weights → JSON)
Julia cold start:             ~2.2s  (or ~0.3s with sysimage)
Julia Phase 2 (HiGHS):        ~0.7s  (2048 LP calls, includes BFS logic)
Julia Phase 3:                ~0.7s  (216 chambers × 5688 pairs)
Total EBE:                    ~3.7s  vs current 13.2s (-julia) or 205s (pure Mathematica)
```

The BFS ordering constraint (sigma₂ cannot be queued until its parent sigma is tested)
means Phase 2 cannot be fully parallelised, but within each BFS level all k=12 neighbour
patterns can be tested in parallel using `Threads.@threads`. Estimated 3–4× speedup on
Phase 2 with threading, bringing it to ~0.2s.

### Phase 3 threading: measured numbers

Synthetic benchmark matching vmmc_2d dimensions (216 chambers × 5688 pairs):

| Configuration | Time | Speedup |
|---|---|---|
| Sequential, `Dict{Vector{Int64}}` | 0.683s | 1× |
| 8 threads, `Dict{Vector{Int64}}` | 0.282s | **2.4×** |
| 8 threads, `Dict{NTuple{6,Int64}}` | ~0.27s | **~2.5×** |

The ceiling is GC contention from `copy(key_buf)` allocations in the Dict inner loop.
Threads compete for GC pauses, limiting efficiency to ~25% of theoretical 8×.

Better approach: replace the inner `Dict` with a fixed-capacity sorted array of
`(key::NTuple{N,Int64}, coeff::Rational{Int64})` pairs. With N_ATOMS=6 and at most
~8 distinct exponent classes per pair, a 16-element stack array covers all cases.
Zero heap allocation in the inner loop. Estimated threading efficiency: **6–8×**
(0.7s → 0.09–0.12s). Implementation: ~50 lines, moderate effort. The absolute saving
(~0.6s) is only significant once Phase 2 is also in Julia and the EBE total drops to ~3s.

### Theoretical best with full Julia EBE (Phase 2 + Phase 3, shared process)

| Phase | Current `-julia` | With Julia EBE + parallel τ-BFS |
|---|---|---|
| τ-BFS | 16s | ~8s (parallelised, see §5) |
| Ergodicity | 2s | 2s |
| Phase 2 BFS | 10s | ~0.7s (same Julia process) |
| Phase 3 | 0.7s | ~0.1s (threaded, no-alloc keys) |
| Export / startup | 2.5s | ~1s (two small JSON handoffs) |
| **Total vmmc\_2d** | **~33s** | **~12s** |

The unmovable floor is τ-BFS. With Mathematica parallelism on τ-BFS and Julia for
EBE, total vmmc_2d wall time reaches approximately **12s**, with τ-BFS still dominant
at ~8s.

---

## 5. Parallelism

### Mathematica τ-BFS parallelism

Each of the 56 orbit reps can be BFS'd completely independently.
`ParallelTable[$dbcBuildStateLeaves[reps[[i]], ...], {i, N}]` distributes them
across Mathematica sub-kernels.

Practical overhead: sub-kernel launch ~2–5s, `Algorithm` closure serialisation ~0.5s.
With 8 sub-kernels and 56 tasks of ~0.3s each:

```
Serial:   56 × 0.3s = 16.8s
Parallel: 5s launch + ceil(56/8) × 0.3s = 5s + 2.1s ≈ 7s
```

Practical gain: **~2× (saving ~9s)**, not 8×, because launch overhead is ~30% of
total runtime. Still the highest-return single change available, implementable as a
`-parallel` flag.

### Julia Phase 3 threading (current)

Measured 2.4× on 8 threads. GC contention is the ceiling (§4 above). Moving to
no-alloc sorted-array keys would yield **6–8×** on Phase 3 specifically.

### Julia Phase 2 threading

Within each BFS level, all k=12 neighbour sigma patterns can be sent to HiGHS in
parallel. BFS levels for vmmc_2d have width 2–8 in the middle, peaking near k=12.
Estimated practical speedup: **3–4×** on Phase 2 (0.65s → ~0.2s).

### GPU

The only GPU-parallelisable phase is Phase 3.

The 216 chambers × 5688 pairs = 1.23 M independent tasks each perform:
- 6-element integer vector arithmetic (`v_coeffs + energy_coeffs[state]`)
- Rational coefficient accumulation into at most ~8 buckets
- Bucket zero-check

With a custom CUDA kernel where each thread handles one (chamber, pair) task and
all state fits in registers (N_ATOMS=6 → 6 Int64s per key = 48 bytes), throughput
on an A100 would be approximately 0.1 ms for the entire Phase 3 — a **7000× speedup**
over the current 0.7s CPU path.

However: Phase 3 is already 0.7s. GPU would save at most 0.7s from a 33s total —
completely negligible. GPU only becomes relevant if the number of genuine chambers
scales to ~10,000+. This requires algorithms with ~15 contradictory pairs
(2^15 ≈ 32K patterns, ~10K genuine chambers after filtering) — possible for large VMMC
with many particle types and distance classes (e.g., 4×4 grid, 4+ particle types).
For all current examples, GPU investment would not pay off.

---

## 6. Robustness Concerns for New Algorithms

### 6.1 Energy linearity is silently assumed, not enforced

`$dbcCompactTerms` extracts `vCoeffs = -Coefficient[Lsym, symParamList[[j]]]`. This
assumes the Boltzmann exponent is **linear in coupling atoms**. If an algorithm uses a
non-linear energy (e.g., `E ∝ J₁ * J₂` or `E ∝ couplingJ[...]²`), the Coefficient
call returns 0 for all atoms. All terms land in the same bucket and appear to cancel,
producing a false PASS.

The guard `$dbcFS` fires when coefficients are non-integer — but if the non-linear
term has integer-valued coefficients at the specific J\* substituted, the check can
still pass accidentally. The Mathematica path handles this correctly because
`$dbcIsExpZero` groups by the actual polynomial exponent expression. **For any algorithm
with non-linear coupling energy, always verify without `-julia` first.**

### 6.2 Degenerate filter only handles pairwise contradictions

The filter identifies pairs (i,j) where `condEL[i] = -condEL[j]` (same hyperplane,
opposite sides) and removes the appropriate degenerate sigma patterns. This handles the
standard case where each hyperplane `L(J)=0` appears exactly twice in the condition list.

It would fail to handle correctly:
- Three or more conditions on the same hyperplane (e.g., additional sentinel conditions).
- Conditions whose degeneracy lies on a lower-dimensional intersection (a 1D edge
  rather than a full hyperplane). For example, `J₁ < J₂`, `J₂ < J₃`, and `J₁ = J₃`
  share a 1D corner — the current pairwise filter would not remove the corner pattern.
  In all tested examples this has not occurred, but there is no guard against it for
  novel algorithms.

### 6.3 Non-strict BFS completeness relies on unproven graph connectivity

The claim that non-strict BFS visits all genuine open chambers rests on the graph of
chambers being connected through degenerate boundary patterns. For vmmc_2d with 3
contradictory pairs this is verified empirically; there is no algebraic proof of
completeness for arbitrary hyperplane arrangements. An algorithm whose Piecewise
conditions produce an unusual arrangement could have disconnected components even under
non-strict BFS. In this case, the checker would silently verify only a subset of chambers
and report PASS without warning.

### 6.4 `DynamicSymParams` completeness is unchecked

If the user's `DynamicSymParams` function omits a coupling atom that appears in the
algorithm, that atom remains as a free symbol. `FindInstance` does not constrain it, so
jStar assigns it an arbitrary rational. Two different atom values could produce the same
rational key in Phase 3, concealing a true violation and returning a false PASS. There is
no guard that checks whether the leaf weights contain free symbols not in `allSymParams`.

### 6.5 maxDepth is a hard error, not a configurable warning

If `RunWithBitsAT` reaches maxDepth before the algorithm returns a state, the checker
exits with an error. For algorithms with probabilistic cluster sizes (Swendsen-Wang,
VMMC on larger grids), maxDepth=22 bits (the default) may be insufficient. The error
message correctly says to increase `-maxDepth`, but does not estimate the required value
— the user must determine this empirically. More importantly, there is no warning mode
that truncates and marks the result as incomplete; the only option is exit.

### 6.6 Integer overflow in v_coeffs

Phase 3 stores `v_coeffs + energy_coeffs[state]` as `Vector{Int64}` keys. For vmmc_2d
with energy coefficients in [-3, 3] and v_coeffs in [-2, 2], overflow is impossible.
For algorithms with many particle types, high-order interactions, or large geometric
prefactors, energy coefficients could approach `Int64` limits. A silent wraparound would
produce incorrect grouping and either false PASSes or spurious FAILs. No range check
exists in the export or Julia code.

### 6.7 `seqBernoulli` requires rational or symbolically reducible weights

`RandomChoice[weights -> elements]` uses `seqBernoulli`, which calls `makeRealVar[]`
for each weight comparison. This works when weights are rational or symbolic expressions
that Mathematica can compare against a threshold. It breaks when weights contain
irrational constants (e.g., `Exp[-1.5]` or `π/4`) that cannot be resolved to True/False
by Mathematica. The fallback to `$dbc$cantHandle` catches this, but the resulting error
message ("RandomReal[] was used in an unsupported way") is not informative about the
root cause.

### 6.8 NormalDistribution truncation is silent above 1%

`RandomVariate[NormalDistribution[μ,σ]]` is discretised to `Floor[nGrid/2]` integers
around μ. A warning is printed when more than 1% of probability mass is discarded. Below
that threshold, no warning is given and the check proceeds on the truncated distribution.
For σ approaching `nGrid/2`, the discarded tail can be non-negligible without triggering
the warning. The PASS/FAIL result certifies the discretised distribution, not the
original continuous one — this distinction is not surfaced in the summary output.

---

## 7. Recommended Optimisations (Ordered by Return)

| Optimisation | Effort | Saving (vmmc\_2d) | Notes |
|---|---|---|---|
| Julia Phase 2 (HiGHS, same process as Phase 3) | Medium | ~9s | Two small JSON round-trips; 3.4× without sysimage, 10× with |
| Mathematica τ-BFS `-parallel` flag | Low–Medium | ~8s | `ParallelTable` over orbit reps; 2× practical gain |
| PackageCompiler.jl sysimage (HiGHS + JSON3) | Low | ~2s | Eliminates JIT overhead; one-time setup |
| Julia Phase 3 no-alloc keys + threading | Medium | ~0.6s | Only worthwhile after Phase 2 is in Julia |
| GPU Phase 3 (CUDA.jl) | High | ~0.7s | Only worthwhile at ~10K+ chambers |

With the first two items implemented, total vmmc\_2d wall time drops from ~33s to
approximately **12s**, with τ-BFS the remaining dominant cost at ~8s.
