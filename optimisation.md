# Julia Acceleration: Diagnosis, Fixes and Architecture

## Verified benchmark results (measured runs)

| Run | τ-BFS | Ergodicity | EBE | Total wall time |
|-----|-------|------------|-----|-----------------|
| `vmmc_2d` (no flags) | 16s | 2s | 205s | 226s |
| `vmmc_2d -julia` (final) | 16s | 2s | 15s | ~33s |
| `single_metropolis -julia` | 4s | 1s | 15s | ~20s |
| **vmmc_2d speedup (EBE only)** | — | — | **14×** | **7×** |

`vmmc_2d -julia` breakdown:
- Mathematica Phase 2 BFS (non-strict): ~10s (2048 patterns visited, 512 feasible, 216 genuine)
- Compact export to JSON: 1.2s (5688 pairs, 4368 leaves, 217 unique weights)
- Julia startup + load: ~1.4s
- Julia Phase 3: 0.7s (216 chambers × 5688 pairs)
- Total EBE: ~15s

---

## Architecture (final correct state)

```
Mathematica:
  Phase 1 — Extract allConds from leaf weights (k=12 for vmmc_2d)
  Phase 2 — NON-STRICT FindInstance BFS:
               - Uses !allConds[[j]] (non-strict) for False conditions
               - Bridges all octants of the hyperplane arrangement
               - Finds 512 feasible patterns (2048 visited)
             DEGENERATE FILTER — removes boundary patterns:
               - Strict pairs (Less/Greater): filter when both False
               - Non-strict pairs (LessEqual/GreaterEqual): filter when both True
               - Mixed pairs: no filter needed (no shared boundary)
               → 216 genuine open chambers for vmmc_2d
  Export  — Compact JSON: unique_weights (sigma-substitution) + feasible_sigmas

Julia (ebe.jl):
  Phase 3 only — for each sigma in feasible_sigmas:
    1. Pre-compute active Case for each of 217 unique weights
    2. For each pair: accumulate T(i→j) and T(j→i) terms
    3. Group by integer key = v_coeffs + energy_coeffs[state]
    4. Check each group sums to zero (Rational{Int64}, no BigInt)
```

---

## Bug 1: Compact export failure — `PiecewiseExpand` produces unmatched conditions

**Symptom**: `Julia compact export: unrecognized leaf weight — falling back to Mathematica`

**Root cause**: Old `$dbcCompactWeight` called `PiecewiseExpand` which creates complex OR/AND
conditions (e.g., `J12_1 < J12_2 || J13_1 >= J13_2`) not present in allConds.

**Fix**: Sigma-substitution — for each unique weight, substitute True/False directly for each
active condition's 2^k' possible assignments. Mathematica auto-evaluates
`Piecewise[{{val, True}}, 0] → val`. 217 unique weights, 1105 cases, 0 failures.

---

## Bug 2: Phase 3 slow — `eval_v` with `Rational{BigInt}`

**Symptom**: Phase 3 taking ~60s in early PoC.

**Root cause**: ~2.74M calls to `eval_v(term.v_coeffs, jStar)` (22μs each = ~60s).

**Fix**: Three changes:
1. Pre-compute active terms per weight per chamber (27 × 217 = 5,859 lookups vs 2.74M)
2. Group by integer coefficient vector (`v_coeffs + energy_coeffs[state]`), not by
   `Rational{BigInt}` dot-product — algebraically exact, no BigInt needed
3. Allocation-free `get_terms` using element-wise comparison

**Result**: Phase 3 time 0.7s for 216 chambers (down from ~60s PoC).

---

## Bug 3: CDDLib Phase 2 finds wrong chambers

**Root cause**: CDDLib's non-strict `>= 0` for True constraints includes degenerate boundary
chambers where contradictory conditions are simultaneously "True". Also, CDDLib's BFS has a
connectivity problem: chambers on opposite sides of a contradictory-pair boundary are separated
by a 2-bit flip, not reachable via CDDLib's 1-bit BFS.

**Fix**: Abandon CDDLib. Use Mathematica's FindInstance BFS.

---

## Bug 4: Strict BFS in Phase 2 misses 7/8 of parameter space (critical correctness bug)

**Symptom**: 27 chambers found instead of 216 genuine ones. Julia reports PASS, but only
1 of 8 octants of parameter space is verified.

**Root cause**: `FindInstance[True, allSymParams, Rationals, 1]` returns a specific interior
point, not the all-zeros boundary. The initial sigma is inside one genuine open chamber. From
there, strict BFS cannot cross any contradictory-pair boundary (both 1-bit paths lead to
infeasible patterns). vmmc_2d has 3 such pairs → 8 disconnected components → 7/8 unverified.

Confirmed empirically: initial jStar = `{J₁₂₁ → -9/5, J₁₂₂ → 3/5, ...}`. At this point:
J₁₂₁ < J₁₂₂ is True (so J₁₂₁ > J₁₂₂ cannot be reached by strict 1-bit BFS).

**Fix**: Use NON-STRICT BFS. The non-strict BFS can hop through degenerate boundaries
(where condEffLhs[i]·J = 0) to reach all octants. Starting from a generic interior point,
the BFS first hops to the degenerate boundary (all-zeros: J₁₂₁ = J₁₂₂), then to the other
octant. All 8 octants are reachable.

**Result**: 512 feasible patterns found → degenerate filter → 216 genuine chambers.

---

## Bug 5: Degenerate filter misses non-strict (LessEqual/GreaterEqual) boundary patterns

**Symptom**: `single_metropolis -julia` reported FAIL (false violations). The algorithm is
correct and should PASS.

**Root cause**: `single_metropolis` uses `dE <= 0` (LessEqual) in Piecewise conditions.
For a LessEqual/GreaterEqual pair, the degenerate boundary occurs when BOTH conditions are
True (L ≤ 0 AND L ≥ 0 → L = 0), not when both are False.

The original filter only checked "both False" (appropriate for Less/Greater pairs). For
LessEqual/GreaterEqual pairs, "both True" sigma patterns correspond to the L=0 boundary —
a measure-zero set where Julia's integer-key Phase 3 correctly identifies non-zero DB
residuals (since it checks for ALL J, not just at L=0).

**Fix**: Classify pairs by strictness and apply the appropriate filter:
- Strict pairs (Less/Greater): filter when `sigma[i]=0 AND sigma[j]=0`
- Non-strict pairs (LessEqual/GreaterEqual): filter when `sigma[i]=1 AND sigma[j]=1`
- Mixed pairs (Less+GreaterEqual, etc.): no degenerate boundary, no filter needed

**Result**: `single_metropolis -julia` finds 48 genuine chambers (from 132 BFS-feasible),
removes 84 degenerate boundary patterns, correctly reports PASS.

---

## Why the old Mathematica non-strict path gave correct results despite degenerate chambers

Mathematica's own Phase 3 evaluates the DB expression at a specific jStar point. For
degenerate boundary patterns:
- Strict pairs: jStar is at the L=0 boundary. At this point, all coupling constants
  involved in L appear equal, making the DB expression coincidentally zero.
- Non-strict pairs: jStar is at the L=0 boundary. At L=0, E(s) = E(t) for affected pairs,
  making exp(-βE(s)) = exp(-βE(t)) and the DB condition trivially satisfied.

Julia's integer-key Phase 3 is MORE rigorous: it checks if the DB expression is zero for
ALL J, not just at a specific boundary point. This is why degenerate sigmas that pass
Mathematica's check can fail Julia's — and why they must be filtered before Phase 3.

---

## Correctness of the complete implementation

The non-strict BFS + degenerate filter correctly identifies all genuine open chambers
because:

1. **Completeness**: Non-strict BFS bridges all octants of the hyperplane arrangement via
   degenerate boundary patterns (both-conditions-boundary). No octant is missed.

2. **No false inclusions**: The degenerate filter removes all boundary patterns:
   - Strict pairs: both-False patterns (on the gap between two strict half-spaces)
   - Non-strict pairs: both-True patterns (on the shared boundary of two non-strict half-spaces)

3. **No false exclusions**: Genuine open chambers have exactly one True condition for each
   strict contradictory pair, and at most one True condition for each non-strict pair at the
   boundary. Interior chambers are never accidentally filtered.

4. **DB sufficiency**: If DB holds in all genuine open chambers, it holds everywhere by
   continuity of the transition probabilities as functions of coupling constants. Boundary
   hyperplanes are measure-zero and don't need separate verification.

---

## Files changed

**`dbc_core.wl`:**
- `$dbcCompactWeight` — sigma-substitution instead of PiecewiseExpand
- `$dbcExportCompact` — exports `feasible_sigmas` (list of genuine sigma vectors)
- `$dbcEBECheckLeavesJulia` — non-strict Phase 2 BFS + degenerate filter (both strict and
  non-strict pair types), exports filtered sigmas, calls Julia for Phase 3 only

**`julia_poc/ebe.jl`** — Phase 3 only:
- `UniqueWeight` struct: `active_cond_idxs` + `cases` (sigma-substitution)
- `get_terms`: allocation-free element-wise comparison
- `run_phase3`: integer coefficient vector keys, pre-computed active terms per chamber

**`julia_poc/Project.toml`** — JSON3, LinearAlgebra only (CDDLib, Polyhedra, HiGHS removed)

---

## Performance profile (measured)

| Phase | Mathematica (no -julia) | Julia (-julia, final) |
|---|---|---|
| τ-BFS | ~16s | ~16s |
| Ergodicity | ~2s | ~2s |
| Phase 2 BFS (Mathematica) | ~10s (512 feasible from 2048) | ~10s (same) |
| Degenerate filter | — | <0.1s |
| Phase 3 DB check | ~195s (Mathematica symbolic, 512 chambers) | ~0.7s (Julia integer-key, 216 genuine) |
| Export + Julia startup | — | ~2.5s |
| **Total** | **~226s** | **~33s** |

The bottleneck in both paths is now Phase 2 (~10s) + τ-BFS (~16s) = ~26s Mathematica floor.
Phase 3 speedup: ~195s → 0.7s = ~280× for vmmc_2d.

## Note on Julia Phase 2 (investigated but not deployed)

Julia could do Phase 2 using HiGHS LP feasibility checks (~1.4s for 2048 LPs). However,
HiGHS JIT compilation costs ~10s on Julia cold start (each `-julia` invocation starts a new
process), making total EBE ~5s longer than the Mathematica Phase 2 approach.

With a Julia sysimage (PackageCompiler.jl) containing HiGHS precompiled to native code,
Julia Phase 2 would be ~1.4s vs Mathematica's ~10s — saving ~9s per run. This is not
implemented due to the complexity of sysimage management, but would be the natural next
optimization step.
