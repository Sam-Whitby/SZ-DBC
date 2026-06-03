# SZ-DBC: Detailed Balance Checker for Lattice MCMC

SZ-DBC verifies that a lattice MCMC algorithm satisfies **detailed balance** — the condition required to correctly sample the Boltzmann distribution. It intercepts all random-number calls made by the algorithm, reconstructs the exact symbolic transition probabilities via BFS, and checks the DB condition algebraically over all coupling-parameter regions.

---

## Quick start

```
wolframscript -file check.wls examples/single_metropolis.wl
```

The algorithm file must define: `$nGrid`, `$particleTypes`, `$seedState`, `Algorithm[state_]`, `energy[state_]`, and `DynamicSymParams[states_]`. See the examples directory for complete templates.

---

## What is detailed balance?

An MCMC algorithm moves from state `s` to `t` with probability `T(s→t)`. It correctly samples the Boltzmann distribution `π(s) ∝ exp(-β E(s))` if and only if:

```
T(s→t) · π(s) = T(t→s) · π(t)    for all pairs (s, t)
```

---

## How it works

### State format

Every state is a sorted list of `{position, type}` pairs:

```mathematica
{{{1,1}, 1}, {{1,2}, 2}, {{2,3}, 3}}
```

Positions are integers 1–`nGrid`; types are positive integers. The checker applies periodic boundary conditions automatically.

### BFS engine

The core of the checker intercepts every call to `RandomReal`, `RandomInteger`, `RandomChoice`, and related functions inside `Algorithm`. Each call is replaced by a **symbolic bit-reader**: one bit is consumed per binary decision, and the current path weight is updated multiplicatively. This produces a complete list of `(nextState, pathWeight)` pairs — one per leaf of the binary decision tree — with weights summing to 1.

For `RandomChoice[list]` with `n = Length[list]`:
- Reads `k = ceil(log2(n))` bits to produce an index `0..n-1`.
- Leaves with `index >= n` are discarded (rejection sampling).
- The path weight is multiplied by `2^k / n` so that each valid leaf carries weight exactly `1/n`.

For weighted `RandomChoice[weights -> elements]`, a sequential Bernoulli decomposition is used (`seqBernoulli`), which is exact.

### Step 1 — State enumeration

All `N! / (n₁! n₂! ...)` states for the given particle type multiset are enumerated on the `nGrid × nGrid` torus. The count is verified against the theoretical combinatorial formula.

### Step 2 — Orbit computation

Translation-only orbits are computed first. D4 rotational symmetry can be verified algebraically in a later step if declared.

### Step 4 — τ-BFS: translational invariance check and leaf capture

Particle positions are augmented with a symbolic offset `{τr, τc}` and the BFS engine runs on these τ-augmented states. If the algorithm uses only **pairwise differences** for all spatial decisions, τ cancels algebraically from every leaf weight and the check passes. This is verified symbolically — no numerical sampling.

When translation invariance passes, the τ-leaves are immediately converted to concrete BFS leaves by substituting `τr→0, τc→0`. This avoids a redundant second BFS pass: setting τ=0 in a τ-BFS leaf gives exactly the same `{bits, nextState, weight}` triple that a separate BFS on the non-augmented representative would produce.

If the τ-BFS encounters a `cantHandle` condition (e.g., `Mod[position, n]` on symbolic values), it falls back to a standard BFS on the original (non-τ-augmented) states.

### Step 5b — Rotational invariance (D4, algebraic via EBE) — optional

If `"D4"` is in `$symmetryGroup` and translation invariance passed, D4 is verified **algebraically** using EBE. This checks `T(s→t) = T(R·s→R·t)` exactly within every feasible parameter region, for two D4 generators (rotation by 90° and left-right reflection), which by group theory is sufficient for all 8 D4 elements.

**D4 is not included in `$symmetryGroup` by default** because for algorithms with many coupling-constant chambers (e.g. VMMC with 512 chambers), the D4 verification costs as much as two full DB checks while providing negligible speedup on the DB check itself. Include `"D4"` only if you specifically need to verify rotational symmetry of the algorithm.

### Step 7 — Ergodicity

A graph BFS is run from `$seedState` over the transition graph derived from the orbit-rep BFS leaves, checking whether every state is reachable.

### Step 8 — Detailed balance: Exhaustive Branch Enumeration (EBE)

The leaf weights are symbolic expressions in the coupling parameters. The Piecewise branch conditions define a **hyperplane arrangement** in parameter space. EBE enumerates all feasible chambers via a BFS starting from a rational interior point. Within each chamber the Piecewise conditions resolve to constants, and the DB equation reduces to a check of the form `sum(rational × exp(-β × rational)) = 0`, verified exactly by grouping exponent classes and checking rational coefficients.

#### Julia acceleration (optional)

When the `-julia` flag is used, EBE uses a two-phase approach:

**Phase 2 (Mathematica)**: A non-strict BFS enumerates all sigma patterns (coupling-parameter sign patterns) where the Piecewise conditions have a definite True/False assignment. The non-strict negation (`>=` for False conditions) bridges disconnected octants of the hyperplane arrangement, ensuring all genuine open chambers are found. After BFS, degenerate boundary patterns — those lying on a measure-zero hyperplane where a pair of contradictory conditions are simultaneously at their shared boundary — are filtered out algebraically. For strict conditions (Less/Greater), degeneracy occurs when both are False; for non-strict conditions (LessEqual/GreaterEqual), when both are True.

**Phase 3 (Julia)**: For each genuine chamber, evaluates leaf weights from the compact structure using sigma-substitution (True/False directly into Piecewise expressions, no PiecewiseExpand), groups DB terms by integer exponent coefficient vector, and checks that each group sums to zero using exact `Rational{Int64}` arithmetic. No BigInt and no specific J* evaluation point needed — grouping by integer coefficient vector is algebraically exact and equivalent to Mathematica's symbolic exponent grouping.

**Fallback**: If the compact export fails (e.g., leaf weights contain non-standard β factors), the checker automatically falls back to the pure-Mathematica EBE path.

---

## Algorithm file format

```mathematica
$nGrid         = 3;           (* lattice side length *)
numBeta        = 1.0;         (* inverse temperature for numerical mode *)
$maxD2         = 2;           (* max squared distance for energy terms *)
$particleTypes = {1, 2, 3};   (* type multiset *)

$seedState = Module[...];     (* canonical starting state *)

$symmetryGroup = {"translation"};   (* "D4" may be added to also verify D4 symmetry *)

energy[state_] := ...         (* bare energy, no beta factor *)

Algorithm[state_List] :=      (* one MCMC step; may call RandomReal/RandomChoice/etc. *)
  Module[..., ...]

DynamicSymParams[states_List] :=   (* returns coupling parameter atoms *)
  <|"couplings" -> {...}, "numericParams" -> {}|>
```

### Supported random primitives

| Primitive | Behaviour |
|---|---|
| `RandomReal[]` | Uniform on [0,1]; used for Metropolis acceptance |
| `RandomReal[{lo,hi}]` | Uniform continuous token |
| `RandomInteger[{lo,hi}]` | Exact rejection sampling; weight corrected |
| `RandomChoice[list]` | Exact rejection sampling; weight corrected |
| `RandomInteger[{lo,hi}, count]` | count independent integers via rejection sampling |
| `RandomInteger[n, count]` | count independent integers in {0..n} |
| `RandomChoice[weights->elems]` | Sequential Bernoulli decomposition |
| `RandomVariate[NormalDistribution[mu,sigma]]` | Discretised to integers in `[Round(mu)-nMax, Round(mu)+nMax]` where `nMax = Floor[nGrid/2]`. Weights are CDF-based and renormalised. A warning is printed when >1% of mass is discarded. For symmetric proposals (`mu=0`) this discretisation does not affect PASS/FAIL classification. |
| `RandomPermutation`, `RandomSample` | Fisher-Yates via `RandomInteger` |

**Note on inverted bounds:** `RandomInteger[{hi, lo}]` with `hi < lo` is now detected as a `cantHandle` error and reported as a likely algorithm bug, rather than silently discarding all paths.

### Coupling atoms

Coupling atoms (`couplingJ[type1, type2, dist]` or `couplingJ[type1, type2]`) are automatically canonicalised so that `couplingJ[b,a,...]` with `b > a` is rewritten to `couplingJ[a,b,...]`. Both 2-argument and 3-argument forms are supported.

### Requirements for τ-BFS to work

All spatial operations must use **pairwise differences** (not absolute positions) for comparisons and energy calculations:

```mathematica
(* correct: difference-based occupancy check *)
Mod[p[[1]] - newPos[[1]], n] === 0

(* wrong: absolute-position Mod inside the algorithm body *)
newPos = Mod[rawPos - 1, n] + 1   (* breaks τ-BFS *)
```

---

## Command-line options

```
wolframscript -file check.wls <algorithm.wl> [options]

  -mode Symbolic        (default) EBE exact DB check
  -mode SZPure          Probabilistic Schwartz-Zippel check
  -mode FullSimplify    FullSimplify each DB expression
  -mode Numerical       Numerical MCMC comparison only
  -szRepeats N          SZ random evaluation points (default 30)
  -ebeMaxK N            Max Piecewise conditions for EBE (default 10000 — effectively unlimited)
  -maxDepth N           BFS bit depth limit (default 22)
  -timeLimit T          Per-state time limit in seconds (default 120)
  -verbose              Print per-rep BFS progress
  -julia                Delegate Phase 3 to Julia subprocess (faster for algorithms
                        with many Piecewise conditions); falls back to Mathematica if export fails
```

---

## Provided examples

| File | Expected result | Time (no -julia) | Time (-julia) |
|---|---|---|---|
| `single_metropolis.wl` | τ PASS, DB PASS | ~35s | ~20s |
| `kawasaki.wl` | τ PASS, DB PASS, Ergodicity FAIL (by design) | ~3s | ~3s |
| `vmmc_2d.wl` | τ PASS, DB PASS | ~226s | ~33s |
| `quadratic_field.wl` | τ FAIL (absolute-position energy), DB PASS | ~5s | ~5s |
| `broken_variable_pool.wl` | DB FAIL — asymmetric pool size (3 or 4) | <5s | <5s |
| `broken_8way_hop.wl` | DB FAIL — asymmetric pool size (7 or 8) | <5s | <5s |
| `broken_biased_direction.wl` | DB FAIL — duplicate direction in proposal pool | <5s | ~15s |
| `broken_metropolis_halfbeta.wl` | DB FAIL — accept probability uses β/2 instead of β | <5s | <5s (fallback) |
| `broken_field_wrong_accept.wl` | DB FAIL — accept uses pair energy only, ignores field | <5s | <5s |

The `vmmc_2d.wl` runtime is dominated by EBE Phase 3 in the default path: VMMC's cluster-building logic generates k=12 Piecewise conditions, producing 216 genuine open chambers (the non-strict Phase 2 BFS visits 2048 sign patterns and finds 512, of which 296 are degenerate boundary patterns that are filtered). Each of the 216 chambers requires checking 5,688 communicating state pairs. The `-julia` flag reduces EBE time from ~205s to ~15s via Julia's integer-coefficient-vector DB check (7× overall speedup on this example).

---

## Known limitations

**EBE always runs exactly; probabilistic SZ is opt-in.** The default `ebeMaxK=10000` means EBE runs for any realistic algorithm. To use the probabilistic Schwartz-Zippel fallback explicitly, pass `-mode SZPure` or lower `-ebeMaxK`.

**D4 check is algebraically exact but expensive relative to its benefit.** For algorithms with many feasible chambers (e.g. VMMC with 512), the D4 EBE verification costs approximately 2× the DB check while the orbit reduction it enables saves less than 10% on DB Phase 3. D4 is therefore not in `$symmetryGroup` by default. Add it back if you need to verify rotational symmetry as a separate correctness claim.

**τ-BFS cantHandle for in-body Mod.** If the algorithm normalises particle positions inside the step function using `Mod[pos, n]` on symbolic values, Mathematica may trigger internal evaluation calls that are intercepted by the BFS override, causing a cantHandle error. The checker falls back to the full state-space BFS without the translation speedup.

**NormalDistribution is discretised.** The checker models `RandomVariate[NormalDistribution[mu, sigma]]` as a discrete distribution on `Floor[nGrid/2]` integers around `mu`. For symmetric proposals (`mu=0`) this does not affect PASS/FAIL classification (the truncated Gaussian is also symmetric). For large `sigma` relative to `nGrid/2`, a warning is printed.

**Open-chamber boundary omission.** The EBE check covers all open chambers of the hyperplane arrangement. Boundaries where two or more conditions are simultaneously tight (e.g., `J1 = J2` exactly) are not tested. For standard Metropolis algorithms (transition probabilities continuous in coupling constants), any violation in a positive-measure region must appear in an adjacent open chamber, so this omission is harmless in practice.

**`-julia` flag**: When enabled, Mathematica runs Phase 2 (non-strict BFS + degenerate filter, same FindInstance cost as the default path), then delegates Phase 3 to a Julia subprocess. Julia groups DB terms by integer exponent coefficient vector using exact `Rational{Int64}` arithmetic — 200–300× faster than Mathematica's symbolic Phase 3 for algorithms with many chambers. For vmmc_2d (216 genuine chambers), `-julia` reduces total runtime from ~226s to ~33s. Fallback to Mathematica if the compact export fails (e.g., non-standard β factors in weights).

**BFS timeout aborts the run.** If a BFS path exceeds the per-state time limit (`-timeLimit`), the checker aborts with an error. Increase `-timeLimit` if needed.

---

## Performance Notes

### Architecture changes from v1

- **Removed D4 from default `$symmetryGroup`**: D4 verification costs 2× the DB check but saves <10% on Phase 3, making it net negative for large chamber counts. Users who need D4 verification can add it explicitly; the algebraic check remains available.
- **τ-BFS now captures leaves**: Setting τ=0 in τ-BFS leaves produces identical BFS output, eliminating the separate BFS pass and saving ~13s on typical algorithms (10–35s faster overall).
- **Julia Phase 3 (optional, `-julia`)**: Phase 3 (the expensive per-chamber DB check) is delegated to Julia. Phase 2 remains in Mathematica using a non-strict BFS that correctly bridges all octants of the hyperplane arrangement. The compact structure uses sigma-substitution (True/False directly into Piecewise) rather than PiecewiseExpand, avoiding unmatched condition expressions. Julia groups DB terms by integer coefficient vector (algebraically exact, no BigInt) for a 200–300× speedup on Phase 3.

### Correctness of Phase 2: non-strict BFS + degenerate filter

The Phase 2 BFS uses non-strict negation (`>=` for False conditions) so that degenerate boundary patterns — lying on measure-zero hyperplane boundaries — are visited and serve as bridges between disconnected octants of the arrangement. After BFS, these patterns are filtered out algebraically:

- **Strict pairs** (Less/Greater conditions): degenerate when both conditions are False (boundary `expr = 0` between `expr < 0` and `expr > 0`)
- **Non-strict pairs** (LessEqual/GreaterEqual conditions): degenerate when both conditions are True (boundary `expr = 0` shared by `expr ≤ 0` and `expr ≥ 0`)

Using strict negation (the v1 approach) for Phase 2 BFS disconnects the graph at contradictory-pair boundaries, causing the BFS to miss 7/8 of the coupling-parameter space for algorithms with 3 such pairs. The non-strict approach is complete: for vmmc_2d it finds all 216 genuine open chambers (versus 512 BFS-feasible patterns minus 296 degenerate ones).

### Bottleneck hierarchy (vmmc_2d as example, time per phase)

| Phase | Time (no -julia) | Time (-julia) | Driver |
|---|---|---|---|
| τ-BFS | ~16s | ~16s | Symbolic bit-reader on 56 reps |
| Ergodicity | ~2s | ~2s | Graph BFS on 504-state transition graph |
| EBE Phase 2 (Mathematica BFS) | ~10s | ~10s | FindInstance BFS, 2048 patterns → 216 genuine |
| EBE Phase 3 | ~195s | ~0.7s | Mathematica symbolic vs Julia integer-key |
| Export + Julia startup | — | ~2.5s | JSON write + Julia JIT (JSON3 only) |
| **Total** | **~226s** | **~33s** | **7× speedup** |

The Mathematica floor (τ-BFS + ergodicity + Phase 2) is ~28s. Phase 3 is the only phase where Julia provides a speedup.

---

## Debugging and tracing

Use `-verbose` to see per-rep BFS progress:

```
wolframscript -file check.wls examples/single_metropolis.wl -verbose
```

This prints indices and timing for each orbit representative's BFS, helping identify which states are expensive.

For Julia-mode debugging, check the stderr output printed by Mathematica: `EBE: N feasible region(s)` shows the BFS result, `EBE: M genuine open chamber(s)` shows the count after degenerate filtering, and `Julia compact export:` shows the export timing. Julia's own stderr (Load/Phase 3 timing) is printed inline.

---

## Citation

If you use SZ-DBC in your research, please cite the work. The checker implements the Schwartz-Zippel test as described in the accompanying paper.
