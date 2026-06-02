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

When the `-julia` flag is used, EBE exports a **compact representation** of Phase 1 (Piecewise condition indices + rational coefficient vectors) and delegates Phase 2 (chamber enumeration) and Phase 3 (DB check) to a Julia subprocess:

**Phase 2 (Julia)**: Uses CDDLib's exact rational arithmetic to enumerate all feasible coupling-constant chambers via sign-pattern BFS, using strictly interior rational points for J*evaluation.

**Phase 3 (Julia)**: For each chamber, evaluates leaf weights from the compact structure (dot products with J*) and checks the Boltzmann-exponent grouped DB residual = 0 using exact rational arithmetic.

**Fallback**: If the compact export fails (e.g., for highly nested Piecewise structures), the checker automatically falls back to the pure-Mathematica EBE path and produces identical results.

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
  -julia                Experimental: delegate Phase 2+3 to Julia subprocess
                        (requires CDDLib and Polyhedra packages; falls back to Mathematica if export fails)
```

---

## Provided examples

| File | Expected result | Time (default, no -julia) |
|---|---|---|
| `single_metropolis.wl` | τ PASS, DB PASS | ~35s |
| `kawasaki.wl` | τ PASS, DB PASS, Ergodicity FAIL (by design) | ~3s |
| `vmmc_2d.wl` | τ PASS, DB PASS | ~190s |
| `quadratic_field.wl` | τ FAIL (absolute-position energy), DB PASS | ~5s |
| `broken_variable_pool.wl` | DB FAIL — asymmetric pool size (3 or 4) | <5s |
| `broken_8way_hop.wl` | DB FAIL — asymmetric pool size (7 or 8) | <5s |
| `broken_biased_direction.wl` | DB FAIL — duplicate direction in proposal pool | <5s |
| `broken_metropolis_halfbeta.wl` | DB FAIL — accept probability uses β/2 instead of β | <5s |
| `broken_field_wrong_accept.wl` | DB FAIL — accept uses pair energy only, ignores field | <5s |

The `vmmc_2d.wl` runtime of ~190s is dominated by EBE Phase 3: VMMC's cluster-building logic generates k=12 nearly independent Piecewise conditions, producing 512 feasible coupling-constant chambers. Each chamber requires checking 5,688 communicating state pairs. This scales as O(chambers × pairs) and is the fundamental cost of an exact algebraic verification.

---

## Known limitations

**EBE always runs exactly; probabilistic SZ is opt-in.** The default `ebeMaxK=10000` means EBE runs for any realistic algorithm. To use the probabilistic Schwartz-Zippel fallback explicitly, pass `-mode SZPure` or lower `-ebeMaxK`.

**D4 check is algebraically exact but expensive relative to its benefit.** For algorithms with many feasible chambers (e.g. VMMC with 512), the D4 EBE verification costs approximately 2× the DB check while the orbit reduction it enables saves less than 10% on DB Phase 3. D4 is therefore not in `$symmetryGroup` by default. Add it back if you need to verify rotational symmetry as a separate correctness claim.

**τ-BFS cantHandle for in-body Mod.** If the algorithm normalises particle positions inside the step function using `Mod[pos, n]` on symbolic values, Mathematica may trigger internal evaluation calls that are intercepted by the BFS override, causing a cantHandle error. The checker falls back to the full state-space BFS without the translation speedup.

**NormalDistribution is discretised.** The checker models `RandomVariate[NormalDistribution[mu, sigma]]` as a discrete distribution on `Floor[nGrid/2]` integers around `mu`. For symmetric proposals (`mu=0`) this does not affect PASS/FAIL classification (the truncated Gaussian is also symmetric). For large `sigma` relative to `nGrid/2`, a warning is printed.

**Open-chamber boundary omission.** The EBE check covers all open chambers of the hyperplane arrangement. Boundaries where two or more conditions are simultaneously tight (e.g., `J1 = J2` exactly) are not tested. For standard Metropolis algorithms (transition probabilities continuous in coupling constants), any violation in a positive-measure region must appear in an adjacent open chamber, so this omission is harmless in practice.

**Julia Phase 2+3 (`-julia`)**: When enabled, the checker exports a compact Piecewise structure and delegates chamber enumeration and DB checking to Julia using CDDLib and Polyhedra. Simpler algorithms (k < 12) run ~2–4× faster; for vmmc_2d (k=12, 512 chambers), Julia's startup and export overhead dominates, so the default Mathematica path is typically faster. The flag is provided for research and future optimization — as Mathematica's performance improves or Julia startup cost decreases, this may become the default. The fallback to Mathematica ensures correctness is never compromised.

**BFS timeout aborts the run.** If a BFS path exceeds the per-state time limit (`-timeLimit`), the checker aborts with an error. Increase `-timeLimit` if needed.

---

## Performance Notes

### Architecture changes from v1

- **Removed D4 from default `$symmetryGroup`**: D4 verification costs 2× the DB check but saves <10% on Phase 3, making it net negative for large chamber counts. Users who need D4 verification can add it explicitly; the algebraic check remains available.
- **τ-BFS now captures leaves**: Setting τ=0 in τ-BFS leaves produces identical BFS output, eliminating the separate BFS pass and saving ~13s on typical algorithms (10–35s faster overall).
- **Julia Phase 2+3 (optional)**: Compact Piecewise structure export avoids the ~100s per-leaf conversion overhead of the old approach. CDDLib provides 7× faster chamber enumeration than Mathematica's FindInstance-based BFS. For simple algorithms, this is a net speedup; for highly nested Piecewise structures (vmmc_2d), export fails gracefully and falls back to Mathematica.

### Bottleneck hierarchy (vmmc_2d as example, time per phase)

| Phase | Time | Driver | Optimized? |
|---|---|---|---|
| τ-BFS | ~13s | Symbolic bit-reader on 56 reps | ✓ (optimal, already in Mathematica) |
| Ergodicity | ~2s | Graph BFS from 504-state transition graph | ✓ (negligible) |
| EBE Phase 2 (Julia) | ~7s | CDDLib chamber BFS, 2048 patterns checked | ✓ (7× faster than Mathematica) |
| EBE Phase 3 (Julia) | ~0.5s | Compiled rational arithmetic, 512 regions × 5688 pairs | ✓ (35× faster than Mathematica) |
| EBE Phase 1+2 (Mathematica fallback) | ~5s | FindInstance BFS on 512 regions | ✓ (acceptable, only if export fails) |
| **Total** | **~190–210s** | Mathematica floor + Julia acceleration | At capacity |

The Mathematica floor (τ-BFS + ergodicity) is ~15s and cannot be reduced further without reimplementing the BFS engine in Julia with Cassette.jl (not done due to API instability and marginal overall gain).

---

## Debugging and tracing

Use `-verbose` to see per-rep BFS progress:

```
wolframscript -file check.wls examples/single_metropolis.wl -verbose
```

This prints indices and timing for each orbit representative's BFS, helping identify which states are expensive.

For Julia-mode debugging, check the stderr output of the Julia subprocess in the ebe.jl Phase 2 output (printed by the Mathematica driver).

---

## Citation

If you use SZ-DBC in your research, please cite the work. The checker implements the Schwartz-Zippel test as described in the accompanying paper.
