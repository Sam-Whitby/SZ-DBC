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

Translation-only orbits are computed first. D4 rotational symmetry is verified algebraically in a later step before being used to reduce the orbit set further.

### Step 4 — Translational invariance (algebraic τ-BFS)

Particle positions are augmented with a symbolic offset `{τr, τc}` and the BFS is re-run on these τ-augmented states. If the algorithm uses only **pairwise differences** for all spatial decisions, τ cancels algebraically from every leaf weight, and the check passes. This is verified symbolically — no numerical sampling.

If the τ-BFS encounters a `cantHandle` condition (e.g., `Mod[position, n]` on symbolic values), it prints an informative message and falls back to the non-translation path.

### Step 5 — BFS from orbit representatives

The BFS engine runs on each orbit representative. Leaves accumulate as `{bits, nextState, weight}` triples.

### Step 5b — Rotational invariance (D4, algebraic via EBE)

If `"D4"` is in `$symmetryGroup` and translation invariance passed, D4 is verified **algebraically** using EBE (Exhaustive Branch Enumeration). This runs the same Phase 1 (condition extraction) and Phase 2 (chamber enumeration) as the DB check, then in Phase 3 verifies `T(s→t) = T(R·s→R·t)` exactly within every feasible parameter region using `$dbcIsExpZero`.

Only two D4 generators are checked: rotation by 90° (`rotIdx=1`) and a left-right reflection (`rotIdx=4`). By group theory, if the condition holds for both generators it holds for all 8 D4 elements.

The EBE chambers computed here are reused by the DB check (Step 8), so no duplicate Phase 2 enumeration occurs. If the new BFS produces conditions not seen in the D4 EBE, the SubsetQ assertion detects this and falls back to running Phase 2 from scratch.

Falls back to probabilistic SZ when `k > ebeMaxK`.

### Step 7 — Ergodicity

A graph BFS is run from `$seedState` over the transition graph derived from the orbit-rep BFS leaves, checking whether every state is reachable.

### Step 8 — Detailed balance: Exhaustive Branch Enumeration (EBE)

The leaf weights are symbolic expressions in the coupling parameters. The Piecewise branch conditions define a **hyperplane arrangement** in parameter space. EBE enumerates all feasible chambers via a BFS starting from a random interior point. Within each chamber the Piecewise conditions resolve to constants, and the DB equation reduces to a check of the form `sum(rational × exp(-β × rational)) = 0`, verified exactly by grouping exponent classes and checking rational coefficients.

If `k > ebeMaxK` (default 50), EBE falls back to the probabilistic Schwartz-Zippel check.

---

## Algorithm file format

```mathematica
$nGrid         = 3;           (* lattice side length *)
numBeta        = 1.0;         (* inverse temperature for numerical mode *)
$maxD2         = 2;           (* max squared distance for energy terms *)
$particleTypes = {1, 2, 3};   (* type multiset *)

$seedState = Module[...];     (* canonical starting state *)

$symmetryGroup = {"translation", "D4"};   (* symmetries to check/exploit *)

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
  -ebeMaxK N            Max Piecewise conditions for EBE (default 50)
  -maxDepth N           BFS bit depth limit (default 22)
  -timeLimit T          Per-state time limit in seconds (default 120)
  -verbose              Print per-rep BFS progress
```

---

## Provided examples

| File | Expected result |
|---|---|
| `single_metropolis.wl` | τ PASS, D4 PASS, DB PASS |
| `kawasaki.wl` | τ PASS, D4 PASS, DB PASS, Ergodicity FAIL (by design) |
| `vmmc_2d.wl` | τ PASS, D4 PASS, DB PASS |
| `quadratic_field.wl` | τ FAIL (absolute-position energy), DB PASS |
| `broken_variable_pool.wl` | DB FAIL — asymmetric pool size (3 or 4) |
| `broken_8way_hop.wl` | DB FAIL — asymmetric pool size (7 or 8) |
| `broken_biased_direction.wl` | DB FAIL — duplicate direction in proposal pool |
| `broken_metropolis_halfbeta.wl` | DB FAIL — accept probability uses β/2 instead of β |
| `broken_field_wrong_accept.wl` | DB FAIL — accept uses pair energy only, ignores field |

---

## Known limitations

**EBE falls back to probabilistic SZ when `k > ebeMaxK`.** Algorithms with many distinct pairwise interaction distances generate many Piecewise conditions. When `k > ebeMaxK` (default 50), EBE switches to probabilistic Schwartz-Zippel. Increase `-ebeMaxK` or accept the probabilistic guarantee.

**D4 check is algebraically exact but can be slow.** The EBE-based D4 check scales as O(feasibleRegions × pairs × generators). `vmmc_2d.wl` (k=12, 512 regions) takes ~240s for D4 EBE. Falls back to probabilistic SZ when k > ebeMaxK.

**τ-BFS cantHandle for in-body Mod.** If the algorithm normalises particle positions inside the step function using `Mod[pos, n]` on symbolic values, Mathematica may trigger internal evaluation calls that are intercepted by the BFS override, causing a cantHandle error. The checker falls back to the full state-space BFS without the translation speedup.

**NormalDistribution is discretised.** The checker models `RandomVariate[NormalDistribution[mu, sigma]]` as a discrete distribution on `Floor[nGrid/2]` integers around `mu`. For symmetric proposals (`mu=0`) this does not affect PASS/FAIL classification (the truncated Gaussian is also symmetric). For large `sigma` relative to `nGrid/2`, a warning is printed.

**Open-chamber boundary omission.** The EBE check covers all open chambers of the hyperplane arrangement. Boundaries where two or more conditions are simultaneously tight (e.g., `J1 = J2` exactly) are not tested. For standard Metropolis algorithms (transition probabilities continuous in coupling constants), any violation in a positive-measure region must appear in an adjacent open chamber, so this omission is harmless in practice.

**BFS timeout aborts the run.** If a BFS path exceeds the per-state time limit (`-timeLimit`), the checker aborts with an error. Increase `-timeLimit` if needed.
