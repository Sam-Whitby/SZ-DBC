# SZ-DBC: Detailed Balance Checker for Lattice MCMC

SZ-DBC verifies that a lattice MCMC algorithm satisfies **detailed balance** — the condition required to correctly sample the Boltzmann distribution. It intercepts all random-number calls made by the algorithm, reconstructs the exact symbolic transition probabilities via BFS, and checks the DB condition algebraically over all coupling-parameter regions.

---

## Quick start

```
wolframscript -file check.wls examples/single_metropolis.wl
```

The algorithm file must define five things: `$nGrid`, `$particleTypes`, `$seedState`, `Algorithm[state_]`, `energy[state_]`, and `DynamicSymParams[states_]`. See the examples directory for complete templates.

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

### BFS engine (RunWithBitsAT)

The core of the checker intercepts every call to `RandomReal`, `RandomInteger`, `RandomChoice`, and related functions inside `Algorithm`. Each call is replaced by a **symbolic bit-reader**: one bit is consumed per binary decision, and the current path weight is updated multiplicatively. This produces a complete list of `(nextState, pathWeight)` pairs — one per leaf of the binary decision tree — with weights that sum to 1 (after the fix described below).

For `RandomChoice[list]` with `n = Length[list]`:
- Reads `k = ceil(log2(n))` bits to produce an index `0..n-1`.
- Leaves with `index >= n` are discarded (rejection sampling).
- The path weight is multiplied by `2^k / n` so that each valid leaf carries weight exactly `1/n`.

For weighted `RandomChoice[weights -> elements]`, a sequential Bernoulli decomposition is used (seqBernoulli), which is exact.

### Step 1 — State enumeration

All `N! / (n₁! n₂! ...)` states for the given particle type multiset are enumerated on the `nGrid × nGrid` torus. The count is verified against the theoretical combinatorial formula.

### Step 2 — Orbit computation

The symmetry group is used to reduce the work. Translation-only orbits are computed first. D4 rotational symmetry is verified in a later step before being used to reduce the orbit set further.

### Step 3/4 — Translational invariance (algebraic τ-BFS)

Particle positions are augmented with a symbolic offset `{τr, τc}` and the BFS is re-run on these τ-augmented states. If the algorithm uses only **pairwise differences** for all spatial decisions, τ cancels algebraically from every leaf weight. The checker verifies this symbolically — no numerical sampling.

If the algorithm contains non-translation-invariant operations (e.g., an absolute position-dependent field), the τ terms appear in the leaf weights, the check fails, and orbits are recomputed without the translation subgroup.

**Limitation**: algorithms that call `Mod[position, n]` on τ-augmented symbolic positions may trigger internal Mathematica evaluation calls. These are intercepted as `cantHandle` errors, causing the τ-BFS to fall back to the non-translation path. This is a known limitation for algorithms that normalise positions inside the step function.

### Step 5 — BFS from orbit representatives

The BFS engine runs on each orbit representative. Leaves accumulate as `{bits, nextState, weight}` triples.

### Step 5b — Rotational invariance (D4, numerical)

If `"D4"` is in `$symmetryGroup` and translation invariance passed, D4 is verified numerically: the transition matrix `T` is evaluated at several random coupling-parameter points and the checker verifies `T(s→t) = T(R·s→R·t)` for all communicating pairs. If confirmed, orbit representatives are reduced by the full D4 group and the BFS is re-run on the smaller set.

**Limitation**: this check is numerical. It is reliable in practice (a random coupling point is unlikely to accidentally satisfy the symmetry condition if the algorithm is truly asymmetric), but not a formal algebraic proof.

### Step 7 — Ergodicity

A graph BFS is run from the `$seedState` over the transition graph derived from the orbit-rep BFS leaves, checking whether every state is reachable.

### Step 8 — Detailed balance: Exhaustive Branch Enumeration (EBE)

The leaf weights from the BFS are polynomial/piecewise expressions in the coupling parameters `{J₁, J₂, ...}`. The Piecewise branch conditions define a **hyperplane arrangement** in parameter space. EBE enumerates all feasible chambers of this arrangement via a chamber-adjacency BFS (starting from a random interior point, navigating by single-condition sign-flips). Within each chamber the Piecewise conditions resolve to constants, and the DB equation reduces to a check of the form `sum(rational × exp(-β × rational)) = 0`, which is verified exactly by grouping exponent classes and checking rational coefficients.

If the number of branch conditions exceeds `-ebeMaxK` (default 50), EBE falls back to the probabilistic Schwartz-Zippel check.

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
| `RandomChoice[weights->elems]` | Sequential Bernoulli decomposition |
| `RandomVariate[NormalDistribution[mu,sigma]]` | Discretised to nearby integers |
| `RandomPermutation`, `RandomSample` | Fisher-Yates via `RandomInteger` |

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
| `broken_variable_pool.wl` | τ PASS, DB FAIL — demonstrates weight bug |
| `broken_8way_hop.wl` | τ PASS, DB FAIL — demonstrates weight bug (k=3 bracket) |

The two `broken_*` files contain algorithms that genuinely violate detailed balance due to asymmetric proposal pool sizes. They were specifically designed to expose the `RandomChoice` weight-correction fix.

---

## Known limitations

**D4 check is numerical.** A formal algebraic proof of rotational symmetry is not yet implemented. The numerical check has never been observed to give a false positive in practice, but is not formally guaranteed.

**τ-BFS cantHandle for in-body Mod.** If the algorithm normalises particle positions inside the step function using `Mod[pos, n]` on symbolic values, Mathematica may trigger internal evaluation calls that are intercepted by the BFS override, causing a cantHandle error. The checker falls back to the full state-space BFS without the translation speedup.
