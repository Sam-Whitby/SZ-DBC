# SZ-DBC: Detailed Balance Checker for Lattice MCMC

SZ-DBC automatically verifies that an MCMC algorithm satisfies **detailed balance** — the condition required for the algorithm to correctly sample the Boltzmann distribution. It is designed for lattice Monte Carlo algorithms where particles sit on a grid.

The checker works by constructing a complete **symbolic transition matrix** (intercepting all random-number calls during a BFS), verifying **translational invariance algebraically**, verifying **rotational invariance numerically** (and reducing orbit representatives accordingly), checking **ergodicity** via reachability, and testing **detailed balance** using Exhaustive Branch Enumeration (EBE) — an exact method that systematically covers all coupling-constant parameter regions using a hyperplane-arrangement BFS.

---

## What is detailed balance?

An MCMC algorithm runs in a Markov chain: at each step it moves from state `s` to a new state `t` with probability `T(s→t)`. The algorithm samples the Boltzmann distribution `π(s) ∝ exp(-β E(s))` correctly if and only if:

```
T(s→t) · π(s) = T(t→s) · π(t)    for all pairs of states (s, t)
```

Verifying this analytically for a complex algorithm is very difficult by hand. This checker does it automatically.

---

## How the checker works

### State representation

Every state is a **sorted list of `{position, type}` pairs**:

```mathematica
{{{1,1}, 1}, {{1,2}, 2}, {{2,3}, 3}}
```

Each particle has a 2D position `{row, col}` (integers 1 to `nGrid`) and an integer type label. The list is always sorted by position. The checker applies periodic boundary condition (PBC) normalisation to all returned states automatically.

### Step 1 — State enumeration and count check

The checker enumerates all `N!/(n₁! n₂! ...)` distinct states for the given particle type multiset on the `nGrid×nGrid` torus and verifies the count against the theoretical formula.

### Step 2 — Initial orbit computation (translation-only)

The checker always starts with **translation-only orbits**, regardless of what is declared in `$symmetryGroup`. For a 3×3 grid with 3 distinguishable particles, this gives 56 orbit representatives.

D4 symmetry (if declared) is verified algebraically after the BFS (Step 5b) before being used. This ensures D4 is never silently assumed.

### Step 3 — Translational invariance (algebraic, τ-BFS)

The checker augments every particle position with a symbolic offset `{τr, τc}` and runs the BFS on these τ-augmented states. In a correctly written algorithm, all spatial operations use **pairwise differences**, so τ cancels algebraically:

```
(pos_i + τ) - (pos_j + τ) = pos_i - pos_j
```

If τ never appears in any leaf weight from any orbit representative, translational invariance is **algebraically certified** for the entire state space. This is exact, not numerical.

Checking from orbit representatives reduces this step from ~135 s to ~3 s.

### Step 4 — BFS from orbit representatives

The BFS engine **intercepts all random-number calls** (via Mathematica's `Block` mechanism) and replaces them with symbolic interval-tracking tokens. At each comparison `u < p`, the BFS splits into two branches — accumulating the probability weight of each branch exactly. The result is a complete list of `{bit_sequence, next_state, symbolic_weight}` leaves.

Leaf weights are **symbolic expressions in coupling constants and β**: they contain `Piecewise` expressions whose conditions are symbolic inequalities (e.g., `couplingJ[1,2,1] < couplingJ[1,2,2]`) and whose values are products of acceptance probabilities.

### Step 5b — Rotational invariance check (D4, numerical)

If `"D4"` is declared in `$symmetryGroup` **and** translation invariance passed, the checker verifies D4 symmetry numerically using the T matrix already built in Step 4.

**Method:** Build T numerically at `nReps` random coupling points (default 10). For each non-identity D4 element R (rotations and reflections) and each communicating pair (s, t): check `|T(s→t) - T(R·s → R·t)| < 10⁻⁶`. The check tests all 7 non-identity D4 elements simultaneously.

**If full D4 passes:** Recompute orbits with the full D4+translation group (72 elements → ~8 orbit reps for nGrid=3). Rerun BFS from the smaller representative set. Use these 8 reps for all subsequent steps.

**If only C4 (rotations) pass but reflections fail:** Recompute orbits with C4+translation (36 elements → ~14 orbit reps). Rerun BFS from ~14 reps.

**If neither passes:** Keep translation-only orbits (56 reps), warn, and continue.

This step is the key to the orbit reduction that speeds up the DB check for symmetric algorithms. It is a clean algebraic/numerical boundary: translation invariance is verified algebraically (τ-cancellation), rotational invariance is verified by numerically comparing T matrix entries across symmetry-related states.

**Note on the numerical check:** The check uses random coupling points (SZ-style). For any physically reasonable algorithm, 10 random points give overwhelming confidence. A false positive (rotationally non-invariant algorithm passing the check) would require the T matrix symmetry to hold at 10 independent random coupling points by coincidence — negligibly unlikely.

### Step 6 — Ergodicity (reachability check)

After the (possibly D4-reduced) orbit-rep BFS is complete, the checker derives the full transition graph by expanding each orbit rep's leaf destinations via the G-action. A BFS from the seed state checks whether every state is reachable.

### Step 7 — Detailed balance check (EBE exact)

This is the heart of the checker. The key challenge is that symbolic leaf weights contain `Piecewise` expressions. Different coupling-constant parameter regions give genuinely different transition matrices, and a DB violation might only manifest in a specific region.

#### Exhaustive Branch Enumeration (EBE) — the default

**Phase 1 — Condition extraction.** Scan all symbolic leaf weights and collect the `k` distinct Piecewise branch conditions. Each is a comparison involving coupling-constant atoms (e.g., `couplingJ[1,2,1] - couplingJ[1,2,2] < 0`).

**Phase 2 — Region enumeration via arrangement BFS.**

The `k` conditions define `k` hyperplanes in the coupling-constant parameter space. These hyperplanes divide the space into `T` non-empty regions (cells of the hyperplane arrangement). We need a rational representative point `J*` for every cell.

A classical result from computational geometry states that the adjacency graph of the cells of any hyperplane arrangement is **connected** under face-adjacency (two cells are adjacent iff they share a codimension-1 face — i.e., differ in exactly one sign bit). This means a BFS from any starting cell reaches all cells.

The arrangement BFS:
1. Find any starting rational point `J*₀` via `FindInstance[True, params, Rationals]`
2. Compute its sign pattern σ₀ (which side of each hyperplane it lies on)
3. BFS: for each dequeued σ, try all `k` single-bit flips
4. For each unvisited flip σ': call `FindInstance[constraints(σ'), params, Rationals]`
   - Feasible → real cell: add to queue and record `J*(σ')`
   - Infeasible → empty region: mark visited, do not queue

**Total `FindInstance` calls ≤ T × k**, where T = number of real cells. This is vastly more efficient than exhaustive enumeration of all `2^k` sign patterns when T << 2^k.

For example, with k=18 (single_metropolis): exhaustive would require 262,144 calls. With T=132 real cells, the BFS uses only 1,932 calls — a **136× reduction**.

**Phase 3 — Exact DB check per region.** For each feasible region with coupling point `J*`:
1. Substitute `J*` into every leaf weight. All `Piecewise` conditions become concrete `True`/`False` and auto-simplify. β remains symbolic.
2. Assemble T using a precomputed G-action table.
3. For each communicating pair (s, t): form `T(s→t)·Exp[-β·E(s)] - T(t→s)·Exp[-β·E(t)]`.
4. Check exactly using `$dbcIsExpZero`: group terms by rational exponent coefficient of β; verify all rational coefficient sums are zero.

**EBE reports an exact certificate:** PASS means detailed balance holds for every possible coupling-constant configuration; FAIL includes the specific coupling values as a counterexample.

#### SZ probabilistic fallback (k > ebeMaxK)

When `k > ebeMaxK` (default 50), falls back to probabilistic SZ: substitute 30 random rational-valued coupling points and check with tolerance `10⁻⁷`. For Metropolis-type algorithms (k=18), any violation appears across essentially all coupling regions, so SZ is reliable despite not being exhaustive.

---

## Files

| File | Purpose |
|------|---------|
| `check.wls` | Command-line entry point |
| `dbc_core.wl` | Core library: BFS engine, orbit computation, τ-check, D4 verification, ergodicity, EBE/SZ check |
| `examples/single_metropolis.wl` | Metropolis on 2D lattice — **PASS** |
| `examples/kawasaki.wl` | Kawasaki exchange dynamics — **DB PASS, ERGODICITY FAIL** (by design) |
| `examples/vmmc_2d.wl` | VMMC cluster moves — **PASS** |
| `examples/quadratic_field.wl` | Metropolis with non-translation-invariant energy — **τ FAIL, DB PASS** |
| `examples/broken_biased_direction.wl` | Biased move direction bug — **DB FAIL** |
| `examples/broken_metropolis_halfbeta.wl` | Wrong β in acceptance — **DB FAIL** |
| `examples/broken_field_wrong_accept.wl` | Ignores field in acceptance — **DB FAIL** |

---

## How to run

### Prerequisites

Wolfram Mathematica or Wolfram Engine (`wolframscript` must be on your PATH).

### Running the checker

```bash
wolframscript -file check.wls examples/single_metropolis.wl
```

### Options

| Flag | Default | Description |
|------|---------|-------------|
| `-mode Symbolic` | default | EBE exact check (falls back to SZ when k > ebeMaxK) |
| `-mode SZPure` | — | Probabilistic SZ check (faster; use when EBE is slow) |
| `-mode FullSimplify` | — | Exact via FullSimplify (very slow; for reference) |
| `-mode Numerical` | — | Numerical MCMC sampling (requires concrete coupling values) |
| `-ebeMaxK N` | 50 | Max branch conditions k for EBE; SZ fallback above this |
| `-szRepeats N` | 30 | Number of SZ evaluation points |
| `-maxDepth N` | 22 | Maximum BFS bit depth per path |
| `-timeLimit T` | 120 | Per-state time limit in seconds |
| `-verbose` | off | Print BFS progress per rep |

### Timings (nGrid=3, 3 distinguishable particles, Apple M-class CPU)

All examples below declare `$symmetryGroup = {"translation", "D4"}` where D4 holds.

| Example | k | Reps (after D4) | EBE regions | LP calls | Total time |
|---|---|---|---|---|---|
| `broken_field_wrong_accept.wl` (nGrid=2) | 2 | 3 (no D4) | 3 | 4 | ~5s |
| `quadratic_field.wl` (nGrid=2, τ fails) | 6 | 12 (no D4) | 12 | 48 | ~4s |
| `kawasaki.wl` | 6 | 8 (D4 ✓) | 12 | 48 | ~9s |
| `broken_biased_direction.wl` | 18 | 56 (no D4) | 132 | 1,932 | ~28s |
| `broken_metropolis_halfbeta.wl` | 18 | 56 (no D4) | 132 | 1,932 | ~26s |
| `single_metropolis.wl` | 18 | 8 (D4 ✓) | 132 | 1,932 | ~49s |
| `vmmc_2d.wl` | 12 | 8 (D4 ✓) | 512 | 2,048 | ~3m 9s |

Notes:
- "Total time" includes kernel startup (~3s), τ-BFS (~3-13s), BFS, D4 check, ergodicity, and DB check
- "LP calls" = `FindInstance` calls in EBE Phase 2 (arrangement BFS). Compare to exhaustive: k=18 → 262,144; k=12 → 4,096
- For `vmmc_2d.wl`, EBE DB check alone is ~141s (was ~215s with translation-only orbits). Use `-mode SZPure` for a fast ~13s probabilistic check
- `broken_biased_direction` and `broken_metropolis_halfbeta` declare translation-only (D4 check not triggered)

### Expected outputs for all examples

**`single_metropolis.wl`** — correct algorithm, D4 verified:
```
Translational : PASS (algebraic)
D4 rotational : D4 (verified)
Ergodicity    : PASS (504/504 reachable)
Detailed bal. : PASS
```

**`kawasaki.wl`** — correct algorithm, non-ergodic by design, D4 verified:
```
Translational : PASS (algebraic)
D4 rotational : D4 (verified)
Ergodicity    : FAIL (6/504 reachable)   ← positions never change; expected
Detailed bal. : PASS
```

**`quadratic_field.wl`** — correct algorithm with non-translation-invariant energy:
```
Translational : FAIL (recomputed without)   ← field breaks τ-invariance; expected
D4 rotational : not checked (τ failed)
Ergodicity    : PASS (12/12 reachable)
Detailed bal. : PASS
```

**`broken_biased_direction.wl`**, **`broken_metropolis_halfbeta.wl`** — broken:
```
Translational : PASS (algebraic)
D4 rotational : not declared
Detailed bal. : FAIL   ← detected exactly (EBE, k=18, 132 regions)
```

**`broken_field_wrong_accept.wl`** — broken:
```
Translational : PASS (algebraic)
D4 rotational : not declared
Detailed bal. : FAIL   ← detected exactly in 3 feasible regions
```

---

## Writing your own algorithm

### Required definitions

```mathematica
$nGrid         = 3;             (* lattice side length *)
$particleTypes = {1, 2, 3};    (* type multiset *)
numBeta        = 1.0;           (* inverse temperature — Numerical mode only *)
$maxD2         = 2;             (* max squared interaction distance *)
$symmetryGroup = {"translation", "D4"};  (* declare D4 if your algorithm is rotationally symmetric *)
$seedState     = ...;           (* canonical starting state *)

Algorithm[state_List] := ...   (* one MCMC step *)
energy[state_]        := ...   (* bare energy — no β factor *)
DynamicSymParams[states_List] := ...  (* symbolic coupling atom list *)
```

### Declaring D4 symmetry

Declare `"D4"` in `$symmetryGroup` if your algorithm satisfies T(R·s → R·t) = T(s→t) for all 90° rotations and reflections R. The checker will:
1. Run the initial BFS from 56 translation-only orbit representatives
2. Verify D4 numerically (builds T matrix, checks symmetry at 10 random coupling points)
3. If verified: reduce to 8 orbit representatives and rerun BFS — speeding up the DB check by ~4–8×

If D4 fails the numerical check, the checker warns and continues with translation-only orbits. No correctness is lost; only the performance benefit is foregone.

Typical cases where D4 holds: isotropic pairwise interactions (energy depends only on |distance|, not direction), symmetric proposal distribution (uniform over all directions or nearest neighbours).

Typical cases where D4 fails: directional external fields, biased move proposals.

### Spatial operations: use pairwise differences

All distance and energy computations must use **pairwise differences** so that the symbolic τ offset cancels:

```mathematica
$myPairD2[p1_, p2_, n_] :=
  With[{dr = p1[[1]] - p2[[1]], dc = p1[[2]] - p2[[2]]},
    With[{dra = Mod[Abs[dr], n], dca = Mod[Abs[dc], n]},
      Min[dra, n-dra]^2 + Min[dca, n-dca]^2]]
```

### Non-translation-invariant energies: normalise positions internally

If the energy includes an absolute-position term (e.g., `fieldH · Σ rowᵢ²`), normalise positions before calling `energy`:

```mathematica
newPosNorm = {Mod[newPos[[1]] - 1, n] + 1, Mod[newPos[[2]] - 1, n] + 1};
With[{newState = SortBy[Append[rest, {newPosNorm, particle[[2]]}], First]},
  With[{dE = energy[newState] - energy[state]}, ...]]
```

### Coupling constants

Use canonical form `couplingJ[a, b, d2]` with `a ≤ b`. The checker installs `couplingJ[b,a,d] → couplingJ[a,b,d]` for `b>a`.

`DynamicSymParams` must return:

```mathematica
DynamicSymParams[states_List] :=
  <|"couplings"      -> {couplingJ[1,1,1], couplingJ[1,2,1], ...},
    "extraSymParams" -> {fieldH},   (* bare symbols, e.g. external fields *)
    "numericParams"  -> {}|>
```

### Random number calls supported

| Call | Behaviour |
|------|-----------|
| `RandomReal[]` | uniform on [0,1] |
| `RandomReal[{lo, hi}]` | uniform on [lo, hi] |
| `RandomInteger[{lo, hi}]` | uniform discrete |
| `RandomChoice[list]` | uniform from list |
| `RandomChoice[weights -> elements]` | weighted choice |
| `RandomVariate[UniformDistribution[{lo,hi}]]` | uniform continuous |
| `RandomPermutation[n]`, `RandomSample[list]` | Fisher-Yates shuffle |

---

## Known limitations

**1. EBE is exact for k ≤ ebeMaxK; SZ is probabilistic above.**
With the arrangement BFS, k=18 now takes ~30s (exact). For k > 50 (default fallback), 30 SZ evaluation points give high practical confidence for Metropolis-type algorithms. Use `-ebeMaxK 100` to push EBE higher if needed.

**2. D4 verification is numerical, not algebraic.**
The τ-invariance check is algebraically exact (τ cancels symbolically). The D4 check is numerical: 10 random coupling points. A false positive — a non-D4 algorithm passing the check — is negligible in practice (probability ~10⁻⁶⁰ or smaller) but is not ruled out by proof. If you need an exact certificate of D4 equivariance, you would need a symbolic check (not currently implemented).

**3. Only the component reachable from `$seedState` is checked.**
If the algorithm's state space has multiple disconnected ergodic components, check each separately with its own `$seedState`.

**4. Non-translation-invariant energies require position normalisation inside the algorithm.**
Failure to normalise is silently wrong and will be detected as a DB violation.

**5. `maxDepth` limits path length.**
If an algorithm draws more than `maxDepth` (default 22) random bits before returning, the checker reports a hard error. Increase `-maxDepth` for algorithms with long random-number sequences.

---

## Architecture

```
check.wls
  │
  ├─ Step 1:  $dbcEnumerateNParticleStates + CheckErgodicity (count only)
  ├─ Step 2:  $dbcAllGElems + $dbcComputeOrbits → translation-only orbit reps (56)
  ├─ Step 3:  $dbcCheckTranslational (τ-BFS from orbit reps) — algebraically exact
  ├─ Step 4:  $dbcBuildStateLeaves × |orbit reps| → repLeaves (symbolic)
  ├─ Step 5b: $dbcCheckRotational (D4 numerical check, if declared)
  │            └─ if PASS: recompute orbits (→ 8 reps), rerun BFS
  ├─ Step 6:  $dbcCheckErgodicityFromLeaves (BFS on derived transition graph)
  └─ Step 7:  $dbcEBECheckLeaves (default) or $dbcSZCheckLeaves (fallback k > 50)
               │
               ├─ Phase 1: Extract k distinct Piecewise conditions
               ├─ Phase 2: Arrangement BFS (O(T×k) LP calls instead of O(2^k))
               │            └─ Start from any rational point, BFS via single-bit flips
               │               Connected-arrangement theorem guarantees completeness
               └─ Phase 3: For each of T feasible regions:
                    ├─ Substitute rational J* into leaf weights (Piecewise auto-resolves)
                    ├─ Assemble T via precomputed G-action table
                    └─ Check T(s→t)·exp(-β E(s)) = T(t→s)·exp(-β E(t)) per pair
                         └─ $dbcIsExpZero: exact rational arithmetic, no floating point

dbc_core.wl
  ├─ SECTION 0:   State format utilities
  ├─ SECTION 0b:  Coupling symmetry (couplingJ canonicalisation)
  ├─ SECTION 1:   Group actions (D4 elements, translations, G-elem application)
  ├─ SECTION 2:   Orbit computation ($dbcComputeOrbits)
  ├─ SECTION 3:   BFS engine (RunWithBitsAT, $dbcBuildStateLeaves)
  ├─ SECTION 4:   Translational invariance check ($dbcCheckTranslational)
  ├─ SECTION 6:   Ergodicity checks
  ├─ SECTION 7:   Schwartz-Zippel DB check
  ├─ SECTION 7b:  Direct SZ check from leaves ($dbcSZCheckLeaves)
  ├─ SECTION 7c:  EBE exact check with arrangement BFS ($dbcEBECheckLeaves)
  ├─ SECTION 7d:  Rotational invariance check ($dbcCheckRotational)
  └─ SECTION 8:   Numerical MCMC check
```

### Key algorithmic improvements in this version

**Arrangement BFS (EBE Phase 2).** The original EBE Phase 2 tested all `2^k` sign patterns via `FindInstance`, discarding the infeasible ones. The arrangement BFS instead exploits the connected-arrangement theorem: the cells of any hyperplane arrangement form a connected graph under face-adjacency. Starting from one cell, BFS via single-bit flips reaches all cells. Total LP calls ≤ T × k (T = real cells, k = conditions), compared to 2^k. For k=18, T=132: 1,932 calls vs 262,144 — a 136× speedup.

**Rotational invariance verification.** Previously, D4 symmetry was either assumed (declared but unchecked) or unused. Now, when `"D4"` is declared, the checker builds the T matrix from the translation-only BFS and numerically verifies T(s→t) = T(R·s → R·t) for all communicating pairs and all 7 non-identity D4 elements. If verified, orbits are recomputed with the full D4+translation group (~8 representatives), BFS is rerun from the smaller set, and the DB check operates on ~8 orbit reps instead of 56 — reducing Phase 3 cost by approximately 4–8×.

**Orbit map correctness fix.** The `$dbcComputeOrbits` function previously stored orbit map entries as `g` where `g·s = img` (s = arbitrary first unvisited state). The orbit expansion requires `h` where `h·rep = img` (rep = canonical representative). These coincide for translation groups (where the canonical rep is always the first-enumerated state) but diverge for D4 (where the canonical rep is the lexicographically smallest D4-image, not necessarily the first-visited state). The fix recomputes the orbit map from the canonical rep after finding it, at a cost of one additional group-action pass per orbit.
