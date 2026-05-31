# SZ-DBC: Detailed Balance Checker for Lattice MCMC

SZ-DBC automatically verifies that an MCMC algorithm satisfies **detailed balance** — the condition required for the algorithm to correctly sample the Boltzmann distribution. It is designed for lattice Monte Carlo algorithms where particles sit on a grid.

The checker works by constructing a complete **symbolic transition matrix** (intercepting all random-number calls during a BFS), verifying **translational invariance algebraically**, checking **ergodicity** via reachability, and then testing **detailed balance** using a Schwartz-Zippel probabilistic evaluation. All steps except the final SZ evaluation are exact.

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

The checker enumerates all `N!/(n₁! n₂! ...)` distinct states for the given particle type multiset on the `nGrid×nGrid` torus and verifies the count against the theoretical formula. This catches misconfigured `$nGrid` or `$particleTypes`.

### Step 2 — G-orbit computation

For a system with declared symmetry group G (translations and/or D4), many states are related by symmetry. The checker computes G-orbits of the state space and identifies one **orbit representative** per orbit, reducing subsequent BFS cost proportionally.

On a 3×3 grid with 3 distinguishable particles:
- Translation-only symmetry (`|G|=9`): 56 orbit representatives
- Translation + D4 symmetry (`|G|=72`): 8 orbit representatives

**If D4 is declared in `$symmetryGroup`, it is used for orbit reduction but is not independently verified.** If your algorithm is not actually D4-symmetric, the transition matrix will be wrong and the DB check unreliable. When in doubt, declare only `"translation"`.

### Step 3 — Translational invariance (algebraic, τ-BFS)

The checker augments every particle position in each orbit representative with a symbolic offset `{τr, τc}` and runs the BFS on these τ-augmented states. In a correctly written algorithm, all spatial operations use **pairwise differences** (e.g., `pos_i - pos_j`), so τ cancels algebraically before it enters any distance or energy:

```
(pos_i + τ) - (pos_j + τ) = pos_i - pos_j
```

If τ never appears in any leaf weight from any orbit representative, translational invariance is **algebraically certified** for the entire state space. This is not a numerical test — it is exact.

Checking from orbit representatives (not all 504 states) reduces this step from ~135 s to ~2 s.

If translation invariance fails, orbit representatives are recomputed without that symmetry before proceeding.

**Requirement:** algorithms must compute PBC via pairwise differences (`Mod[pos_i - pos_j, n]`), never via absolute positions (`Mod[pos_i, n]`). For non-translation-invariant energies (e.g., external fields), positions must be normalised before computing energy — see the notes below.

### Step 4 — BFS from orbit representatives

The BFS engine works by **intercepting all random-number calls** made by the algorithm (via Mathematica's `Block` mechanism) and replacing them with symbolic interval-tracking tokens. At each comparison `u < p` involving a random number, the BFS splits into two branches — one where the condition is true, one where it is false — accumulating the probability weight of each branch exactly. The result is a complete list of `{bit_sequence, next_state, symbolic_weight}` leaves covering every possible execution path.

Leaf weights are **symbolic expressions in the coupling constants and β**: they contain `Piecewise` expressions whose conditions are symbolic inequalities (e.g., `couplingJ[1,2,1] < couplingJ[1,2,2]`) and whose values are products of acceptance probabilities such as `1 - Exp[β(J₁ - J₂)]`.

Before any BFS, the checker installs a canonicalisation rule so that `couplingJ[b,a,d]` with `b>a` rewrites to `couplingJ[a,b,d]`, ensuring all weights use a canonical form.

If any BFS path reaches `maxDepth` (default 22) random bits without the algorithm returning, the checker reports an error — increase `-maxDepth` in this case.

### Step 5 — Ergodicity (reachability check)

After the orbit-rep BFS is complete, the checker derives the full transition graph of all states by expanding each orbit rep's leaf destinations via the G-action. A BFS from the seed state then checks whether every state is reachable.

Note that an algorithm can satisfy detailed balance but fail ergodicity — the two conditions are independent. Kawasaki exchange dynamics, for example, satisfies detailed balance but never moves particles, so only the `N! = 6` type-permutation states of the seed positions are reachable from the seed.

### Step 6 — Detailed balance check (Schwartz-Zippel)

The checker evaluates the detailed balance condition numerically at 30 randomly chosen parameter points `{β, J₁, J₂, ...}`.

For each evaluation point:
1. All symbolic leaf weights are evaluated to floating-point numbers. The Piecewise conditions (e.g., `J₁ < J₂`) become concrete comparisons once coupling values are substituted, so every weight collapses to a float.
2. The full transition matrix for all states is assembled in a single pass via a precomputed G-action integer table — no per-evaluation symbolic operations.
3. The condition `T(s→t) · exp(-β E(s)) = T(t→s) · exp(-β E(t))` is checked numerically for every communicating pair.

**Why this is sound:** for a correct algorithm, the detailed balance expression is algebraically zero as a function of the parameters. After substituting any random coupling values (which resolves all Piecewise branches), the remaining expression in β involves only products of matching exponentials that cancel exactly, giving floating-point zero to machine precision (~10⁻¹⁶). The checker tests against tolerance 10⁻⁷, giving nine orders of magnitude of safety margin. For an incorrect algorithm, the residual expression is a non-zero function of the parameters, and a random evaluation reliably detects it.

The checker reports at most **10 unique violating pairs** and stops early as soon as 10 are found, making it fast even for severely broken algorithms.

---

## Files

| File | Purpose |
|------|---------|
| `check.wls` | Command-line entry point |
| `dbc_core.wl` | Core library: BFS engine, orbit computation, τ-check, ergodicity, SZ check |
| `examples/single_metropolis.wl` | Metropolis on 2D lattice — **PASS** |
| `examples/kawasaki.wl` | Kawasaki exchange dynamics — **DB PASS, ERGODICITY FAIL** (by design) |
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
| `-mode Symbolic` | default | SZ probabilistic check |
| `-mode FullSimplify` | — | Exact symbolic check via FullSimplify (slow) |
| `-mode Numerical` | — | Numerical MCMC sampling (requires concrete coupling values) |
| `-szRepeats N` | 30 | Number of SZ evaluation points |
| `-maxDepth N` | 22 | Maximum BFS bit depth per path |
| `-verbose` | off | Print BFS progress per rep |

### Expected outputs for all examples

**`single_metropolis.wl`** — correct algorithm:
```
Translational : PASS (algebraic)
Ergodicity    : PASS (504/504 reachable)
Detailed bal. : PASS
```

**`kawasaki.wl`** — correct algorithm, non-ergodic by design:
```
Translational : PASS (algebraic)
Ergodicity    : FAIL (6/504 reachable)   ← positions never change; expected
Detailed bal. : PASS
```

**`quadratic_field.wl`** — correct algorithm with non-translation-invariant energy:
```
Translational : FAIL (recomputed without)   ← field breaks τ-invariance; expected
Ergodicity    : PASS (12/12 reachable)
Detailed bal. : PASS
```

**`broken_biased_direction.wl`**, **`broken_metropolis_halfbeta.wl`**, **`broken_field_wrong_accept.wl`** — intentionally broken:
```
Translational : PASS (algebraic)
Detailed bal. : FAIL   ← detected correctly
```

---

## Writing your own algorithm

### Required definitions

```mathematica
$nGrid         = 3;             (* lattice side length *)
$particleTypes = {1, 2, 3};    (* type multiset *)
numBeta        = 1.0;           (* inverse temperature — Numerical mode only *)
$maxD2         = 2;             (* max squared interaction distance *)
$symmetryGroup = {"translation"};  (* declared symmetries *)
$seedState     = ...;           (* canonical starting state *)

Algorithm[state_List] := ...   (* one MCMC step *)
energy[state_]        := ...   (* bare energy — no β factor *)
DynamicSymParams[states_List] := ...  (* symbolic coupling atom list *)
```

### State format

`state` is a sorted list of `{{{row,col}, type}, ...}` pairs. Positions are integers in `{1,...,nGrid}`. The checker applies PBC normalisation to every state returned by `Algorithm` — you do not need to wrap positions in the return value.

### Spatial operations: use pairwise differences

All distance and energy computations must use **pairwise differences** so that the symbolic τ offset cancels. Never apply `Mod` to an absolute coordinate stored in the state. The correct pattern:

```mathematica
(* Min-image squared distance — τ-safe *)
$myPairD2[p1_, p2_, n_] :=
  With[{dr = p1[[1]] - p2[[1]], dc = p1[[2]] - p2[[2]]},
    With[{dra = Mod[Abs[dr], n], dca = Mod[Abs[dc], n]},
      Min[dra, n-dra]^2 + Min[dca, n-dca]^2]]
```

### Non-translation-invariant energies: normalise positions internally

If the energy function is not translation-invariant (e.g., it includes a term like `fieldH · Σ rowᵢ²`), the algorithm **must** normalise particle positions before calling `energy` internally. The checker normalises the returned state for orbit purposes, but internal calls to `energy` from within `Algorithm` see whatever positions the algorithm computes. For a boundary-crossing move, an un-normalised row like 3 on `nGrid=2` gives `3² = 9`, not the correct `1² = 1`:

```mathematica
(* After computing newPos = particle[[1]] + dir, normalise before energy: *)
newPosNorm = {Mod[newPos[[1]] - 1, n] + 1, Mod[newPos[[2]] - 1, n] + 1};
With[{newState = SortBy[Append[rest, {newPosNorm, particle[[2]]}], First]},
  With[{dE = energy[newState] - energy[state]}, ...]]
```

This is required for any energy with absolute-position dependence. Pairwise-only energies are safe without normalisation because `$myPairD2` uses differences.

### Coupling constants

The energy function must use coupling atoms in canonical form `couplingJ[a, b, d2]` with `a ≤ b`. The checker installs the rule `couplingJ[b,a,d] → couplingJ[a,b,d]` for `b>a` before any BFS.

`DynamicSymParams` must return an Association listing all coupling atoms:

```mathematica
DynamicSymParams[states_List] :=
  <|"couplings"      -> {couplingJ[1,1,1], couplingJ[1,2,1], ...},
    "extraSymParams" -> {fieldH},   (* bare symbols, e.g. external fields *)
    "numericParams"  -> {}|>
```

Only list canonical (`a ≤ b`) atoms. The checker substitutes random values for everything in `"couplings"` and `"extraSymParams"` during the SZ evaluation.

### Random number calls supported

The BFS intercepts the following Mathematica functions:

| Call | Behaviour |
|------|-----------|
| `RandomReal[]` | uniform on [0,1] |
| `RandomReal[{lo, hi}]` | uniform on [lo, hi] |
| `RandomInteger[{lo, hi}]` | uniform discrete |
| `RandomChoice[list]` | uniform from list |
| `RandomChoice[weights -> elements]` | weighted choice |
| `RandomVariate[UniformDistribution[{lo,hi}]]` | uniform continuous |
| `RandomPermutation[n]`, `RandomSample[list]` | Fisher-Yates shuffle |

Any call not on this list causes the BFS to return `$dbc$cantHandle[msg]` and the checker reports an error.

---

## Known limitations

**1. SZ check is probabilistic, not exact.**
With 30 evaluation points, the probability of a false negative (declaring a broken algorithm correct) is negligible in practice. Use `-mode FullSimplify` for an exact certificate (much slower; can take hours for complex algorithms).

**2. D4 equivariance is assumed, not verified.**
Declaring `"D4"` in `$symmetryGroup` reduces orbit computation cost significantly (8 representatives instead of 56 for `nGrid=3`) but equivariance is not independently checked. An algorithm that declares D4 but breaks it will produce a wrong transition matrix and potentially a spurious PASS. Declare only `"translation"` if unsure.

**3. Only the component reachable from `$seedState` is checked.**
If the algorithm's state space has multiple disconnected ergodic components, each component must be checked separately with its own `$seedState`.

**4. Non-translation-invariant energies require position normalisation inside the algorithm.**
See the notes above. Failure to normalise is silently wrong — it causes incorrect ΔE in boundary-crossing moves, which the checker will correctly identify as a DB violation.

**5. `maxDepth` limits path length.**
If an algorithm draws more than `maxDepth` (default 22) random bits before returning, the checker reports a hard error. Increase `-maxDepth` for algorithms with long random-number sequences.

---

## Architecture

```
check.wls
  │
  ├─ Step 1: $dbcEnumerateNParticleStates + CheckErgodicity (count only)
  ├─ Step 2: $dbcAllGElems + $dbcComputeOrbits → orbit reps
  ├─ Step 3: $dbcCheckTranslational (τ-BFS from orbit reps)
  │            └─ $dbcBuildStateLeaves × |orbit reps|, τ-augmented positions
  ├─ Step 4: $dbcBuildStateLeaves × |orbit reps| → repLeaves (symbolic)
  ├─ Step 5: $dbcCheckErgodicityFromLeaves
  │            └─ BFS on derived transition graph
  └─ Step 6: $dbcSZCheckLeaves (fast numerical SZ)
               ├─ Precompute G-action integer table (once)
               ├─ For each of 30 SZ points:
               │    ├─ Evaluate leaf weights numerically (resolve Piecewise)
               │    ├─ Assemble SparseArray T matrix via G-table
               │    └─ Check T(s→t)·exp(-β E(s)) = T(t→s)·exp(-β E(t)) per pair
               └─ Report up to 10 unique violating pairs, stop early

dbc_core.wl
  ├─ SECTION 0:  State format utilities
  ├─ SECTION 0b: Coupling symmetry (couplingJ canonicalisation)
  ├─ SECTION 1:  Group actions (D4 elements, translations, G-elem application)
  ├─ SECTION 2:  Orbit computation (enumerate states, compute G-orbits)
  ├─ SECTION 3:  BFS engine (RunWithBitsAT, $dbcBuildStateLeaves)
  ├─ SECTION 4:  Translational invariance check ($dbcCheckTranslational)
  ├─ SECTION 6:  Ergodicity checks (CheckErgodicity, $dbcCheckErgodicityFromLeaves)
  ├─ SECTION 7:  Schwartz-Zippel DB check ($dbcSZCheckOne, CheckDetailedBalanceSZ)
  ├─ SECTION 7b: Direct SZ check from leaves ($dbcSZCheckLeaves)
  └─ SECTION 8:  Numerical MCMC check (RunNumericalMCMCAT, BoltzmannWeightsAT)
```

### Timings (nGrid=3, 3 distinguishable particles, translation-only, Apple M-class CPU)

| Step | Time | What it does |
|------|------|-------------|
| Kernel startup | ~3 s | wolframscript JVM + kernel init |
| State enumeration + orbit computation | ~0.05 s | 504 states, 56 orbit reps |
| τ-BFS from orbit reps | ~3 s | algebraic τ-freedom check |
| BFS from orbit reps | ~3 s | build symbolic leaf tree |
| Ergodicity check | ~1 s | reachability BFS on derived graph |
| SZ DB check | ~10 s | 30 evaluations × 10,000 pairs |
| **Total** | **~20 s** | |
