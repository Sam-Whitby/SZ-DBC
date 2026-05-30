# SZ-DBC: Detailed Balance Checker for Lattice MCMC

SZ-DBC automatically verifies that an MCMC algorithm satisfies **detailed balance** — the condition required for the algorithm to correctly sample the Boltzmann distribution. It is designed for lattice Monte Carlo algorithms where particles sit on a grid.

The checker produces algebraic or probabilistic proofs of correctness, not just numerical tests. It verifies translational and rotational symmetry exactly, checks detailed balance across the entire reachable state space in a single run, and verifies ergodicity by constructing the full transition graph.

---

## What is detailed balance?

An MCMC algorithm runs in a Markov chain: at each step, it moves from a state `s` to a new state `t` with some probability `T(s→t)`. The algorithm samples the Boltzmann distribution `π(s) ∝ exp(-β E(s))` correctly if and only if it satisfies **detailed balance**:

```
T(s→t) · π(s) = T(t→s) · π(t)    for all pairs of states (s, t)
```

This says that the flow of probability from `s` to `t` equals the reverse flow from `t` to `s`, weighted by their equilibrium probabilities. Verifying this analytically for a complex algorithm (especially one with cluster moves, like VMMC) is very difficult by hand. This checker does it automatically.

---

## How the checker works

### State representation

Every state is a **sorted list of `{position, type}` pairs**:

```mathematica
{{{1,1}, 1}, {{1,2}, 2}, {{2,3}, 3}}
```

Each particle has a 2D position `{row, col}` (integers from 1 to `nGrid`) and an integer type label. The list is sorted by position. This representation allows the checker to augment positions with symbolic offsets (used for translational invariance verification) while keeping the arithmetic transparent.

### Step 1 — State enumeration and combinatorial count

The checker enumerates all `N!/(n₁! n₂! ...)` distinct states for the given particle type multiset on the `nGrid × nGrid` torus (where `n₁, n₂, ...` are the multiplicities of each type) using `Permutations`. This count is verified against the theoretical formula as a configuration sanity check — if there is a mismatch, `$nGrid` or `$particleTypes` is misconfigured.

### Step 2 — G-orbit computation

For a system with symmetry group G (translations × D4), many states are related by symmetry. The checker computes the G-orbits of the state space and identifies one **orbit representative** per orbit. On a `3×3` grid with 3 distinguishable particles, there are 504 states but only 8 orbit representatives (a 63× reduction). This step is O(|G| × |orbit count|) and takes < 0.1 s.

### Step 3 — Translational invariance (algebraic, τ-BFS from orbit reps)

The checker augments every particle position in each orbit representative with a symbolic offset `{τr, τc}`:

```
{{1+τr, 1+τc}, 1}, {{1+τr, 2+τc}, 2}, ...
```

It then runs `$dbcBuildStateLeaves` for each orbit rep with these τ-augmented positions and scans every transition probability for τ using `FreeQ`. Because all spatial computations in a correctly written algorithm use **pairwise differences** (e.g., `pos_i - pos_j`), the symbol τ cancels algebraically before it ever enters a distance, energy, or acceptance probability:

```
(pos_i + τ) - (pos_j + τ) = pos_i - pos_j
```

If τ never appears in any transition probability from any orbit rep, translational invariance is **algebraically certified** for the entire reachable state space. This is not a numerical test — it is exact. Checking from orbit representatives (rather than all 504 states) is equivalent under G-symmetry and reduces this step from ~135 s to ~2 s.

**Requirement on the algorithm**: the algorithm must compute periodic boundary conditions by wrapping *differences*, not absolute positions. `Mod[pos_i - pos_j, n]` is correct; `Mod[pos_i, n] - Mod[pos_j, n]` is not.

If translation invariance fails, the orbit representatives are recomputed using only the D4 subgroup (removing the 9× translation speedup) before proceeding.

### Step 4 — BFS from orbit representatives

The BFS engine works by **intercepting all random number calls** made by the algorithm (via Mathematica's `Block` mechanism) and replacing them with symbolic interval-tracking tokens. At each random number comparison `u < p`, the BFS splits into two branches — one where the condition is true, one where it is false — accumulating the probability of each branch as an exact symbolic weight. The result is a complete list of `{bit_sequence, next_state, symbolic_weight}` leaves covering all possible execution paths.

Transition probabilities are **symbolic functions of the coupling constants and β** — they are never evaluated numerically during the BFS. The weights contain `Piecewise` expressions whose conditions are symbolic inequalities (e.g., `couplingJ[1,2,1] < couplingJ[1,2,2]`) and whose values are products of acceptance probabilities such as `1 - Exp[β(J₁ - J₂)]`.

**Coupling symmetry**: before running any BFS, the checker installs a canonicalisation rule so that `couplingJ[b,a,d]` with `b>a` automatically rewrites to `couplingJ[a,b,d]`. This ensures all leaf weights and energy expressions use a canonical (a≤b) form, so the Schwartz-Zippel substitution covers every atom that actually appears.

**Completeness**: every possible execution path up to `maxDepth` random bits is explored. For `vmmc_2d.wl` with 3 particles, the maximum path depth is ~17 bits, well within the default `maxDepth = 22`. A warning is printed if any paths hit the depth limit.

### Step 5 — D4 rotational invariance (generator testing, numerical)

For each orbit representative `s₀`, the checker runs two additional BFS runs:
- One from `rot90(s₀)` — the state rotated 90° clockwise
- One from `reflect(s₀)` — the state reflected across the main diagonal

It then checks, for each transition `T(s₀ → t)`, whether:

```
T(rot90(s₀) → rot90(t)) ≈ T(s₀ → t)
T(reflect(s₀) → reflect(t)) ≈ T(s₀ → t)
```

The comparison uses **SZ-style numerical evaluation**: both symbolic expressions are evaluated at 3 randomly chosen parameter points `{β, J₁, J₂, ...}`. If they agree numerically at all points (within floating-point tolerance), D4 invariance is verified for that transition. This is more robust than structural equality (`===`) which can fail for mathematically identical but structurally different Piecewise expressions.

By group theory, checking these two generators (rotation and reflection) is sufficient to certify all 8 elements of D4. This requires only 2 extra BFS runs per representative.

### Step 6 — Ergodicity (reachability from orbit-BFS transition graph)

After the orbit-rep BFS is complete, the checker derives the **full transition graph** of all 504 states by expanding each orbit rep's leaf destinations via the G-action: for each orbit rep `s₀` and each leaf destination `t`, every state `s = g·s₀` in the orbit can reach `g·t` in one step.

A BFS is then run on this derived graph starting from the first enumerated state (the seed). If all 504 states are reachable, the algorithm is **ergodic** — any state can eventually reach any other state. If some states are not reachable, a failure is reported with the reachable count.

This check uses the already-computed orbit BFS data with no additional algorithm calls, completing in ~2 s.

### Step 7 — Detailed balance check (SZPure, default)

The checker uses a **Schwartz-Zippel probabilistic zero-test**. Rather than building a large symbolic transition matrix and applying `FullSimplify` (which can take hours), it evaluates the detailed balance condition numerically at 30 randomly chosen parameter points `{β, J₁, J₂, ...}`.

For each evaluation point:
1. The symbolic leaf weights (products of acceptance probabilities) evaluate immediately to floating-point numbers, because all Piecewise conditions (e.g., `J₁ < J₂`) become concrete comparisons once coupling values are substituted.
2. The full transition matrix for all 504 states is assembled via a precomputed G-action integer table — no symbolic operations needed.
3. The detailed balance condition `T(s→t) · exp(-β E(s)) = T(t→s) · exp(-β E(t))` is checked numerically for every pair.

If the algorithm satisfies detailed balance, this expression is identically zero as a function of the parameters, and all 30 evaluations return zero. If there is a violation, it appears as a nonzero polynomial in the parameters, which a random evaluation detects with high probability (by the Schwartz-Zippel lemma). With 30 evaluations and coupling constants drawn from a bounded rational range, the probability of a false negative is negligibly small.

**Why this is faster than FullSimplify**: numerical substitution makes all Piecewise conditions collapse immediately to True/False, reducing the weight to a simple float product. `FullSimplify` on the raw symbolic sum can take hours.

---

## Files

| File | Purpose |
|------|---------|
| `check.wls` | Command-line entry point. Load and run this with `wolframscript`. |
| `dbc_core.wl` | Core library: BFS engine, symmetry checks, SZ checker, numerical MCMC. |
| `examples/vmmc_2d.wl` | Example algorithm: Virtual Move Monte Carlo on a 2D square lattice. |

---

## How to use

### Prerequisites

- Wolfram Mathematica / Wolfram Engine (for `wolframscript`)

### Running the checker

```bash
wolframscript -file check.wls examples/vmmc_2d.wl
```

With verbose BFS progress output:

```bash
wolframscript -file check.wls examples/vmmc_2d.wl -verbose
```

Using `FullSimplify` instead of SZPure (much slower, exact result):

```bash
wolframscript -file check.wls examples/vmmc_2d.wl -mode FullSimplify
```

Numerical MCMC check only (fast sanity check, not a proof — requires concrete coupling values):

```bash
wolframscript -file check.wls examples/vmmc_2d.wl -mode Numerical
```

Adjust SZ evaluation count (default 30):

```bash
wolframscript -file check.wls examples/vmmc_2d.wl -szRepeats 50
```

### Expected output (vmmc_2d.wl, nGrid=3, 3 particles)

```
================================================================
  SZ-DBC  —  vmmc_2d.wl
================================================================
  nGrid      : 3
  particles  : {1, 2, 3}
  symmetry   : {translation, D4}
  mode       : Symbolic

Step 1: State enumeration
  States found : 504  (0.01s)
  State count  : OK  (found 504, theoretical 504)

Step 2: G-orbit computation
  |G| = 72
  Orbit reps : 8  (0.04s)

Step 4: Translational invariance (τ-BFS from orbit reps)
  States checked : 8 orbit reps  (2.0s)
  PASS  — translation invariance algebraically verified

Step 5: BFS from orbit reps
  BFS done: 8 reps  (1.9s)

Step 6: D4 generator testing
  PASS  — D4 invariance verified (8 reps × 2 generators)
  Time: 4.0s

Step 7: Ergodicity (reachability check)
  Reachable : 504/504  (1.8s)
  PASS  — all states reachable from seed

Step 8: Detailed balance (SZPure direct)
  Time: 5.3s
  PASS  — detailed balance satisfied for all 504 states

================================================================
  SUMMARY
================================================================
  Translational : PASS (algebraic)
  D4 rotational : PASS
  State count   : OK (504/504)
  Ergodicity    : PASS (504/504 reachable)
  Detailed bal. : PASS
================================================================
```

Total runtime: approximately **15 seconds** on a modern laptop. Exit code 0 = all checks passed; exit code 1 = failure.

---

## Writing your own algorithm

Your algorithm file must define the following globals:

```mathematica
$nGrid         = 3;             (* lattice side length *)
$particleTypes = {1, 2, 3};    (* type multiset for this run *)
numBeta        = 1.0;           (* inverse temperature *)
$maxD2         = 2;             (* max squared interaction distance *)
$symmetryGroup = {"translation", "D4"};  (* declared symmetries *)
$seedState     = ...;           (* canonical starting state *)

Algorithm[state_List] := ...   (* MCMC step *)
energy[state_]        := ...   (* bare energy (no β) *)

DynamicSymParams[states_List] := ...  (* coupling atom list *)
```

### State format rules

- `state` is a sorted list of `{{{row,col}, type}, ...}` pairs
- Positions are integers in `{1, ..., nGrid}` (1-indexed)
- The checker applies PBC normalisation to returned states — you do **not** need to `Mod` positions in the algorithm
- All spatial operations must use **pairwise differences** (so that τ cancels). Never compute `Mod[absolute_position, nGrid]` directly on a stored coordinate

### Coupling constants

The energy and acceptance probability functions must use coupling atoms in **canonical form**: `couplingJ[a, b, d2]` with `a ≤ b`. The checker automatically installs the rule `couplingJ[b, a, d] → couplingJ[a, b, d]` for `b > a` before any BFS, so algorithms that don't canonicalise the argument order will still work correctly. However, `DynamicSymParams` must list only the canonical (a ≤ b) atoms — the checker substitutes values only for those.

### Random number calls supported

The BFS intercepts all standard Mathematica random functions:

- `RandomReal[]` — uniform on [0,1]
- `RandomReal[{lo, hi}]` — uniform on [lo, hi]  
- `RandomInteger[{lo, hi}]` — uniform discrete
- `RandomChoice[list]` — uniform from list
- `RandomChoice[weights -> elements]` — weighted choice
- `RandomVariate[UniformDistribution[{lo,hi}]]` — uniform continuous
- `RandomPermutation[n]`, `RandomSample[list]`

### Algorithm requirements

- `Algorithm[state]` must return a valid state (same type multiset, positions in any integer range — checker normalises)
- The energy function `energy[state]` must be a symbolic expression in the coupling atoms listed by `DynamicSymParams` — β is **not** included (it is multiplied in by the checker)
- `DynamicSymParams[states]` should return `<|"couplings" -> {couplingJ[1,1,1], couplingJ[1,2,1], ...}, "numericParams" -> {}|>` listing every canonical (a ≤ b) coupling atom that appears in the energy

### Geometry helpers (from vmmc_2d.wl)

```mathematica
(* Min-image squared distance — τ-safe *)
$pairD2[p1_, p2_, n_] :=
  With[{dr0 = p1[[1]] - p2[[1]], dc0 = p1[[2]] - p2[[2]]},
    With[{dra = Mod[Abs[dr0], n], dca = Mod[Abs[dc0], n]},
      Min[dra, n-dra]^2 + Min[dca, n-dca]^2]]

(* Find neighbour at given displacement from pos — τ-safe *)
$findNeighborByDir[pos_, dir_, state_, n_] :=
  SelectFirst[state, Function[p,
    Mod[p[[1,1]] - pos[[1]], n] === dir[[1]] &&
    Mod[p[[1,2]] - pos[[2]], n] === dir[[2]]], None]
```

These use differences (`pos_i - pos_j`) so that τ cancels before `Mod` is applied.

---

## Timings (nGrid=3, N=3 distinguishable particles, Apple M-class CPU)

| Step | Time | What it does |
|------|------|-------------|
| Kernel startup | ~2–3s | wolframscript JVM + Mathematica kernel init |
| State enumeration + orbit compute | ~0.05s | Enumerate 504 states, compute 8 orbit reps |
| τ-BFS from orbit reps | ~2s | τ-freedom check at 8 states (not 504) |
| BFS from 8 orbit reps | ~2s | Exact symbolic path tree for each rep |
| D4 testing | ~4s | 16 extra BFS runs, numerical comparison |
| Ergodicity check | ~2s | Reachability BFS on derived transition graph |
| SZ direct check | ~5s | G-table precompute + 30 numerical evaluations |
| **Total** | **~15s** | |

The key speedup vs the previous version is the τ-BFS optimisation: running from 8 orbit representatives instead of all 504 states reduces that step from ~135 s to ~2 s (63× faster).

---

## Architecture overview

```
check.wls
  │
  ├─ loads dbc_core.wl (core library)
  ├─ loads algorithm file (e.g., vmmc_2d.wl)
  │
  ├─ Step 1: $dbcEnumerateNParticleStates + CheckErgodicity (count)
  ├─ Step 2: $dbcAllGElems + $dbcComputeOrbits → orbit reps
  ├─ Step 4: $dbcCheckTranslational (τ-BFS from orbit reps)
  │            └─ $dbcBuildStateLeaves × |orbit reps|
  ├─ Step 5: $dbcBuildStateLeaves × |orbit reps| → repLeaves
  ├─ Step 6: $dbcVerifyD4 (SZ-numerical comparison)
  │            └─ 2 extra BFS runs per rep
  ├─ Step 7: $dbcCheckErgodicityFromLeaves
  │            └─ BFS on derived transition graph
  └─ Step 8: $dbcSZCheckLeaves (numerical SZ, no symbolic matrix)
               ├─ Build G-action integer table (once)
               ├─ For each of 30 SZ points:
               │    ├─ Evaluate leaf weights numerically
               │    ├─ Assemble SparseArray T matrix via G-table
               │    └─ Check T(s→t)π(s) = T(t→s)π(t) for all pairs
               └─ Report violations or PASS
```

---

## Outstanding issues and known limitations

### 1. Numerical mode requires concrete coupling values

The `-mode Numerical` path runs an actual MCMC simulation to check that sample frequencies match the Boltzmann distribution. This requires `energy[state]` to evaluate to a number. For `vmmc_2d.wl` (which keeps couplings symbolic), the numerical mode reports `Indeterminate`. Set specific values for all `couplingJ[a,b,d]` atoms before using this mode.

### 2. SZ check is probabilistic, not exact

With 30 evaluation points, the probability of a false negative (declaring a wrong algorithm correct) is extremely small in practice. Use `-mode FullSimplify` for an exact certificate (warning: much slower for complex algorithms).

### 3. Only the ergodic component from the seed state is checked

The checker verifies detailed balance for all states reachable from `$seedState`. If the algorithm's state space has multiple disconnected ergodic components (e.g., different particle number sectors), each component must be checked separately with its own seed state.

### 4. D4 symmetry is only tested, not assumed

The checker always runs the D4 generator tests regardless of `$symmetryGroup`. If the algorithm fails these tests, a warning is shown, but the DB check still runs on the orbit rep leaves (which may give wrong results if D4 is genuinely broken — a future improvement would fall back to a full BFS in this case).

### 5. maxDepth paths are dropped with a warning

If an execution path requires more than `maxDepth` (default 22) random bits, it is silently dropped and a warning is printed. This can happen for algorithms that draw many random numbers per step. Increase `-maxDepth` if warnings appear.

### 6. Translation failure recomputes orbits (minor cost)

If the τ-BFS finds a violation, orbit representatives are recomputed without the translation subgroup before proceeding. The D4-only orbit computation has more representatives (up to nGrid² × 8 / 8 = nGrid²) which may increase subsequent BFS cost.

---

## Changes from v1 (initial release)

| Issue | Fix |
|-------|-----|
| τ-BFS ran all 504 states (~135s) | Now runs from 8 orbit reps only (~2s, 63× faster) |
| Coupling atom asymmetry bug (31% of leaves silently zeroed, SZ result random-seed dependent) | Added canonical symmetry rule: `couplingJ[b,a,d]→couplingJ[a,b,d]` for b>a, installed before BFS |
| Ergodicity check always passed (only compared combinatorial counts) | New reachability BFS on derived transition graph |
| D4 check used structural equality `===` (fragile) | Replaced with SZ-style numerical comparison |
| `ClearAll` on expression strings was a no-op | Now correctly clears the coupling head symbol |
| Dropped paths at maxDepth were silent | Now prints a warning with path count |

---

## Comparison with previous version

| Feature | v1 | v2 (current) |
|---------|-----|--------------|
| τ-BFS cost | ~135s (all 504 states) | ~2s (8 orbit reps) |
| Total runtime | ~155s | ~15s |
| Coupling symmetry | Bug: a>b atoms silently zeroed | Fixed: symmetry rule installed |
| Ergodicity check | Always passed (count only) | Real reachability BFS |
| D4 comparison | Structural `===` (fragile) | Numerical SZ (robust) |
| `ClearAll` | No-op on strings | Correctly clears symbol |
| maxDepth drops | Silent | Warning printed |
