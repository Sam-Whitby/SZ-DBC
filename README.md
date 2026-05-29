# SZ-DBC: Detailed Balance Checker for Lattice MCMC

SZ-DBC automatically verifies that an MCMC algorithm satisfies **detailed balance** — the condition required for the algorithm to correctly sample the Boltzmann distribution. It is designed for lattice Monte Carlo algorithms where particles sit on a grid.

The checker produces algebraic or probabilistic proofs of correctness, not just numerical tests. It can verify translational and rotational symmetry exactly, and checks detailed balance across the entire reachable state space in a single run.

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

### Step 1 — Translational invariance (algebraic, τ-BFS)

The checker augments every particle position with a symbolic offset `{τr, τc}`:

```
{{1+τr, 1+τc}, 1}, {{1+τr, 2+τc}, 2}, ...
```

It then runs the algorithm on this symbolic state and follows all possible execution paths. Because all spatial computations in a correctly written algorithm use **pairwise differences** (e.g., `pos_i - pos_j`), the symbol τ cancels algebraically before it ever enters a distance, energy, or acceptance probability:

```
(pos_i + τ) - (pos_j + τ) = pos_i - pos_j
```

After the BFS completes, the checker scans every transition probability for τ using `FreeQ`. If τ never appears, translational invariance is **algebraically certified** for the entire reachable state space in a single BFS run. This is not a numerical test — it is exact.

**Requirement on the algorithm**: the algorithm must compute periodic boundary conditions by wrapping *differences*, not absolute positions. `Mod[pos_i - pos_j, n]` is correct; `Mod[pos_i, n] - Mod[pos_j, n]` is not (and would also be incorrect physics for the minimum-image convention).

### Step 2 — State enumeration and ergodicity

The checker enumerates all `N!/(n₁! n₂! ...)` distinct states for the given particle type multiset on the `nGrid × nGrid` torus, where `n₁, n₂, ...` are the multiplicities of each type. It checks that the BFS discovers exactly this many states (ergodicity check). If the count is wrong, the algorithm is not ergodic for this component.

### Step 3 — G-orbit computation and BFS from representatives

For a system with symmetry group G (translations × D4), many states are related by symmetry. The checker computes the G-orbits of the state space and identifies one **orbit representative** per orbit. It then runs the full BFS — tracking all possible execution paths and their exact symbolic probability weights — from each representative only. On a `3×3` grid with 3 distinguishable particles, there are 504 states but only 8 orbit representatives (a 63× reduction).

The BFS works by **intercepting all random number calls** made by the algorithm (via Mathematica's `Block` mechanism) and replacing them with symbolic interval-tracking tokens. At each random number comparison `u < p`, the BFS splits into two branches — one where the condition is true, one where it is false — accumulating the probability of each branch as an exact symbolic weight. The result is a complete tree of all possible execution paths and their transition probabilities.

### Step 4 — D4 rotational invariance (generator testing)

For each orbit representative `s₀`, the checker runs two additional BFS runs:
- One from `rot90(s₀)` — the state rotated 90° clockwise
- One from `reflect(s₀)` — the state reflected across the main diagonal

It then checks, for each transition `T(s₀ → t)`, whether:

```
T(rot90(s₀) → rot90(t)) = T(s₀ → t)
T(reflect(s₀) → reflect(t)) = T(s₀ → t)
```

By group theory, checking these two generators (rotation and reflection) is sufficient to certify all 8 elements of D4. This requires only 2 extra BFS runs per representative — 16 total for 8 representatives — rather than checking all 8 × 8 = 64 combinations.

### Step 5 — Detailed balance check (SZPure, default)

The checker uses a **Schwartz-Zippel probabilistic zero-test**. Rather than building a large symbolic transition matrix and applying `FullSimplify` (which can take hours), it evaluates the detailed balance condition numerically at 30 randomly chosen parameter points `{β, J₁, J₂, ...}`.

For each evaluation point:
1. The symbolic leaf weights (products of acceptance probabilities) evaluate immediately to floating-point numbers, because all Piecewise conditions (e.g., `J₁ < J₂`) become concrete comparisons once coupling values are substituted.
2. The full transition matrix for all 504 states is assembled via the precomputed G-action integer table — no symbolic operations needed.
3. The detailed balance condition `T(s→t) · exp(-β E(s)) = T(t→s) · exp(-β E(t))` is checked numerically for every pair.

If the algorithm satisfies detailed balance, this expression is identically zero as a function of the parameters, and all 30 evaluations return zero. If there is a violation, it appears as a nonzero polynomial in the parameters, which a random evaluation detects with high probability (by the Schwartz-Zippel lemma). With 30 evaluations and coupling constants drawn from a bounded rational range, the probability of a false negative is negligibly small.

**Why this is faster than FullSimplify**: The transition probability from one state to another is a symbolic Piecewise expression — a product of acceptance weights, each of which is a Piecewise function of the coupling constants (e.g., `1 - exp(β(J₁ - J₂))` if `J₁ < J₂`, else 0). Summing these products over many execution paths produces an expression that is symbolically enormous. `FullSimplify` on such expressions can take hours. Numerical substitution makes all Piecewise conditions collapse immediately to True/False, reducing the weight to a simple float product. This is the key architectural insight.

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

Numerical MCMC check only (fast sanity check, not a proof):

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

Step 1: Translational invariance (τ-BFS)
  PASS  — translation invariance algebraically verified

Step 2: State enumeration
  States found : 504  (0.01s)
  Ergodicity   : PASS  (found 504, theoretical 504)

Step 5: G-orbits & BFS
  |G| = 72
  Orbit reps : 8
  BFS done: 8 reps  (2.19s)

Step 6: D4 generator testing
  PASS  — D4 invariance verified (8 reps × 2 generators)
  Time: 4.48s

Step 7: Detailed balance (SZPure direct)
  Time: 8.77s
  PASS  — detailed balance satisfied for all 504 states

================================================================
  SUMMARY
================================================================
  Translational : PASS (algebraic)
  D4 rotational : PASS
  Ergodicity    : PASS (504/504 states)
  Detailed bal. : PASS
================================================================
```

Total runtime: approximately 15–20 seconds on a modern laptop. Exit code 0 = all checks passed; exit code 1 = failure.

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

$checkerAbstractParams = {"physLen", "epsLJ"};  (* params kept symbolic *)

Algorithm[state_List] := ...   (* MCMC step *)
energy[state_]        := ...   (* bare energy (no β) *)

DynamicSymParams[states_List] := ...  (* coupling atom list *)
```

### State format rules

- `state` is a sorted list of `{{{row,col}, type}, ...}` pairs
- Positions are integers in `{1, ..., nGrid}` (1-indexed)
- The checker applies PBC normalisation to returned states — you do **not** need to `Mod` positions in the algorithm
- All spatial operations must use **pairwise differences** (so that τ cancels). Never compute `Mod[absolute_position, nGrid]` directly on a stored coordinate

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
- The energy function `energy[state]` must be a symbolic expression in the coupling atoms listed by `DynamicSymParams` and β is **not** included (it is multiplied in by the checker)
- `DynamicSymParams[states]` should return `<|"couplings" -> {couplingJ[1,1,1], couplingJ[1,2,1], ...}, "numericParams" -> {}|>` listing every coupling atom that appears in the energy

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
| τ-BFS (translation check) | ~2–4s | State-space BFS with symbolic τ offset |
| State enumeration | <0.1s | List all 504 states |
| BFS from 8 orbit reps | ~2s | Exact symbolic path tree for each rep |
| D4 testing | ~4.5s | 16 extra BFS runs (2 generators × 8 reps) |
| SZ direct check | ~9s | G-table precompute + 30 numerical evaluations |
| **Total** | **~15–20s** | |

The 72× orbit speedup (D4 + translations: 8 reps instead of 504) reduces the BFS work from ~2 minutes to ~2 seconds. The SZ numerical path avoids the FullSimplify symbolic bottleneck (which took >11 minutes).

---

## Outstanding issues and known limitations

### 1. Runtime is ~15–20s, not 1s

The dominant costs are:
- **BFS**: ~0.25s per orbit representative. For 8 reps this is 2s. The BFS is inherently serial (each bit sequence is run sequentially to build the tree).
- **D4 testing**: 16 extra BFS runs (~4.5s). This is the biggest single cost and is proportional to |G| × 2/|G|_D4 = 2 × (number of reps).
- **G-action table precomputation**: inside the SZ check, 36K calls to `$dbcApplyGElem` to build the integer action table once before the 30 SZ evaluations.
- **Kernel startup**: wolframscript always takes 2–3s to initialise Mathematica, regardless of what the script does.

Getting to 1s would require either running inside a persistent kernel (removing startup cost) or further optimising the BFS itself (possibly via compiled Mathematica functions or reduced bit-depth).

### 2. SZ check is probabilistic, not exact

With 30 evaluation points, the probability of a false negative (declaring a wrong algorithm correct) is extremely small in practice — each violation that exists in an open parameter region is detected with probability at least `1 - (max_degree / range)^30` per Schwartz-Zippel. However, it is not an exact proof in the way FullSimplify is. Use `-mode FullSimplify` for an exact certificate (warning: much slower for complex algorithms).

### 3. Only the ergodic component from the seed state is checked

The checker verifies detailed balance for all states reachable from `$seedState`. If the algorithm's state space has multiple disconnected ergodic components (e.g., different particle number sectors), each component must be checked separately with its own seed state. The ergodicity count tells you whether the full expected component was explored.

### 4. D4 symmetry is only tested, not assumed

The checker always runs the D4 generator tests regardless of `$symmetryGroup`. If the algorithm fails these tests, the speedup is not applied and a warning is shown, but the DB check still runs on the orbit rep leaves (which may give wrong results if D4 is genuinely broken — a future improvement would be to fall back to a full BFS in this case).

### 5. τ-BFS verifies translation invariance algebraically but may be slow for large systems

The τ-BFS runs a full state-space BFS (all 504 states for the current example). For larger grids or more particles, this could become expensive. An optimisation would be to run the τ-BFS only from the orbit representatives (sufficient if D4 invariance has already been confirmed, since D4 images of the reps cover the full space).

### 6. Algorithm must not branch on absolute position

If the algorithm contains a conditional like `If[row > nGrid/2, ...]` where `row` includes the τ offset, Mathematica cannot evaluate the condition symbolically and the τ-BFS will stall or throw a `$dbc$cantHandle` error. This is correctly detected and reported — it means the algorithm depends on absolute position and is not translation-invariant, which is the correct diagnosis.

### 7. Coupling atoms must be declared correctly in DynamicSymParams

The SZ check substitutes random values for all coupling atoms listed by `DynamicSymParams`. If an atom is used in the energy or acceptance weights but not listed, it will remain symbolic during the SZ evaluation and the numerical check will fail or give incorrect results. Make sure every `couplingJ[...]` symbol that appears in the energy is listed.

---

## Architecture overview

```
check.wls
  │
  ├─ loads dbc_core.wl (core library)
  ├─ loads algorithm file (e.g., vmmc_2d.wl)
  │
  ├─ Step 1: $dbcCheckTranslational   → τ-BFS, FreeQ check
  ├─ Step 2: $dbcEnumerateNParticleStates + CheckErgodicity
  ├─ Step 5: $dbcComputeOrbits + $dbcBuildStateLeaves (BFS per rep)
  │            └─ BuildTreeAT → RunWithBitsAT (intercepts RandomReal etc.)
  ├─ Step 6: $dbcVerifyD4 (2 BFS runs per rep)
  └─ Step 7: $dbcSZCheckLeaves (numerical SZ, no symbolic matrix)
               ├─ Build G-action integer table (once)
               ├─ For each of 30 SZ points:
               │    ├─ Evaluate leaf weights numerically
               │    ├─ Assemble SparseArray T matrix via G-table
               │    └─ Check T(s→t)π(s) = T(t→s)π(t) for all pairs
               └─ Report violations or PASS
```

---

## Comparison with previous version (Check_Detailed_Balance)

| Feature | Old (Check_Detailed_Balance) | New (SZ-DBC) |
|---------|------------------------------|--------------|
| State format | Flat occupancy vector | `{position, type}` pairs |
| Default check | FullSimplify (exact, slow) | SZPure numerical (fast) |
| Translation invariance | Assumed/declared | Algebraically verified (τ-BFS) |
| D4 invariance | Assumed/declared | Generator-tested (exact) |
| Parallelisation | Multi-kernel (expensive startup) | None needed (fast enough) |
| SZ bottleneck | Symbolic Piecewise explosion | Eliminated (leaf-level eval) |
| Runtime (3×3, N=3) | >11 min (SZ path broken) | ~15–20s |
