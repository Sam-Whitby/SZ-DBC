# Critical audit — Julia-native checker (`julia_poc/` ≡ DB_Julia)

An independent, deliberately critical look at the Julia detailed-balance
checker, and a reconciliation with the two design notes (`julia_bfs.md`,
`optimisation.md`) written before it existed. Written 2026-06-04.

---

## 1. What the prior design notes got right, wrong, or now-outdated

### `julia_bfs.md` (τ-BFS feasibility PoC)

- **Right:** the core thesis — Julia can replicate the τ-BFS exactly with
  `TauPos` (→ `TauNum`), a `BitSeqRNG`, and `DeltaE` vectors, with no symbolic
  algebra library — is borne out. The "always-read policy", the τ=0 leaf
  substitution, and the τ-violation-via-explicit-tracking all made it into the
  production code.
- **Understated the threshold problem.** The PoC only ever handled
  Metropolis-style *clamped* thresholds `min(1, exp(−βΔE))`. VMMC needs far
  more: a **ratio of `(1−exp)` sums capped by `Min`**, which does not reduce to
  a sum of exponentials until two acceptance factors are multiplied. None of
  this is anticipated in `julia_bfs.md`. It required a whole extra component —
  the exact rational-function engine (`Val = num / ∏(1−exp) binomials`) with
  denominator clearing — that the note's "DeltaE-vector storage" sketch does
  not cover. This was the single largest piece of unforeseen work.
- **Contains a latent bug it did not catch.** The PoC's `rand_choice!` used
  `k = ceil(log2(n + 0.5))`, which is **wrong for powers of two** (`n=8` →
  `k=4` instead of 3) and would have mis-weighted `broken_8way_hop`. The
  production code uses `nbits(n) = ndigits(n-1; base=2)` (Mathematica's
  `IntegerLength[n-1,2]`). Worth flagging that the PoC "passed" only because it
  never exercised a power-of-two pool.
- **Effort estimate** ("vmmc_2d: 2–4 days") was in the right ballpark.

### `optimisation.md` (Mathematica + Julia Phase 2/3)

- **Accurate** for the architecture it describes (compact-JSON export +
  `ebe.jl` subprocess, non-strict BFS + degenerate filter). The Julia-native
  checker reproduces that Phase-2 logic faithfully and finds the **same 216
  genuine chambers** for VMMC and **48** for single-Metropolis.
- **Now falsified:** "The Mathematica floor (τ-BFS + ergodicity ≈ 16 s) cannot
  be reduced by Julia." The Julia-native τ-BFS *with translation-orbit
  reduction* does the whole VMMC τ-BFS over the 56 orbit reps in ~16 s and the
  rest of the pipeline in ~10 s — i.e. the "floor" it called irreducible is now
  the entire runtime, and is itself reducible further (see §3).
- **Superseded:** the JSON round-trip, the `feasible_sigmas` export, the
  fallback-to-Mathematica path — all unnecessary in the native pipeline. The
  "PackageCompiler sysimage → 0.3 s" observation still stands as the obvious way
  to remove JIT warmup.

---

## 2. Soundness — where the checker is exact, and where it is not

The pipeline is exact (`Rational{BigInt}` throughout) **except** at the points
below. None causes a wrong verdict on the nine examples, but an honest user
should know the boundaries.

1. **τ-detection is sound only modulo the translation discipline.** Mathematica
   is *genuinely* symbolic: τ propagates through every expression, so any
   τ-dependence in any weight is detected automatically. The Julia checker
   detects τ-dependence only for values that flow through `TauNum` arithmetic
   and the provided geometry helpers. A translation that extracts `tau0(p.r)`
   early and branches on the raw integer would hide that dependence → a **false
   τ-PASS**. This is documented as a contract ("positions stay `TauNum`"), and
   it does not affect the *detailed-balance* verdict (which is computed at τ=0
   regardless), but it is a real weakening relative to Mathematica and is the
   checker's most important caveat.

2. **The `Min`-condition derivation makes an unasserted positivity assumption.**
   `min_condition` reduces `Min[a,b]`'s switch to a linear hyperplane by
   clearing the operands' denominators — which is only direction-preserving if
   those `(1−exp(−βL))` denominators are **positive** in the region where the
   `Min` is reached. For VMMC this holds (the `Min` sits under the guard
   `eInit<eFwd`, i.e. `L>0`), but the code **checks the balanced-binomial
   structure, not the positivity**. A `Min` outside the VMMC pattern could in
   principle get a mis-oriented hyperplane. This is bounded (only `ThMin`, only
   with denominators) and should either be asserted or documented as a hard
   limit; right now it is an implicit assumption.

3. **Chamber enumeration uses a floating-point LP (HiGHS, `eps=1e-6`).** This is
   the *only* non-exact step. The initial witness `J*_a = 100^a` is exact, and
   the hyperplanes have small integer coefficients, so mis-classifying a chamber
   would require a feasibility error on a chamber thinner than `eps` — not
   possible for these arrangements, but a theoretical gap. A rational LP /
   Fourier–Motzkin step would close it at a performance cost.

4. **Open chambers only.** Detailed balance is verified on the open chambers of
   the hyperplane arrangement; measure-zero boundaries (e.g. `J1=J2` exactly)
   are not separately checked. Sound because transition probabilities are
   continuous in the couplings, so any positive-measure violation appears in an
   adjacent open chamber. Inherited from the Mathematica design; correct but
   worth stating.

5. **Supported exact-weight class is bounded.** Weights must be rational
   functions whose denominators are products of `(1−exp(−βL))` binomials.
   Anything else (a different threshold algebra) is a **hard error**, not a
   guess — which is the right failure mode, but means "faithful translation" is
   only available for algorithms in this class.

---

## 3. Performance — honest state and the remaining headroom

VMMC went from **355 s → ~26 s** via four profiling-driven changes (per-pair
condition projection; lazy weight caching + unique-weight scan; translation-
orbit reduction; threshold hash-consing). The simple examples run in 6–11 s,
**most of which is Julia JIT warmup** (~4–5 s per process).

Remaining, *not* done (correctness was prioritised over these):

- **Lattice coordinates are carried as `Rational{BigInt}`.** Positions and
  distances are tiny integers; the residual VMMC bottleneck is the symbolic BFS,
  dominated by this. A narrower type (`Rational{Int64}`, or an integer position
  with a small τ-tracker) is the obvious next win — deferred only because it
  touches the foundational numeric type and a silent overflow there would be a
  wrong verdict, so it needs overflow-checked arithmetic to stay sound.
- **No memoisation of structurally-identical cluster sub-trees** in the VMMC
  BFS; the same partial cluster is re-expanded along different paths.
- **No PackageCompiler sysimage**, so every CLI invocation pays JIT. This is the
  difference between "≈10 s" and "≈1 s" for the simple examples.

The user's suggested "merge Phase 1 and Phase 2 (resolve inequalities during the
BFS)" was considered and **rejected on profiling grounds**: it would re-run the
BFS once per chamber (216× for VMMC), whereas the chosen design BFSs once and
evaluates thresholds per chamber lazily — strictly less work given that control
flow here is coupling-independent (only the *weights*, never the branch
structure, depend on the couplings).

---

## 4. Smaller criticisms

- **Determinism:** confirmed — no RNG anywhere; BFS over fixed bit trees;
  canonical `sort` on atoms/conditions; HiGHS LP is deterministic. Dict
  iteration order is used only to build order-independent sets. Good.
- **Hash-consing relies on `objectid` stability**, which holds because
  `_TH_CACHE` keeps interned nodes alive for the run and is reset at the start
  of `build_transitions`. Calling the `th_*` builders outside that lifecycle
  would be unsafe; this is an internal invariant, not part of the public API,
  but is undocumented in the code.
- **`enumerate_states` is O(P(S,N))** (all k-permutations) — fine for the small
  systems these checks target, exponential beyond them. The whole approach is
  exhaustive by design and does not scale to large lattices; it is a *verifier*
  for small representative systems, not a production sampler check.
- **Error messages are good** (`CantHandle` carries a specific reason), but a
  few internal `cant(...)` sites ("bad binomial numerator") would be opaque to a
  user translating a new algorithm.
- The SZ-DBC repo now carries **two engines** (`dbc_core.wl` Mathematica +
  `julia_poc/dbc.jl`); they are independent and must be kept in sync by hand if
  the semantics change. DB_Julia is the single-engine version.
