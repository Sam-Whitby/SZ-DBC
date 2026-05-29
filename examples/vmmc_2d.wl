(* ================================================================
   vmmc_2d.wl — VMMC on a 2D square lattice (checker-compatible)
   ================================================================
   State format: sorted list of {{row, col}, type} pairs.
     row, col ∈ {1,..,nGrid},  type ∈ {1, 2, ...}.
   All spatial queries use PAIRWISE DIFFERENCES so the checker
   can verify translational invariance algebraically (τ approach).

   Required checker globals (set below or by user):
     $nGrid          — lattice side length (integer)
     numBeta         — inverse temperature (numeric)
     $maxD2          — max squared interaction distance
     $particleTypes  — type multiset for this check run
     $seedState      — canonical seed state for the checker

   Checker reads:
     Algorithm[state_]        — MCMC step
     energy[state_]           — bare energy (no β)
     $symmetryGroup           — {"translation","D4"} for speedup
     DynamicSymParams[states] — symbolic coupling atoms per component
   ================================================================ *)


(* ================================================================
   SYSTEM PARAMETERS — edit these for your system
   ================================================================ *)

$nGrid         = 3;              (* lattice side length *)
numBeta        = 1.0;            (* inverse temperature *)
$maxD2         = 2;              (* include nearest and next-nearest neighbours *)
physLen        = 1.;             (* physical diameter (sets displacement range) *)
epsLJ          = 1.;             (* LJ energy scale *)

$checkerAbstractParams = {"physLen", "epsLJ"};

$particleTypes  = {1, 2, 3};    (* type multiset for this component check *)

(* Seed state: first $particleTypes placed at lexicographic-first positions *)
$seedState = Module[{pos = Flatten[Table[{r,c},{r,$nGrid},{c,$nGrid}],1], types = Sort[$particleTypes]},
  SortBy[MapThread[{#1,#2}&, {Take[pos, Length[types]], types}], First]];

$symmetryGroup = {"translation", "D4"};

$checkerAbstractParams = {"physLen", "epsLJ"};


(* ================================================================
   GEOMETRY HELPERS
   ================================================================
   All spatial operations use PAIRWISE DIFFERENCES (not absolute
   positions wrapped by Mod) so they are compatible with the
   checker's symbolic τ-offset for translation-invariance proofs. *)

(* Minimum-image squared distance between positions p1, p2 on nGrid×nGrid torus.
   Uses the difference p1-p2 (τ cancels if positions include a τ offset). *)
$pairD2[p1_, p2_, n_] :=
  With[{dr0 = p1[[1]] - p2[[1]], dc0 = p1[[2]] - p2[[2]]},
    With[{dra = Mod[Abs[dr0], n], dca = Mod[Abs[dc0], n]},
      Min[dra, n - dra]^2 + Min[dca, n - dca]^2]]

(* Find any particle in `state` whose position differs from `pos` by `dir`
   (modulo nGrid).  Returns the full {{r,c},type} pair, or None.
   Uses Mod on the DIFFERENCE — τ cancels automatically. *)
$findNeighborByDir[pos_, dir_, state_, n_] :=
  SelectFirst[state, Function[p,
    Mod[p[[1,1]] - pos[[1]], n] === dir[[1]] &&
    Mod[p[[1,2]] - pos[[2]], n] === dir[[2]]], None]

(* All candidate neighbours of position pos within squared distance $maxD2. *)
$candidateNeighbors[pos_, state_, n_] :=
  Module[{rMax = Ceiling[Sqrt[$maxD2]]},
    Select[state, Function[p,
      With[{d2 = $pairD2[p[[1]], pos, n]},
        p[[1]] =!= pos && 0 < d2 <= $maxD2]]]]


(* ================================================================
   PAIR ENERGY (Lennard-Jones style, user-defined coupling)
   ================================================================ *)

(* couplingJ[ta, tb, d2] — coupling energy between types ta and tb at
   squared distance d2.  Left unbound during symbolic BFS.
   User may override with concrete values for testing. *)

(* $pairEnergy: compute total pairwise energy of a state *)
$pairEnergyTotal[state_, n_] :=
  Module[{pts = state, total = 0},
    Do[
      With[{d2 = $pairD2[pts[[i,1]], pts[[j,1]], n]},
        If[0 < d2 <= $maxD2,
          total += couplingJ[pts[[i,2]], pts[[j,2]], d2]]],
      {i, Length[pts]}, {j, i+1, Length[pts]}];
    total]

energy[state_] := $pairEnergyTotal[state, $nGrid]


(* ================================================================
   VMMC: Whitelam-Geissler cluster builder
   ================================================================ *)

$nStep = Max[1, Round[physLen]];
$displacements = DeleteCases[
  Flatten[Table[{dx,dy},{dx,-$nStep,$nStep},{dy,-$nStep,$nStep}],1],
  {0,0}];

$pairE[ti_, tj_, pi_, pj_, n_] :=
  With[{d2 = $pairD2[pi, pj, n]},
    If[0 < d2 <= $maxD2, couplingJ[ti, tj, d2], 0]]

$vmmcBuildCluster[state_, n_, seedParticle_, dir_] :=
  Module[{cluster = {seedParticle},
          inCluster = <|seedParticle[[1]] -> True|>,
          queue = {seedParticle},
          frustrated = False,
          p, pPos, pType, pPost, pRev,
          cands, q, qPos, qType,
          eInit, eFwd, eRev, wFwd, wRev, r1, r2},
    While[queue =!= {} && !frustrated,
      p     = First[queue]; queue = Rest[queue];
      pPos  = p[[1]];
      pType = p[[2]];
      pPost = pPos + dir;    (* no Mod — checker canonicalises *)
      pRev  = pPos - dir;

      (* Candidates: occupied neighbours within $maxD2, not in cluster *)
      cands = SortBy[
        Select[state, Function[q,
          !KeyExistsQ[inCluster, q[[1]]] &&
          (0 < $pairD2[q[[1]], pPos,  n] <= $maxD2 ||
           0 < $pairD2[q[[1]], pPost, n] <= $maxD2 ||
           0 < $pairD2[q[[1]], pRev,  n] <= $maxD2)]],
        (* Canonical ordering: sort by distance tuple so G-related states
           produce identical random-number sequences during BFS. *)
        {$pairD2[#[[1]], pPos, n] &, $pairD2[#[[1]], pPost, n] &,
         $pairD2[#[[1]], pRev,  n] &, #[[2]] &}];

      Do[
        q     = cands[[k]];
        qPos  = q[[1]];
        qType = q[[2]];

        eInit = $pairE[pType, qType, pPos,  qPos, n];
        eFwd  = $pairE[pType, qType, pPost, qPos, n];
        eRev  = $pairE[pType, qType, pRev,  qPos, n];

        wFwd = Piecewise[{{1, eFwd === Infinity}, {1 - Exp[\[Beta](eInit-eFwd)], eInit < eFwd}}, 0];
        wRev = Piecewise[{{1, eRev === Infinity}, {1 - Exp[\[Beta](eInit-eRev)], eInit < eRev}}, 0];

        r1 = RandomReal[];
        If[r1 <= wFwd,
          r2 = RandomReal[];
          If[r2 > Piecewise[{
              {1, eFwd === Infinity && eRev === Infinity},
              {1-Exp[\[Beta](eInit-eRev)], eFwd === Infinity && eInit < eRev},
              {0, eFwd === Infinity},
              {1, eRev === Infinity && eInit < eFwd},
              {Min[(1-Exp[\[Beta](eInit-eRev)])/(1-Exp[\[Beta](eInit-eFwd)]),1],
               eInit < eFwd && eInit < eRev},
              {0, eInit < eFwd}}, 0],
            frustrated = True; Break[],
            AppendTo[cluster, q];
            inCluster[qPos] = True;
            AppendTo[queue, q]]],
        {k, Length[cands]}]
    ];
    If[frustrated, None, cluster]]


(* ================================================================
   Algorithm: one VMMC step
   ================================================================
   1. Choose seed particle uniformly.
   2. Choose displacement uniformly from $displacements.
   3. Build cluster (Whitelam-Geissler).
   4. Translate cluster by dir; reject on hard-core overlap.

   NOTE: positions are NOT wrapped with Mod after displacement.
   The checker canonicalises returned states (applies PBC).
   This keeps positions in a form that allows τ-cancellation. *)

Algorithm[state_List] :=
  Module[{n = $nGrid, seedParticle, dir, cluster, clusterPos,
          newState, inCluster, destPos, overlap},
    If[state === {}, Return[state]];

    seedParticle = RandomChoice[state];
    dir          = RandomChoice[$displacements];

    cluster = $vmmcBuildCluster[state, n, seedParticle, dir];
    If[cluster === None, Return[state]];

    clusterPos = <|#[[1]] -> True & /@ cluster|>;
    newState   = Select[state, !KeyExistsQ[clusterPos, #[[1]]] &];

    (* Check hard-core overlap: for each cluster particle, its new position
       must not coincide (mod nGrid) with any non-cluster particle. *)
    overlap = False;
    Do[
      destPos = cluster[[i,1]] + dir;   (* no Mod — checker canonicalises *)
      (* Check overlap: Mod of difference cancels τ *)
      If[AnyTrue[newState, Function[p,
          Mod[p[[1,1]] - destPos[[1]], n] === 0 &&
          Mod[p[[1,2]] - destPos[[2]], n] === 0]],
        overlap = True; Break[]],
      {i, Length[cluster]}];
    If[overlap, Return[state]];

    (* Move cluster particles *)
    Do[
      AppendTo[newState, {cluster[[i,1]] + dir, cluster[[i,2]]}],
      {i, Length[cluster]}];

    (* SortBy canonical position (τ cancels: #[[1,1]]-τr gives integer) *)
    SortBy[newState, First]
  ]


(* ================================================================
   Symbolic parameters for the checker
   ================================================================ *)

DynamicSymParams[states_List] :=
  Module[{types, n, d2Vals, couplingAtoms},
    types  = Sort @ DeleteDuplicates @ Flatten[#[[2]] & /@ # & /@ states, 1];
    n      = $nGrid;
    d2Vals = Sort @ DeleteDuplicates @ Select[
      Flatten @ Table[
        Min[dr, n-dr]^2 + Min[dc, n-dc]^2,
        {dr, 0, Floor[n/2]}, {dc, 0, Floor[n/2]}],
      0 < # <= $maxD2 &];
    couplingAtoms = Flatten @ Table[
      If[a <= b, Table[couplingJ[a, b, d2], {d2, d2Vals}], Nothing],
      {a, types}, {b, types}];
    <|"couplings"     -> couplingAtoms,
      "numericParams" -> {}|>]
