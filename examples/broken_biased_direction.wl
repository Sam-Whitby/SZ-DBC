(* ================================================================
   broken_biased_direction.wl  —  Single-particle move with biased direction
   ================================================================
   Identical to single_metropolis.wl except that the displacement
   {0,1} ("right") appears TWICE in the direction list, so it is
   chosen with probability 2/9 while every other direction has 1/9.

   Why this breaks detailed balance:
     For any pair (s, t) where t = s + {0,1} (move right):
       T(s→t)  uses proposal probability 2/9
       T(t→s)  uses proposal probability 1/9  (the reverse is "left")
     The Metropolis accept/reject correctly gives π(t)/π(s), so:
       T(s→t)/T(t→s) = 2 × π(t)/π(s)   ≠   π(t)/π(s)
     Detailed balance is violated by a constant factor of 2 for
     every pair reachable by a rightward move.

   τ-check passes because the direction list is fixed (no τ) and
   the acceptance probability uses τ-free pairwise ΔE.

   Expected result:  τ PASS,  DB FAIL
   ================================================================ *)


$nGrid         = 3;
numBeta        = 1.0;
$maxD2         = 2;

$particleTypes = {1, 2, 3};

$seedState = Module[
  {pos = Flatten[Table[{r, c}, {r, $nGrid}, {c, $nGrid}], 1],
   types = Sort[$particleTypes]},
  SortBy[MapThread[{#1, #2} &, {Take[pos, Length[types]], types}], First]];

$symmetryGroup = {"translation"};


(* ---- Geometry ---- *)

(* {0,1} appears twice → P(right) = 2/9, all others P = 1/9.
   This is the only change from single_metropolis.wl. *)
$biasedDisps = {{0,1},{0,1},{0,-1},{1,0},{-1,0},{1,1},{1,-1},{-1,1},{-1,-1}};

$bdPairD2[p1_, p2_, n_] :=
  With[{dr = p1[[1]] - p2[[1]], dc = p1[[2]] - p2[[2]]},
    With[{dra = Mod[Abs[dr], n], dca = Mod[Abs[dc], n]},
      Min[dra, n - dra]^2 + Min[dca, n - dca]^2]]

$bdPairEnergy[state_, n_] :=
  Module[{total = 0},
    Do[
      With[{d2 = $bdPairD2[state[[i, 1]], state[[j, 1]], n]},
        If[0 < d2 <= $maxD2,
          total += couplingJ[state[[i, 2]], state[[j, 2]], d2]]],
      {i, Length[state]}, {j, i + 1, Length[state]}];
    total]

energy[state_] := $bdPairEnergy[state, $nGrid]


(* ---- Algorithm ---- *)

Algorithm[state_List] :=
  Module[{n = $nGrid, particle, dir, newPos, rest},
    particle = RandomChoice[state];
    dir      = RandomChoice[$biasedDisps];    (* biased: right appears twice *)
    newPos   = particle[[1]] + dir;
    rest     = DeleteCases[state, particle, 1];
    If[AnyTrue[rest, Function[p,
        Mod[p[[1, 1]] - newPos[[1]], n] === 0 &&
        Mod[p[[1, 2]] - newPos[[2]], n] === 0]],
      Return[state]];
    With[{newState = SortBy[Append[rest, {newPos, particle[[2]]}], First]},
      With[{dE = energy[newState] - energy[state]},
        If[RandomReal[] < Piecewise[{{1, dE <= 0}}, Exp[-\[Beta] * dE]],
          newState, state]]]]


(* ---- Symbolic parameters ---- *)

DynamicSymParams[states_List] :=
  Module[{types, n, d2Vals, atoms},
    types  = Sort @ DeleteDuplicates @ Flatten[#[[2]] & /@ # & /@ states, 1];
    n      = $nGrid;
    d2Vals = Sort @ DeleteDuplicates @ Select[
      Flatten @ Table[Min[dr, n - dr]^2 + Min[dc, n - dc]^2,
        {dr, 0, Floor[n/2]}, {dc, 0, Floor[n/2]}],
      0 < # <= $maxD2 &];
    atoms  = Flatten @ Table[
      Table[couplingJ[a, b, d2], {d2, d2Vals}],
      {a, types}, {b, a, Max[types]}];
    <|"couplings" -> atoms, "numericParams" -> {}|>]
