(* ================================================================
   single_metropolis.wl  —  Single-particle Metropolis on a 2D periodic lattice
   ================================================================
   Each step:
     1. Choose one particle uniformly at random.
     2. Choose a displacement uniformly from the 8 nearest/next-nearest
        directions (|dx|,|dy| <= 1, not both zero).
     3. Reject immediately if the target site is occupied.
     4. Accept the move with probability min(1, exp(-β ΔE)).

   Proposal is symmetric in forward/reverse: both directions are
   equiprobable, and occupancy rejection is applied equally.
   Energy uses pairwise couplingJ[type_a, type_b, d²].

   Expected result:  τ PASS,  DB PASS
   ================================================================ *)


$nGrid         = 3;
numBeta        = 1.0;
$maxD2         = 2;

$particleTypes = {1, 2, 3};

$seedState = Module[
  {pos = Flatten[Table[{r, c}, {r, $nGrid}, {c, $nGrid}], 1],
   types = Sort[$particleTypes]},
  SortBy[MapThread[{#1, #2} &, {Take[pos, Length[types]], types}], First]];

$symmetryGroup = {"translation", "D4"};   (* D4 verified numerically in Step 5b *)


(* ---- Geometry ---- *)

$smDisps = DeleteCases[
  Flatten[Table[{dx, dy}, {dx, -1, 1}, {dy, -1, 1}], 1],
  {0, 0}];

$smPairD2[p1_, p2_, n_] :=
  With[{dr = p1[[1]] - p2[[1]], dc = p1[[2]] - p2[[2]]},
    With[{dra = Mod[Abs[dr], n], dca = Mod[Abs[dc], n]},
      Min[dra, n - dra]^2 + Min[dca, n - dca]^2]]

$smPairEnergy[state_, n_] :=
  Module[{total = 0},
    Do[
      With[{d2 = $smPairD2[state[[i, 1]], state[[j, 1]], n]},
        If[0 < d2 <= $maxD2,
          total += couplingJ[state[[i, 2]], state[[j, 2]], d2]]],
      {i, Length[state]}, {j, i + 1, Length[state]}];
    total]

energy[state_] := $smPairEnergy[state, $nGrid]


(* ---- Algorithm ---- *)

Algorithm[state_List] :=
  Module[{n = $nGrid, particle, dir, newPos, rest},
    particle = RandomChoice[state];
    dir      = RandomChoice[$smDisps];
    newPos   = particle[[1]] + dir;          (* no Mod — checker normalises *)
    rest     = DeleteCases[state, particle, 1];

    (* Hard-core rejection: new position occupied?
       Use Mod of difference so τ cancels algebraically. *)
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
