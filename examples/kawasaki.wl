(* ================================================================
   kawasaki.wl  —  Kawasaki exchange dynamics on a 2D periodic lattice
   ================================================================
   Each step:
     1. Form the list of all pairs of particles with distinct types.
     2. Choose one pair uniformly at random.
     3. Swap the types at the two sites (positions unchanged).
     4. Accept with probability min(1, exp(-β ΔE)).

   Because positions never change, pairwise distances are invariant.
   ΔE depends only on how the type assignment at each pair of sites
   changes — purely symbolic in couplingJ, τ-free.

   The proposal is symmetric: |distinct-type pairs| is the same in
   any state with fixed type multiset, so P(propose s→t) = P(propose t→s).

   Kawasaki exchange never moves particles — only swaps types between
   fixed sites.  Starting from any seed, only the N! = 3! = 6 states
   that are type-permutations of the seed's positions are reachable.
   Ergodicity over the full 504-state position+type space therefore
   FAILS by design; this is correct physics, not a checker bug.

   Expected result:  τ PASS,  DB PASS,  ERGODICITY FAIL (by design)
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


(* ---- Geometry / energy ---- *)

$kwPairD2[p1_, p2_, n_] :=
  With[{dr = p1[[1]] - p2[[1]], dc = p1[[2]] - p2[[2]]},
    With[{dra = Mod[Abs[dr], n], dca = Mod[Abs[dc], n]},
      Min[dra, n - dra]^2 + Min[dca, n - dca]^2]]

$kwPairEnergy[state_, n_] :=
  Module[{total = 0},
    Do[
      With[{d2 = $kwPairD2[state[[i, 1]], state[[j, 1]], n]},
        If[0 < d2 <= $maxD2,
          total += couplingJ[state[[i, 2]], state[[j, 2]], d2]]],
      {i, Length[state]}, {j, i + 1, Length[state]}];
    total]

energy[state_] := $kwPairEnergy[state, $nGrid]


(* ---- Algorithm ---- *)

Algorithm[state_List] :=
  Module[{n = $nGrid, pairs, pair, p1, p2, newState, dE},

    (* All pairs (i,j) with i<j and distinct types *)
    pairs = Select[
      Flatten[Table[{state[[i]], state[[j]]},
               {i, Length[state]}, {j, i + 1, Length[state]}], 1],
      #[[1, 2]] =!= #[[2, 2]] &];

    If[pairs === {}, Return[state]];

    pair = RandomChoice[pairs];
    p1 = pair[[1]];  p2 = pair[[2]];

    (* Swap types; positions stay fixed *)
    newState = SortBy[
      Join[Select[state, # =!= p1 && # =!= p2 &],
           {{p1[[1]], p2[[2]]}, {p2[[1]], p1[[2]]}}],
      First];

    dE = energy[newState] - energy[state];
    If[RandomReal[] < Piecewise[{{1, dE <= 0}}, Exp[-\[Beta] * dE]],
      newState, state]]


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
