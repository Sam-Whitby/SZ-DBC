(* ================================================================
   broken_metropolis_halfbeta.wl  —  Single-particle move with wrong accept
   ================================================================
   Identical to single_metropolis.wl except that the Metropolis
   acceptance probability uses Exp[-β·ΔE/2] instead of Exp[-β·ΔE].

   Why this breaks detailed balance:
     For a pair (s,t) with ΔE = E(t)-E(s) > 0:
       T(s→t) = (1/N)(1/8) · Exp[-β·ΔE/2]
       T(t→s) = (1/N)(1/8) · 1            (reverse move is downhill)
     DB requires T(s→t)/T(t→s) = π(t)/π(s) = Exp[-β·ΔE]
     But the actual ratio is Exp[-β·ΔE/2] ≠ Exp[-β·ΔE] for ΔE ≠ 0.

     This is a subtle implementation error that could arise from
     accidentally writing β/2 instead of β (e.g., confusing a factor-
     of-2 convention in the energy definition).

   τ-check passes because the acceptance probability depends only on
   the pairwise ΔE (τ-free).

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

$mhDisps = DeleteCases[
  Flatten[Table[{dx, dy}, {dx, -1, 1}, {dy, -1, 1}], 1],
  {0, 0}];

$mhPairD2[p1_, p2_, n_] :=
  With[{dr = p1[[1]] - p2[[1]], dc = p1[[2]] - p2[[2]]},
    With[{dra = Mod[Abs[dr], n], dca = Mod[Abs[dc], n]},
      Min[dra, n - dra]^2 + Min[dca, n - dca]^2]]

$mhPairEnergy[state_, n_] :=
  Module[{total = 0},
    Do[
      With[{d2 = $mhPairD2[state[[i, 1]], state[[j, 1]], n]},
        If[0 < d2 <= $maxD2,
          total += couplingJ[state[[i, 2]], state[[j, 2]], d2]]],
      {i, Length[state]}, {j, i + 1, Length[state]}];
    total]

energy[state_] := $mhPairEnergy[state, $nGrid]


(* ---- Algorithm ---- *)

Algorithm[state_List] :=
  Module[{n = $nGrid, particle, dir, newPos, rest},
    particle = RandomChoice[state];
    dir      = RandomChoice[$mhDisps];
    newPos   = particle[[1]] + dir;
    rest     = DeleteCases[state, particle, 1];
    If[AnyTrue[rest, Function[p,
        Mod[p[[1, 1]] - newPos[[1]], n] === 0 &&
        Mod[p[[1, 2]] - newPos[[2]], n] === 0]],
      Return[state]];
    With[{newState = SortBy[Append[rest, {newPos, particle[[2]]}], First]},
      With[{dE = energy[newState] - energy[state]},
        (* BUG: uses β/2 instead of β in the exponent.
           The correct line would be:
             Piecewise[{{1, dE <= 0}}, Exp[-\[Beta] * dE]]  *)
        If[RandomReal[] < Piecewise[{{1, dE <= 0}}, Exp[-\[Beta] * dE / 2]],
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
