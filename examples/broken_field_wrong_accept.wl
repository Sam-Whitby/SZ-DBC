(* ================================================================
   broken_field_wrong_accept.wl  —  Field algorithm with wrong acceptance
   ================================================================
   Based on quadratic_field.wl (nGrid=2, types {1,2}).

   The DECLARED energy includes a quadratic external field:
     energy[s] = Σ couplingJ(pair) + fieldH · Σ row_i²

   But the ALGORITHM'S Metropolis acceptance uses ONLY the pairwise
   energy change, completely ignoring the field contribution:
     ΔE_accept = $bfPairEnergy(new) - $bfPairEnergy(old)   ← wrong

   Why this breaks detailed balance:
     For a pair (s,t) where particle moves from row r to row r+dr
     with dr ≠ 0 and pairwise ΔE = 0 (equal-energy pairwise config):
       T(s→t) = (1/N)(1/4) · 1         (ΔE_pair = 0, always accept)
       T(t→s) = (1/N)(1/4) · 1         (ΔE_pair = 0, always accept)
     DB requires π(t)/π(s) = Exp[-β·ΔE_field]
       where ΔE_field = fieldH·(2r·dr + dr²) ≠ 0
     But T(s→t)/T(t→s) = 1 ≠ Exp[-β·ΔE_field].

   Note: τ-check PASSES because the acceptance uses only pairwise ΔE
   (τ-free). The translation orbit expansion for T is valid.
   The DB check uses the full declared energy (including field),
   so it correctly detects the mismatch.

   Expected result:  τ PASS,  DB FAIL
   ================================================================ *)


$nGrid         = 2;
numBeta        = 1.0;
$maxD2         = 2;

$particleTypes = {1, 2};

$seedState = Module[
  {pos = Flatten[Table[{r, c}, {r, $nGrid}, {c, $nGrid}], 1],
   types = Sort[$particleTypes]},
  SortBy[MapThread[{#1, #2} &, {Take[pos, Length[types]], types}], First]];

$symmetryGroup = {"translation"};


(* ---- Geometry ---- *)

$bfDisps = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};

$bfPairD2[p1_, p2_, n_] :=
  With[{dr = p1[[1]] - p2[[1]], dc = p1[[2]] - p2[[2]]},
    With[{dra = Mod[Abs[dr], n], dca = Mod[Abs[dc], n]},
      Min[dra, n - dra]^2 + Min[dca, n - dca]^2]]

$bfPairEnergy[state_, n_] :=
  Module[{total = 0},
    Do[
      With[{d2 = $bfPairD2[state[[i, 1]], state[[j, 1]], n]},
        If[0 < d2 <= $maxD2,
          total += couplingJ[state[[i, 2]], state[[j, 2]], d2]]],
      {i, Length[state]}, {j, i + 1, Length[state]}];
    total]

(* Declared energy: pairwise + quadratic field.
   This is what the checker uses for the DB condition. *)
energy[state_] :=
  $bfPairEnergy[state, $nGrid] +
  fieldH * Total[#[[1, 1]]^2 & /@ state]


(* ---- Algorithm ---- *)

Algorithm[state_List] :=
  Module[{n = $nGrid, particle, dir, newPos, rest, dE},
    particle = RandomChoice[state];
    dir      = RandomChoice[$bfDisps];
    newPos   = particle[[1]] + dir;
    rest     = DeleteCases[state, particle, 1];
    If[AnyTrue[rest, Function[p,
        Mod[p[[1, 1]] - newPos[[1]], n] === 0 &&
        Mod[p[[1, 2]] - newPos[[2]], n] === 0]],
      Return[state]];
    With[{newState = SortBy[Append[rest, {newPos, particle[[2]]}], First]},
      (* BUG: acceptance uses only pairwise ΔE, ignoring the field.
         The correct acceptance should use:
           dE = energy[newState] - energy[state]
         Instead we use only the pair contribution: *)
      With[{dE = $bfPairEnergy[newState, n] - $bfPairEnergy[state, n]},
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
    <|"couplings"      -> atoms,
      "extraSymParams" -> {fieldH},
      "numericParams"  -> {}|>]
