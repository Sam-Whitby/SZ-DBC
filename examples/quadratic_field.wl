(* ================================================================
   quadratic_field.wl  —  Single-particle Metropolis with a quadratic field
   ================================================================
   Energy = pairwise couplingJ + fieldH * Σ row_i²

   A LINEAR field (Σ row_i) would leave ΔE translation-invariant
   because the gradient is constant.  A QUADRATIC field (Σ row_i²)
   makes ΔE depend on the particle's absolute row:

     ΔE_field = fieldH * ((r+dr)² - r²) = fieldH * (2r·dr + dr²)

   With τ-augmented positions r → r+τr:
     ΔE_field = fieldH * (2(r+τr)·dr + dr²)

   τr appears in the Boltzmann factor → τ-check FAILS.

   After the τ-failure the checker recomputes orbits with the identity
   group only (all 12 states are their own representative) and runs
   the full DB check, which should PASS for a correctly implemented
   Metropolis accept/reject using the full energy.

   Uses nGrid=2 so the post-failure full-state BFS remains fast.

   Expected result:  τ FAIL,  DB PASS (after orbit recompute)
   ================================================================ *)


$nGrid         = 2;
numBeta        = 1.0;
$maxD2         = 2;

$particleTypes = {1, 2};

$seedState = Module[
  {pos = Flatten[Table[{r, c}, {r, $nGrid}, {c, $nGrid}], 1],
   types = Sort[$particleTypes]},
  SortBy[MapThread[{#1, #2} &, {Take[pos, Length[types]], types}], First]];

$symmetryGroup = {"translation"};   (* translation is declared, but will fail *)


(* ---- Geometry ---- *)

$qfDisps = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};

$qfPairD2[p1_, p2_, n_] :=
  With[{dr = p1[[1]] - p2[[1]], dc = p1[[2]] - p2[[2]]},
    With[{dra = Mod[Abs[dr], n], dca = Mod[Abs[dc], n]},
      Min[dra, n - dra]^2 + Min[dca, n - dca]^2]]

$qfPairEnergy[state_, n_] :=
  Module[{total = 0},
    Do[
      With[{d2 = $qfPairD2[state[[i, 1]], state[[j, 1]], n]},
        If[0 < d2 <= $maxD2,
          total += couplingJ[state[[i, 2]], state[[j, 2]], d2]]],
      {i, Length[state]}, {j, i + 1, Length[state]}];
    total]

(* Total energy: pairwise interactions + quadratic external field.
   fieldH is a symbolic coupling constant (declared in DynamicSymParams).
   The quadratic field Σ row_i² is NOT translation-invariant:
   ΔE_field contains the absolute row of the moved particle. *)
energy[state_] :=
  $qfPairEnergy[state, $nGrid] +
  fieldH * Total[#[[1, 1]]^2 & /@ state]


(* ---- Algorithm ---- *)

Algorithm[state_List] :=
  Module[{n = $nGrid, particle, dir, newPos, newPosNorm, rest},
    particle = RandomChoice[state];
    dir      = RandomChoice[$qfDisps];
    newPos   = particle[[1]] + dir;
    rest     = DeleteCases[state, particle, 1];

    If[AnyTrue[rest, Function[p,
        Mod[p[[1, 1]] - newPos[[1]], n] === 0 &&
        Mod[p[[1, 2]] - newPos[[2]], n] === 0]],
      Return[state]];

    (* Normalize to canonical grid coordinates before computing energy.
       The quadratic field uses absolute row², which differs for
       out-of-range positions — must use the PBC-wrapped value. *)
    newPosNorm = {Mod[newPos[[1]] - 1, n] + 1, Mod[newPos[[2]] - 1, n] + 1};
    With[{newState = SortBy[Append[rest, {newPosNorm, particle[[2]]}], First]},
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
    (* fieldH is declared in extraSymParams so it is assigned a random
       numeric value during the SZ check but is NOT subject to the
       coupling symmetry rule (it has no type indices). *)
    <|"couplings"     -> atoms,
      "extraSymParams" -> {fieldH},
      "numericParams"  -> {}|>]
