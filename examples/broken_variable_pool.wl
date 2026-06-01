(* ================================================================
   broken_variable_pool.wl  —  BROKEN single-particle 4-hop
   ================================================================
   Each step:
     1. Choose a particle uniformly at random.
     2. Collect all 4-connected neighbors that are unoccupied.
     3. Choose one of those positions to hop to.
     4. No Metropolis step (energy = 0, π uniform).

   Why this is wrong:
     With 2 particles on a 3×3 grid, if the two particles are
     adjacent (share a 4-neighbor) the chosen particle has 3 valid
     moves; if not adjacent it has 4 valid moves.  Both 3 and 4 use
     k = 2 BFS bits (IntegerLength[2,2]=IntegerLength[3,2]=2), so the
     BFS assigns weight 1/4 to each leaf regardless of true pool size.
     True probability: 1/3 (adjacent) or 1/4 (non-adjacent).

   Consequence:
     For the transition adjacent-state → non-adjacent-state and its
     reverse, the BFS computes (1/2)·(1/4) for both directions and
     reports PASS, but the true probabilities are (1/2)·(1/3) ≠
     (1/2)·(1/4), a genuine DB violation.

   Expected result (before fix): τ PASS, DB PASS  ← false positive!
   Expected result (after fix):  τ PASS, DB FAIL  ← correct
   ================================================================ *)

$nGrid         = 3;
numBeta        = 1.0;
$maxD2         = 0;          (* zero energy: π uniform, T must be symmetric *)
$particleTypes = {1, 2};

$seedState = Module[
  {pos = Flatten[Table[{r, c}, {r, $nGrid}, {c, $nGrid}], 1],
   types = Sort[$particleTypes]},
  SortBy[MapThread[{#1, #2} &, {Take[pos, Length[types]], types}], First]];

$symmetryGroup = {"translation"};  (* τ cancels by difference arithmetic *)


(* ---- Energy ---- *)

energy[state_] := 0


(* ---- Algorithm ---- *)

(* 4-connectivity displacement vectors *)
$bvpDisps = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};

Algorithm[state_List] :=
  Module[{n = $nGrid, particle, rest, restPos, validDisps, disp, newPos, newState},
    particle = RandomChoice[state];
    rest     = DeleteCases[state, particle, 1];
    restPos  = #[[1]] & /@ rest;
    (* Filter to unoccupied 4-neighbours: pool size = 3 or 4 *)
    validDisps = Select[$bvpDisps, Function[d,
      !AnyTrue[restPos, Function[p,
        Mod[p[[1]] - particle[[1, 1]] - d[[1]], n] === 0 &&
        Mod[p[[2]] - particle[[1, 2]] - d[[2]], n] === 0]]]];
    If[validDisps === {}, Return[state]];
    disp     = RandomChoice[validDisps];    (* BUG: pool 3 → weight 1/4 ≠ 1/3 *)
    newPos   = particle[[1]] + disp;        (* unnormalised; checker applies PBC *)
    newState = SortBy[Append[rest, {newPos, particle[[2]]}], First];
    newState]


(* ---- Symbolic parameters ---- *)

DynamicSymParams[states_List] :=
  <|"couplings" -> {}, "numericParams" -> {}|>
