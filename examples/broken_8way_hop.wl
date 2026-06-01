(* ================================================================
   broken_8way_hop.wl  —  BROKEN single-particle 8-hop (no rejection)
   ================================================================
   Each step:
     1. Choose a particle uniformly at random.
     2. Collect all 8-connected neighbours that are unoccupied.
     3. Choose one of those positions to hop to.
     4. No Metropolis step (energy = 0, pi uniform).

   Why this is wrong:
     With 2 particles on a 4x4 grid, if the two particles are
     8-adjacent the chosen particle has 7 valid moves; if not
     8-adjacent it has 8.  Both 7 and 8 use k = 3 BFS bits
     (IntegerLength[6,2] = IntegerLength[7,2] = 3), so the BFS
     assigns weight 1/8 to each leaf regardless of true pool size.
     True probability: 1/7 (adjacent) or 1/8 (non-adjacent).

   This triggers the same RandomChoice weight bug as
   broken_variable_pool.wl but in the k=3 bracket rather than k=2,
   confirming the bug affects every non-power-of-2 pool size.

   Note: must use nGrid >= 4 so that non-8-adjacent pairs exist.
   On the 3x3 torus every pair of sites is within 8-distance 1
   (max |dr|=1, max |dc|=1), making the pool uniformly 7 and the
   algorithm accidentally correct.

   Expected result (before fix): tau PASS, DB PASS  -- false positive!
   Expected result (after fix):  tau PASS, DB FAIL  -- correct
   ================================================================ *)

$nGrid         = 4;
numBeta        = 1.0;
$maxD2         = 0;          (* zero energy: pi uniform, T must be symmetric *)
$particleTypes = {1, 2};

$seedState = Module[
  {pos = Flatten[Table[{r, c}, {r, $nGrid}, {c, $nGrid}], 1],
   types = Sort[$particleTypes]},
  SortBy[MapThread[{#1, #2} &, {Take[pos, Length[types]], types}], First]];

$symmetryGroup = {"translation"};


(* ---- Energy ---- *)

energy[state_] := 0


(* ---- Algorithm ---- *)

(* 8-connectivity displacement vectors *)
$b8Disps = DeleteCases[
  Flatten[Table[{dx, dy}, {dx, -1, 1}, {dy, -1, 1}], 1],
  {0, 0}];

Algorithm[state_List] :=
  Module[{n = $nGrid, particle, rest, restPos, validDisps, disp, newPos, newState},
    particle = RandomChoice[state];
    rest     = DeleteCases[state, particle, 1];
    restPos  = #[[1]] & /@ rest;
    (* Filter to unoccupied 8-neighbours: pool size = 7 or 8 *)
    validDisps = Select[$b8Disps, Function[d,
      !AnyTrue[restPos, Function[p,
        Mod[p[[1]] - particle[[1, 1]] - d[[1]], n] === 0 &&
        Mod[p[[2]] - particle[[1, 2]] - d[[2]], n] === 0]]]];
    If[validDisps === {}, Return[state]];
    disp     = RandomChoice[validDisps];    (* BUG: pool 7 -> weight 1/8 != 1/7 *)
    newPos   = particle[[1]] + disp;
    newState = SortBy[Append[rest, {newPos, particle[[2]]}], First];
    newState]


(* ---- Symbolic parameters ---- *)

DynamicSymParams[states_List] :=
  <|"couplings" -> {}, "numericParams" -> {}|>
