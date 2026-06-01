#!/usr/bin/env wolframscript
(* ================================================================
   check.wls — SZ-DBC command-line checker
   ================================================================
   Usage:
     wolframscript -file check.wls <algorithm-file.wl> [options]

   Options:
     -mode Symbolic    (default) EBE exact DB check
     -mode SZPure                Probabilistic SZ check (legacy)
     -mode FullSimplify          Use FullSimplify instead of EBE
     -mode Numerical             Numerical MCMC check only
     -szRepeats N      Random evaluation points for SZ fallback (default 30)
     -ebeMaxK N        Max branch conditions for EBE; falls back to SZ above (default 50)
     -maxDepth N       Maximum BFS bit depth (default 22)
     -timeLimit T      Per-state time limit in seconds (default 120)
     -verbose          Print extra BFS progress info

   The algorithm file must define:
     Algorithm[state_]      — MCMC step, state = {{{r,c},type},...}
     energy[state_]         — bare energy (no β factor)
     $nGrid                 — integer lattice side length
     $particleTypes         — type multiset, e.g. {1,2,3}
     $seedState             — canonical seed state
     $symmetryGroup         — list of "translation", "D4" (optional)
     numBeta                — inverse temperature (for Numerical mode)
     DynamicSymParams[states_] — returns <|"couplings"->{...},...|>
   ================================================================ *)

(* ---- Parse command line ---- *)
argv = Rest[$ScriptCommandLine];
If[Length[argv] < 1 || argv[[1]] === "-help",
  Print["Usage: wolframscript -file check.wls <algorithm-file.wl> [options]"];
  Exit[0]];
algFile = argv[[1]];
If[!FileExistsQ[algFile],
  Print["ERROR: file not found: ", algFile]; Exit[1]];

parseOpt[flag_,default_] :=
  With[{pos = Position[argv,flag]},
    If[pos==={}, default, argv[[pos[[1,1]]+1]]]];
parseBool[flag_] := MemberQ[argv, flag];

mode       = parseOpt["-mode", "Symbolic"];
szRepeats  = ToExpression[parseOpt["-szRepeats", "30"]];
maxDepth   = ToExpression[parseOpt["-maxDepth",  "22"]];
timeLimit  = ToExpression[parseOpt["-timeLimit", "120"]];
ebeMaxK    = ToExpression[parseOpt["-ebeMaxK",   "50"]];
verbose    = parseBool["-verbose"];


(* ---- Load core library ---- *)
coreFile = FileNameJoin[{DirectoryName[$InputFileName], "dbc_core.wl"}];
If[!FileExistsQ[coreFile],
  Print["ERROR: dbc_core.wl not found at ", coreFile]; Exit[1]];
Get[coreFile];


(* ---- Load algorithm file ---- *)
Get[algFile];


(* ---- Validate required globals ---- *)
checkDefined[sym_String] :=
  If[!ValueQ[ToExpression[sym]],
    Print["ERROR: ", sym, " not defined in ", algFile]; Exit[1]];
Scan[checkDefined, {"$nGrid","$particleTypes","$seedState","Algorithm","energy"}];

nGrid       = $nGrid;
seedState   = $seedState;
symGroup    = If[ValueQ[$symmetryGroup], $symmetryGroup, {}];
wantTrans   = MemberQ[symGroup, "translation"];
wantD4      = MemberQ[symGroup, "D4"];
numBetaVal  = If[ValueQ[numBeta], numBeta, 1.0];
useFS       = (mode === "FullSimplify");
doNumerical = (mode === "Numerical");
(* "Symbolic" and "SZPure" both use the symbolic BFS; Symbolic then
   applies EBE, SZPure applies the legacy probabilistic check.      *)


(* ---- Coupling atoms: clear any pre-existing numeric assignments,
        then install the canonicalisation symmetry rule so that
        couplingJ[b,a,d] with b>a → couplingJ[a,b,d] throughout. ---- *)
(* ValueQ[DynamicSymParams] would return False even when the function is defined,
   because ValueQ only tests OwnValues, not DownValues.  A function f[x_]:=...
   has DownValues but no OwnValue.  Use DownValues to detect definition. *)
symParamData = If[Length[DownValues[DynamicSymParams]] > 0,
  DynamicSymParams[{seedState}],
  <|"couplings"->{}, "numericParams"->{}|>];
symCouplings  = symParamData["couplings"];
(* extraSymParams: bare symbols for field-like parameters (e.g. fieldH).
   Unlike coupling atoms, these are not canonicalised by the symmetry rule. *)
extraSymParams = Lookup[symParamData, "extraSymParams", {}];
symParamList  = Join[symCouplings, extraSymParams];

(* Clear the coupling head symbol(s) to remove any residual assignments *)
Scan[ClearAll, DeleteDuplicates[Head /@ symCouplings]];
(* Clear extra parameter symbols directly *)
Scan[ClearAll, extraSymParams];

(* Add symmetry rule: couplingJ[b,a,d]/;b>a := couplingJ[a,b,d] *)
$dbcAddCouplingSymmetry[symCouplings];

(* Validate: all coupling atom arguments must be machine integers.
   The symmetry rule and SZ substitution only work correctly for integer
   type labels {1,2,3,...}.  Reject non-integer labels immediately. *)
If[Length[symCouplings] > 0 &&
   !AllTrue[symCouplings, MatchQ[#, _[_Integer, _Integer, _Integer]]&],
  Print["ERROR: coupling atom arguments must all be integers."];
  Print["  Particle types must be positive integers {1, 2, 3, ...}."];
  Print["  Got: ", symCouplings];
  Exit[1]];

(* Declare extra assumptions for DB check: all params real *)
extraAssumptions = Map[Element[#, Reals]&, symParamList];


(* ---- Helpers ---- *)
sep  = "================================================================";
sep2 = "----------------------------------------------------------------";

printHeader[title_] := (Print[sep]; Print["  ", title]; Print[sep])
printSection[title_] := (Print[sep2]; Print["  ", title]; Print[sep2])


printHeader["SZ-DBC  —  " <> FileNameTake[algFile]];
Print["  nGrid      : ", nGrid];
Print["  particles  : ", $particleTypes];
Print["  symmetry   : ", symGroup];
Print["  mode       : ", mode];
Print["  seed state : ", seedState];
Print[""];


(* ================================================================
   STEP 1 — Enumerate N-particle state space + combinatorial count
   ================================================================ *)

printSection["Step 1: State enumeration"];
t1 = AbsoluteTime[];
allStates = $dbcEnumerateNParticleStates[seedState, nGrid];
Print["  States found : ", Length[allStates], "  (", Round[AbsoluteTime[]-t1,0.01], "s)"];

(* Combinatorial sanity check: validates $nGrid and $particleTypes *)
ergCount = CheckErgodicity[allStates, nGrid];
Print["  State count  : ", If[ergCount["ergodic"],"OK","MISMATCH"],
      "  (found ", ergCount["found"],
      ", theoretical ", ergCount["theoretical"], ")"];
If[!ergCount["ergodic"],
  Print["  WARNING: state count mismatch. Check $particleTypes and $nGrid."]];


(* ================================================================
   STEP 2 — G-orbit computation
   ================================================================ *)

printSection["Step 2: G-orbit computation"];
t2 = AbsoluteTime[];
(* Always start with translation-only orbits.  D4/C4 is verified
   numerically after the BFS (Step 5b) and orbits are recomputed
   to the verified subgroup before the DB check. *)
allGElems = $dbcAllGElems[nGrid, Intersection[symGroup, {"translation"}]];
Print["  |G| (initial) = ", Length[allGElems], " (translation only)"];
{reps, repOfState, repToOrbitMap} = $dbcComputeOrbits[allStates, allGElems, nGrid];
Print["  Orbit reps : ", Length[reps], "  (", Round[AbsoluteTime[]-t2,0.01], "s)"];


(* ================================================================
   STEP 3 — NUMERICAL MODE (early exit)
   ================================================================ *)

If[doNumerical,
  printSection["Numerical MCMC check"];
  Print["  Running ", 100000, " steps at β=", numBetaVal, " ..."];
  t3 = AbsoluteTime[];
  counts   = RunNumericalMCMCAT[allStates, Algorithm, numBetaVal, nGrid];
  bweights = BoltzmannWeightsAT[allStates, energy, numBetaVal];
  nRun = Total[Values[counts]];
  maxDev = Max[Abs[N[Values[counts]]/nRun - N[Values[bweights]]]];
  Print["  Max deviation from Boltzmann: ", ScientificForm[maxDev, 3]];
  Print["  Time: ", Round[AbsoluteTime[]-t3, 0.01], "s"];
  Print[sep];
  Exit[0]];


(* ================================================================
   STEP 4 — Translational invariance check (τ-BFS from orbit reps)
   ================================================================
   Runs $dbcBuildStateLeaves from each orbit rep with τ-augmented
   positions.  This covers all distinct states (via G-symmetry)
   at ~1/63 of the cost of the full state-space BFS.
   If translation invariance fails, orbit reps are recomputed
   without the translation subgroup.
   ================================================================ *)

tauResult = <|"tauFree"->True|>;

If[wantTrans && mode =!= "Numerical",
  printSection["Step 4: Translational invariance (τ-BFS from orbit reps)"];
  t4 = AbsoluteTime[];
  tauResult = $dbcCheckTranslational[reps, Algorithm, nGrid, maxDepth, timeLimit];
  Print["  States checked : ", Length[reps], " orbit reps  (",
        Round[AbsoluteTime[]-t4, 0.01], "s)"];
  If[tauResult["tauFree"] === $Failed,
    Print["  ERROR during τ-BFS: ", tauResult["error"]];
    Print["  Continuing without translational symmetry speedup."];
    wantTrans = False;
    allGElems = $dbcAllGElems[nGrid, DeleteCases[symGroup,"translation"]];
    {reps, repOfState, repToOrbitMap} = $dbcComputeOrbits[allStates, allGElems, nGrid];
    Print["  Recomputed orbit reps without translation: ", Length[reps]],
  If[tauResult["tauFree"],
    Print["  PASS  — translation invariance algebraically verified"],
    Print["  FAIL  — τ appears in transition probabilities"];
    Scan[Print["    sample violation: ", Short[#,5]]&, tauResult["violations"]];
    Print["  Recomputing orbit reps without translation symmetry."];
    wantTrans = False;
    allGElems = $dbcAllGElems[nGrid, DeleteCases[symGroup,"translation"]];
    {reps, repOfState, repToOrbitMap} = $dbcComputeOrbits[allStates, allGElems, nGrid];
    Print["  Orbit reps after recompute: ", Length[reps]]]]];


(* ================================================================
   STEP 5 — BFS from orbit representatives
   ================================================================ *)

printSection["Step 5: BFS from orbit reps"];
t5 = AbsoluteTime[];
normFn = $dbcNormState[#, nGrid]&;
repLeaves = <||>;
err = None;
Do[
  s      = reps[[ri]];
  leaves = $dbcBuildStateLeaves[s, Algorithm, maxDepth, timeLimit, nGrid, normFn];
  If[MatchQ[leaves, $dbc$cantHandle[_]],
    err = leaves[[1]]; Break[]];
  repLeaves[s] = leaves;
  If[verbose, Print["  rep ", ri, "/", Length[reps], " done"]],
  {ri, Length[reps]}];

If[err =!= None,
  Print["  ERROR in BFS: ", err]; Exit[1]];

Print["  BFS done: ", Length[reps], " reps  (",
      Round[AbsoluteTime[]-t5, 0.01], "s)"];


(* ================================================================
   STEP 5b — Rotational invariance check and orbit reduction
   ================================================================
   If "D4" is declared in $symmetryGroup AND translation invariance
   passed, verify D4 numerically using the T matrix built from the
   translation-only BFS.  If any rotation/reflection passes, recompute
   orbits with the verified subgroup and rerun BFS from the smaller
   representative set before the DB check.

   The check is SZ-style: builds T at nReps random coupling points
   and verifies T(s→t) = T(R·s→R·t) for all communicating pairs.
   This is reliable in practice; see documentation for discussion of
   the algebraic vs. numerical trade-off.
   ================================================================ *)

(* Helper: run BFS from a given rep list, update repLeaves in place *)
$runBFS[repsIn_List] :=
  Module[{leavesNew = <||>, errNew = None},
    Do[
      With[{sRep = repsIn[[ri]],
            lv = $dbcBuildStateLeaves[repsIn[[ri]], Algorithm, maxDepth,
                                       timeLimit, nGrid, normFn]},
        If[MatchQ[lv, $dbc$cantHandle[_]], errNew = lv[[1]]; Break[]];
        leavesNew[sRep] = lv],
      {ri, Length[repsIn]}];
    {leavesNew, errNew}];

d4Result = <|"d4Pass" -> False, "label" -> "not checked"|>;
rotCheck  = <|"pass" -> False, "ebe" -> False, "chambers" -> {}|>;

If[wantD4 && tauResult["tauFree"] === True,
  printSection["Step 5b: Rotational invariance check (D4)"];
  t5b = AbsoluteTime[];
  (* nReps used only when EBE falls back to numerical check *)
  d4CheckReps = Max[5, Min[10, Length[symParamList] + 2]];
  (* Check via D4 generators {rot90CW=1, LR-reflect=4}.
     By group theory: if T(s→t)=T(R·s→R·t) for both generators,
     it holds for all 8 D4 elements (generators span the group).
     C4 is verified by generator 1 alone; reflections require 4. *)
  rotCheck = $dbcEBECheckD4[
    repLeaves, repToOrbitMap, allStates, nGrid, symParamList,
    {1, 4}, ebeMaxK, d4CheckReps];
  Print["  Time: ", Round[AbsoluteTime[] - t5b, 0.01], "s"];

  Which[
    rotCheck["pass"],
      (* Full D4 verified *)
      (* rotCheck["violations"] is empty; determine if EBE or numerical was used *)
      If[rotCheck["ebe"],
        Print["  PASS — full D4 symmetry verified (algebraically exact via EBE)"],
        Print["  PASS — full D4 symmetry verified (", d4CheckReps, " random coupling points)"]];
      d4Result = <|"d4Pass" -> True, "label" -> If[rotCheck["ebe"], "D4 (exact)", "D4 (numerical)"]|>;
      allGElems = $dbcAllGElems[nGrid, symGroup];
      {reps, repOfState, repToOrbitMap} = $dbcComputeOrbits[allStates, allGElems, nGrid];
      Print["  |G| = ", Length[allGElems], "  →  ", Length[reps],
            " orbit reps after D4 reduction"];
      printSection["Step 5c: BFS from reduced orbit reps (D4)"];
      tReBFS = AbsoluteTime[];
      {repLeaves, err} = $runBFS[reps];
      If[err =!= None, Print["  ERROR in BFS: ", err]; Exit[1]];
      Print["  BFS done: ", Length[reps], " reps  (",
            Round[AbsoluteTime[] - tReBFS, 0.01], "s)"],

    (* C4 rotations only? Check elements 1-3 separately *)
    Select[rotCheck["violations"], #["rot"] <= 3 &] === {},
      If[rotCheck["ebe"],
        Print["  PASS — C4 (rotations) verified (algebraically exact); reflections failed"],
        Print["  PASS — C4 (rotations) verified; reflections failed"]];
      Print["  Using C4+translation orbits"];
      d4Result = <|"d4Pass" -> True, "label" -> If[rotCheck["ebe"], "C4 (exact)", "C4 (numerical)"]|>;
      allGElems = $dbcAllGElems[nGrid, Union[DeleteCases[symGroup, "D4"], {"C4"}]];
      {reps, repOfState, repToOrbitMap} = $dbcComputeOrbits[allStates, allGElems, nGrid];
      Print["  |G| = ", Length[allGElems], "  →  ", Length[reps],
            " orbit reps after C4 reduction"];
      printSection["Step 5c: BFS from reduced orbit reps (C4)"];
      tReBFS = AbsoluteTime[];
      {repLeaves, err} = $runBFS[reps];
      If[err =!= None, Print["  ERROR in BFS: ", err]; Exit[1]];
      Print["  BFS done: ", Length[reps], " reps  (",
            Round[AbsoluteTime[] - tReBFS, 0.01], "s)"],

    True,
      (* Neither C4 nor full D4 *)
      Print["  FAIL — D4 symmetry not verified; keeping translation-only orbits"];
      d4Result = <|"d4Pass" -> False, "label" -> "FAIL"|>;
      Scan[Print["    violation: rotElem=", #["rot"],
                 "  s=", Short[#["s"], 3],
                 "  T(s,t)=", #["Tst"],
                 "  T(Rs,Rt)=", #["TRsRt"]] &,
           Take[rotCheck["violations"], UpTo[3]]]]];



(* ================================================================
   STEP 7 — Ergodicity (reachability from orbit-BFS transition graph)
   ================================================================ *)

printSection["Step 7: Ergodicity (reachability check)"];
t7 = AbsoluteTime[];
ergoResult = $dbcCheckErgodicityFromLeaves[repLeaves, repToOrbitMap, allStates, nGrid, seedState];
Print["  Reachable : ", ergoResult["reached"], "/", ergoResult["total"],
      "  (", Round[AbsoluteTime[]-t7, 0.01], "s)"];
If[ergoResult["ergodic"],
  Print["  PASS  — all states reachable from seed"],
  Print["  FAIL  — only ", ergoResult["reached"], " of ", ergoResult["total"],
        " states are reachable."];
  Print["  The algorithm may not be ergodic for this system size."]];


(* ================================================================
   STEP 8 — Detailed balance check
   ================================================================
   SZPure (default): evaluate leaf weights numerically at random
     {β, J...} points — Piecewise conditions resolve immediately,
     no symbolic matrix needed.  Fast O(reps × leaves × |G|).
   FullSimplify: build full symbolic matrix then FS each pair.
   ================================================================ *)

printSection["Step 8: Detailed balance (" <>
  Which[useFS, "FullSimplify",
        mode === "SZPure", "SZPure (probabilistic)",
        True, "EBE (exact)"] <> ")"];
t8 = AbsoluteTime[];

If[useFS,
  fullMatrix = $dbcExpandOrbitsToMatrix[repLeaves, repToOrbitMap, nGrid];
  Print["  T-matrix entries : ", Length[fullMatrix]];
  symEnergyFn = Function[st, Block[{$nGrid=nGrid}, energy[st] /. r_Real :> Rationalize[r]]];
  violations = CheckDetailedBalanceSZ[
    fullMatrix, allStates, symEnergyFn, extraAssumptions,
    "FullSimplify"->True, "SZRepeats"->szRepeats],

If[mode === "SZPure",
  violations = $dbcSZCheckLeaves[
    repLeaves, repToOrbitMap,
    Function[st, Block[{$nGrid=nGrid}, energy[st]]],
    allStates, nGrid, symParamList, szRepeats],

  (* Default: EBE exact check — enumerates all coupling-constant
     branch regions and checks DB exactly with $dbcIsExpZero.
     Falls back to probabilistic SZ when k > ebeMaxK.
     If D4 EBE ran (d4Result["ebe"]==True), reuse its feasible
     regions to skip the expensive Phase 2 chamber enumeration.   *)
  violations = $dbcEBECheckLeaves[
    repLeaves, repToOrbitMap,
    Function[st, Block[{$nGrid=nGrid}, energy[st]]],
    allStates, nGrid, symParamList,
    ebeMaxK, szRepeats, 1*^-7,
    If[TrueQ[d4Result["d4Pass"]] && TrueQ[rotCheck["ebe"]],
       rotCheck["chambers"], {}]]]];

Print["  Time: ", Round[AbsoluteTime[]-t8, 0.01], "s"];
If[violations === {},
  Print["  PASS  — detailed balance satisfied for all ", Length[allStates], " states"],
  Print["  FAIL  — ", Length[violations], " violating pair(s) found:"];
  Scan[Function[v,
    Print["    s=",v["pair"][[1]]," t=",v["pair"][[2]]];
    Print["      T(s->t)=",v["tij"],"  T(t->s)=",v["tji"]];
    Print["      E(s)=",v["ei"],"  E(t)=",v["ej"]]],
    Take[violations, UpTo[5]]];
  If[Length[violations] >= 10,
    Print["    (stopped after 10 unique violating pairs — likely many more)"]]];


(* ================================================================
   SUMMARY
   ================================================================ *)

Print[sep];
Print["  SUMMARY"];
Print[sep];
Print["  Translational : ", If[wantTrans,
  If[tauResult["tauFree"],"PASS (algebraic)","FAIL"],
  If[MemberQ[symGroup,"translation"],"FAIL (recomputed without)","not checked"]]];
Print["  D4 rotational : ", Which[
  !wantD4, "not declared",
  d4Result["label"] === "not checked", "declared but not checked (τ failed)",
  True, d4Result["label"]]];
Print["  State count   : ", If[ergCount["ergodic"],"OK","MISMATCH"],
      " (", ergCount["found"], "/", ergCount["theoretical"], ")"];
Print["  Ergodicity    : ", If[ergoResult["ergodic"],"PASS","FAIL"],
      " (", ergoResult["reached"], "/", ergoResult["total"], " reachable)"];
Print["  Detailed bal. : ", If[violations==={},"PASS","FAIL"]];
Print[sep];

allPassed = (violations === {} &&
  (tauResult["tauFree"] === True || !MemberQ[symGroup, "translation"]) &&
  ergoResult["ergodic"] &&
  ergCount["ergodic"]);

If[allPassed, Exit[0], Exit[1]]
