#!/usr/bin/env wolframscript
(* ================================================================
   check.wls — SZ-DBC command-line checker
   ================================================================
   Usage:
     wolframscript -file check.wls <algorithm-file.wl> [options]

   Options:
     -mode Symbolic    (default) Full symbolic DB check (SZPure)
     -mode FullSimplify           Use FullSimplify instead of SZ
     -mode Numerical              Numerical MCMC check only
     -szRepeats N      SZ random evaluation points (default 30)
     -maxDepth N       Maximum BFS bit depth (default 22)
     -timeLimit T      Per-state time limit in seconds (default 120)
     -verbose          Print extra BFS progress info

   The algorithm file must define:
     Algorithm[state_]      — MCMC step, state = {{{r,c},type},...}
     energy[state_]         — total energy (no β factor)
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


(* ---- Helpers ---- *)
sep  = "================================================================";
sep2 = "----------------------------------------------------------------";

printHeader[title_] := (Print[sep]; Print["  ", title]; Print[sep])
printSection[title_] := (Print[sep2]; Print["  ", title]; Print[sep2])


(* ================================================================
   STEP 0 — Clear abstract parameters before symbolic BFS
   ================================================================ *)

symParamData = If[ValueQ[DynamicSymParams],
  DynamicSymParams[{seedState}],
  <|"couplings"->{}, "numericParams"->{}|>];
symCouplings = symParamData["couplings"];
symParamList = symCouplings;  (* atoms to keep abstract *)

(* Clear abstract param values so BFS remains symbolic *)
Scan[ClearAll, Map[ToString, symCouplings]];

(* Declare extra assumptions for DB check: all couplings real *)
extraAssumptions = Map[Element[#, Reals]&, symCouplings];

printHeader["SZ-DBC  —  " <> FileNameTake[algFile]];
Print["  nGrid      : ", nGrid];
Print["  particles  : ", $particleTypes];
Print["  symmetry   : ", symGroup];
Print["  mode       : ", mode];
Print["  seed state : ", seedState];
Print[""];


(* ================================================================
   STEP 1 — Translational invariance check (τ-BFS)
   ================================================================ *)

tauResult = <|"tauFree"->True|>;   (* default if no translation claimed *)

If[wantTrans && mode =!= "Numerical",
  printSection["Step 1: Translational invariance (τ-BFS)"];
  tauResult = $dbcCheckTranslational[seedState, Algorithm, nGrid, maxDepth, timeLimit];
  If[tauResult["tauFree"] === $Failed,
    Print["  ERROR during τ-BFS: ", tauResult["error"]];
    Print["  Continuing without translational symmetry speedup."];
    wantTrans = False,
  If[tauResult["tauFree"],
    Print["  PASS  — translation invariance algebraically verified"],
    Print["  FAIL  — τ appears in transition probabilities (translation-dependence detected)"];
    Print["  First violations:"];
    Scan[Print["    ", #]&, tauResult["violations"]];
    Print["  Continuing without translational symmetry."];
    wantTrans = False]]];


(* ================================================================
   STEP 2 — Enumerate N-particle state space
   ================================================================ *)

printSection["Step 2: State enumeration"];
t2 = AbsoluteTime[];
allStates = $dbcEnumerateNParticleStates[seedState, nGrid];
Print["  States found : ", Length[allStates], "  (", Round[AbsoluteTime[]-t2,0.01], "s)"];


(* ================================================================
   STEP 3 — Ergodicity check (state count)
   ================================================================ *)

ergResult = CheckErgodicity[allStates, nGrid];
Print["  Ergodicity   : ", If[ergResult["ergodic"],"PASS","FAIL"],
      "  (found ", ergResult["found"],
      ", theoretical ", ergResult["theoretical"], ")"];
If[!ergResult["ergodic"],
  Print["  WARNING: state count mismatch. Check $particleTypes and $nGrid."]];


(* ================================================================
   STEP 4 — NUMERICAL MODE (early exit)
   ================================================================ *)

If[doNumerical,
  printSection["Numerical MCMC check"];
  Print["  Running ", 100000, " steps at β=", numBetaVal, " ..."];
  t4 = AbsoluteTime[];
  counts   = RunNumericalMCMCAT[allStates, Algorithm, numBetaVal, nGrid];
  bweights = BoltzmannWeightsAT[allStates, energy, numBetaVal];
  nRun = Total[Values[counts]];
  maxDev = Max[Abs[N[Values[counts]]/nRun - N[Values[bweights]]]];
  Print["  Max deviation from Boltzmann: ", ScientificForm[maxDev, 3]];
  Print["  Time: ", Round[AbsoluteTime[]-t4, 0.01], "s"];
  Print[sep];
  Exit[0]];


(* ================================================================
   STEP 5 — G-orbit computation & BFS from orbit reps
   ================================================================ *)

printSection["Step 5: G-orbits & BFS"];
t5 = AbsoluteTime[];
allGElems = $dbcAllGElems[nGrid, If[wantTrans, symGroup, DeleteCases[symGroup,"translation"]]];
Print["  |G|   = ", Length[allGElems]];

{reps, repOfState, repToOrbitMap} = $dbcComputeOrbits[allStates, allGElems, nGrid];
Print["  Orbit reps : ", Length[reps]];

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
   STEP 6 — D4 generator testing
   ================================================================ *)

d4Result = <|"d4Pass"->True,"fails"->{}|>;

If[wantD4,
  printSection["Step 6: D4 generator testing"];
  t6 = AbsoluteTime[];
  d4Result = $dbcVerifyD4[repLeaves, Algorithm, nGrid, maxDepth, timeLimit];
  If[d4Result["d4Pass"],
    Print["  PASS  — D4 invariance verified (", Length[reps], " reps × 2 generators)"],
    Print["  FAIL  — D4 invariance violated:"];
    Scan[Print["    ", #]&, d4Result["fails"]]];
  Print["  Time: ", Round[AbsoluteTime[]-t6, 0.01], "s"]];


(* ================================================================
   STEP 7 — Detailed balance check
   ================================================================
   SZPure (default): evaluate leaf weights numerically at random
     {β, J...} points — Piecewise conditions resolve immediately,
     no symbolic matrix needed.  Fast O(reps × leaves × |G|).
   FullSimplify: build full symbolic matrix then FS each pair.
   ================================================================ *)

printSection["Step 7: Detailed balance (" <> If[useFS,"FullSimplify","SZPure direct"] <> ")"];
t7 = AbsoluteTime[];

If[useFS,
  (* FS path: needs symbolic T matrix *)
  fullMatrix = $dbcExpandOrbitsToMatrix[repLeaves, repToOrbitMap, nGrid];
  Print["  T-matrix entries : ", Length[fullMatrix]];
  symEnergyFn = Function[st, Block[{$nGrid=nGrid}, energy[st] /. r_Real :> Rationalize[r]]];
  violations = CheckDetailedBalanceSZ[
    fullMatrix, allStates, symEnergyFn, extraAssumptions,
    "FullSimplify"->True, "SZRepeats"->szRepeats],

  (* SZPure fast path: numerical leaf evaluation *)
  violations = $dbcSZCheckLeaves[
    repLeaves, repToOrbitMap,
    Function[st, Block[{$nGrid=nGrid}, energy[st]]],
    allStates, nGrid, symCouplings, szRepeats]];

Print["  Time: ", Round[AbsoluteTime[]-t7, 0.01], "s"];
If[violations === {},
  Print["  PASS  — detailed balance satisfied for all ", Length[allStates], " states"],
  Print["  FAIL  — ", Length[violations], " violation(s):"];
  Scan[Function[v,
    Print["    s=",v["pair"][[1]]," t=",v["pair"][[2]]];
    Print["      T(s->t)=",v["tij"],"  T(t->s)=",v["tji"]];
    Print["      E(s)=",v["ei"],"  E(t)=",v["ej"]]],
    violations]];


(* ================================================================
   SUMMARY
   ================================================================ *)

Print[sep];
Print["  SUMMARY"];
Print[sep];
Print["  Translational : ", If[wantTrans,
  If[tauResult["tauFree"],"PASS (algebraic)","FAIL"],
  "not checked"]];
Print["  D4 rotational : ", If[wantD4,
  If[d4Result["d4Pass"],"PASS","FAIL"],
  "not checked"]];
Print["  Ergodicity    : ", If[ergResult["ergodic"],"PASS","FAIL"],
      " (", ergResult["found"], "/", ergResult["theoretical"], " states)"];
Print["  Detailed bal. : ", If[violations==={},"PASS","FAIL"]];
Print[sep];

If[violations==={}&&(d4Result["d4Pass"])&&(tauResult["tauFree"]===True||!wantTrans),
  Exit[0], Exit[1]]
