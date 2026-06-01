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
     -ebeMaxK N        Max branch conditions for EBE; falls back to SZ above (default 12)
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
ebeMaxK    = ToExpression[parseOpt["-ebeMaxK",   "12"]];
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
allGElems = $dbcAllGElems[nGrid, symGroup];
Print["  |G| = ", Length[allGElems]];
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
   STEP 6 — D4 generator testing  [DISABLED — see note below]
   ================================================================
   D4 verification is currently disabled pending a robust implementation
   that handles orbit recomputation on failure (analogous to the τ-failure
   fallback in Step 4).

   When "D4" is declared in $symmetryGroup:
     - D4 IS used for orbit computation (Step 2) for performance: 8 reps
       instead of ~56 translation-only reps.
     - D4 equivariance is NOT independently verified here.
     - The DB check (Step 8) trusts the D4 orbit expansion implicitly.

   Consequence: if your algorithm declares D4 but actually breaks it,
   the DB check result may be unreliable (transition matrix is wrong).
   If you are unsure whether your algorithm is D4-symmetric, declare
   $symmetryGroup = {"translation"} only — this uses translation-only
   orbits (~56 reps for nGrid=3, slower) and verifies only translation.
   ================================================================ *)

d4Result = <|"d4Pass" -> True, "fails" -> {}|>;

If[wantD4,
  printSection["Step 6: D4 (assumed, not verified)"];
  Print["  D4 declared in $symmetryGroup — used for orbit reduction (", Length[reps], " reps)."];
  Print["  D4 equivariance is assumed, not independently verified."]];



(* ================================================================
   STEP 7 — Ergodicity (reachability from orbit-BFS transition graph)
   ================================================================ *)

printSection["Step 7: Ergodicity (reachability check)"];
t7 = AbsoluteTime[];
ergoResult = $dbcCheckErgodicityFromLeaves[repLeaves, repToOrbitMap, allStates, nGrid];
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
     Falls back to probabilistic SZ when k > ebeMaxK.             *)
  violations = $dbcEBECheckLeaves[
    repLeaves, repToOrbitMap,
    Function[st, Block[{$nGrid=nGrid}, energy[st]]],
    allStates, nGrid, symParamList,
    ebeMaxK, szRepeats]]];

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
Print["  D4 rotational : ", If[wantD4, "assumed (not verified)", "not declared"]];
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
