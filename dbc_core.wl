(* ================================================================
   SZ-DBC — Detailed Balance Checker (core library)
   ================================================================
   State format: sorted list of {{row,col}, type} pairs.
   Positions are integers in {1,..,nGrid}.

   Key features:
   - TRANSLATIONAL INVARIANCE verified algebraically via τ-BFS:
     checker augments all positions with a free symbol {τr,τc};
     τ cancels in pairwise differences so transition probabilities
     are τ-free iff the algorithm is translation-invariant.
   - D4 ROTATIONAL INVARIANCE verified by comparing BFS outputs
     for each orbit rep with its rotation and reflection.
   - SCHWARTZ-ZIPPEL probabilistic DB check (default, no FS).
   - FULL-SIMPLIFY available as an option.

   Load with:  Get["path/to/dbc_core.wl"]
   ================================================================ *)

$dbcDir = DirectoryName[$InputFileName];
$dbcCurrentNGrid = 1;


(* ================================================================
   SECTION 0 — STATE FORMAT UTILITIES  (position-type pairs)
   ================================================================ *)

(* Normalise a state: apply PBC to all positions, sort by position.
   Works for plain integer positions or τ-augmented positions. *)
$dbcNormState[state_List, nGrid_Integer] :=
  SortBy[
    Map[{{Mod[#[[1,1]]-1, nGrid]+1, Mod[#[[1,2]]-1, nGrid]+1}, #[[2]]} &, state],
    First]

(* Free symbols for translational τ-BFS *)
τr; τc;

$dbcAddTau[state_List] :=
  Map[{{#[[1,1]] + τr, #[[1,2]] + τc}, #[[2]]} &, state]

(* Normalise a τ-augmented state: PBC on integer part, sort by canonical {row,col}.
   Subtracting τr/τc extracts the integer part since τ cancels algebraically. *)
$dbcNormTauState[state_List, nGrid_Integer] :=
  SortBy[
    Map[
      {{Mod[#[[1,1]]-τr-1, nGrid]+1 + τr,
        Mod[#[[1,2]]-τc-1, nGrid]+1 + τc},
       #[[2]]} &, state],
    {#[[1,1]]-τr, #[[1,2]]-τc} &]


(* ================================================================
   SECTION 1 — GROUP ACTIONS ON STATES
   ================================================================ *)

$dbcApplyRot90[state_List, nGrid_Integer] :=
  $dbcNormState[Map[{{#[[1,2]], nGrid+1-#[[1,1]]}, #[[2]]} &, state], nGrid]

$dbcApplyReflect[state_List, nGrid_Integer] :=
  $dbcNormState[Map[{{#[[1,2]], #[[1,1]]}, #[[2]]} &, state], nGrid]

$dbcApplyTrans[state_List, dr_Integer, dc_Integer, nGrid_Integer] :=
  $dbcNormState[
    Map[{{Mod[#[[1,1]]+dr-1, nGrid]+1, Mod[#[[1,2]]+dc-1, nGrid]+1}, #[[2]]} &, state],
    nGrid]

(* D4 element: rotIdx 0-7 (id, rot90CW, rot180, rot270, reflLR, reflUD, reflDiag, reflAnti) *)
$dbcD4Pos[0, r_, c_, n_] := {r, c}
$dbcD4Pos[1, r_, c_, n_] := {c, n+1-r}
$dbcD4Pos[2, r_, c_, n_] := {n+1-r, n+1-c}
$dbcD4Pos[3, r_, c_, n_] := {n+1-c, r}
$dbcD4Pos[4, r_, c_, n_] := {r, n+1-c}
$dbcD4Pos[5, r_, c_, n_] := {n+1-r, c}
$dbcD4Pos[6, r_, c_, n_] := {c, r}
$dbcD4Pos[7, r_, c_, n_] := {n+1-c, n+1-r}

$dbcApplyGElem[<|"rot"->ri_, "dr"->dr_, "dc"->dc_|>, state_List, nGrid_Integer] :=
  $dbcNormState[Map[
    With[{rc2 = $dbcD4Pos[ri, #[[1,1]], #[[1,2]], nGrid]},
      {{Mod[rc2[[1]]+dr-1, nGrid]+1, Mod[rc2[[2]]+dc-1, nGrid]+1}, #[[2]]}] &,
    state], nGrid]

$dbcAllGElems[nGrid_Integer, symGroup_List] :=
  Module[{rotIdxs, drDcs},
    rotIdxs = If[MemberQ[symGroup, "D4"], Range[0,7], {0}];
    drDcs   = If[MemberQ[symGroup, "translation"],
      Flatten[Table[{dr,dc},{dr,0,nGrid-1},{dc,0,nGrid-1}],1],
      {{0,0}}];
    Flatten[Table[<|"rot"->ri,"dr"->drdc[[1]],"dc"->drdc[[2]]|>,
      {ri,rotIdxs},{drdc,drDcs}],1]]


(* ================================================================
   SECTION 2 — ORBIT COMPUTATION
   ================================================================ *)

$dbcEnumerateNParticleStates[seedState_List, nGrid_Integer] :=
  Module[{types, allPos},
    types  = Sort[#[[2]] & /@ seedState];
    allPos = Flatten[Table[{r,c},{r,nGrid},{c,nGrid}],1];
    DeleteDuplicates @ Map[
      Function[sites, SortBy[MapThread[{#1,#2}&,{sites,types}], First]],
      Permutations[allPos, {Length[types]}]]]

(* Group allStates into G-orbits.  Returns {reps, repOfState, repToOrbitMap}. *)
$dbcComputeOrbits[allStates_List, allGElems_List, nGrid_Integer] :=
  Module[{stateSet, unvisited, repOfState, reps, repToOrbitMap,
          s, imgPairs, minImg, orbitMap},
    stateSet      = Association[# -> True & /@ allStates];
    unvisited     = Association[# -> True & /@ allStates];
    repOfState    = <||>; reps = {}; repToOrbitMap = <||>;
    While[Length[unvisited] > 0,
      s        = First[Keys[unvisited]];
      imgPairs = {$dbcApplyGElem[#, s, nGrid], #} & /@ allGElems;
      minImg   = First @ SortBy[imgPairs, First][[All,1]];
      Scan[Function[ip,
          repOfState[ip[[1]]] = minImg;
          unvisited = Delete[unvisited, Key[ip[[1]]]]],
        imgPairs];
      AppendTo[reps, minImg];
      orbitMap = <||>;
      Scan[Function[ip,
          With[{img=ip[[1]], g=ip[[2]]},
            If[KeyExistsQ[stateSet,img] && !KeyExistsQ[orbitMap,img],
              orbitMap[img] = g]]],
        imgPairs];
      repToOrbitMap[minImg] = orbitMap];
    {reps, repOfState, repToOrbitMap}]

(* Expand per-rep BFS leaves to the full T matrix via G-action. *)
$dbcExpandOrbitsToMatrix[repLeaves_Association, repToOrbitMap_Association,
                          nGrid_Integer] :=
  Module[{fullMatrix = <||>, repMatrix, orbitMap, gElem, sp, spNew},
    Do[
      repMatrix = <||>;
      Do[With[{ns=leaf[[2]], w=leaf[[3]]},
        repMatrix[{repKey,ns}] = Lookup[repMatrix,Key[{repKey,ns}],0] + w],
        {leaf, repLeaves[repKey]}];
      orbitMap = repToOrbitMap[repKey];
      Do[
        gElem = orbitMap[sKey];
        KeyValueMap[Function[{repSp, tVal},
          sp    = repSp[[2]];
          spNew = $dbcApplyGElem[gElem, sp, nGrid];
          fullMatrix[{sKey, spNew}] = tVal],
          repMatrix],
        {sKey, Keys[orbitMap]}],
      {repKey, Keys[repLeaves]}];
    fullMatrix]


(* ================================================================
   SECTION 3 — BFS ENGINE
   ================================================================ *)

$OutOfBits      = Symbol["$OutOfBits"];
$dbc$tag        = Symbol["$dbc$tag"];
$dbc$outOfRange = Symbol["$dbc$outOfRange"];
$dbc$cantHandle = Symbol["$dbc$cantHandle"];

(* UpValues: $dbc$irand[j, at] comparison → acceptTestI call *)
$dbc$irand /: Less[$dbc$irand[j_, at_], p_]         :=  at[j, p]
$dbc$irand /: LessEqual[$dbc$irand[j_, at_], p_]    :=  at[j, p]
$dbc$irand /: Greater[$dbc$irand[j_, at_], p_]      := !at[j, p]
$dbc$irand /: GreaterEqual[$dbc$irand[j_, at_], p_] := !at[j, p]
$dbc$irand /: Less[p_, $dbc$irand[j_, at_]]         := !at[j, p]
$dbc$irand /: LessEqual[p_, $dbc$irand[j_, at_]]    := !at[j, p]
$dbc$irand /: Greater[p_, $dbc$irand[j_, at_]]      :=  at[j, p]
$dbc$irand /: GreaterEqual[p_, $dbc$irand[j_, at_]] :=  at[j, p]

(* Continuous token arithmetic UpValues *)
$dbc$contToken /: Plus[before___, $dbc$contToken["Uniform", lo_, hi_, sb_], after___] :=
  With[{offset = Plus @@ {before, after}},
    If[!FreeQ[offset, $dbc$contToken],
      Throw[$dbc$cantHandle["Adding two random tokens: unsupported"], $dbc$tag],
      $dbc$contToken["Uniform", lo+offset, hi+offset, sb]]]

$dbc$contToken /: Mod[$dbc$contToken["Uniform", lo_, hi_, sb_], L_] :=
  With[{span = hi-lo},
    If[TrueQ[span === L],
      $dbc$contToken["Uniform", 0, L, sb],
      Throw[$dbc$cantHandle["Mod[token," <> ToString[L] <> "]: span " <>
                            ToString[span] <> " ≠ modulus"], $dbc$tag]]]

$dbc$contToken /: Floor[$dbc$contToken["Uniform", 0, L_Integer, sb_]] :=
  sb[Table[1/L,{L}], Range[0,L-1]]
$dbc$contToken /: Round[$dbc$contToken["Uniform", 0, L_Integer, sb_]] :=
  sb[Table[1/L,{L}], Range[0,L-1]]
$dbc$contToken /: Ceiling[$dbc$contToken["Uniform", 0, L_Integer, sb_]] :=
  sb[Table[1/L,{L}], Range[1,L]]

(* ----------------------------------------------------------------
   RunWithBitsAT
   Run alg[state] against a fixed bit sequence.
   Uses interval-tracking for RandomReal[]: each call creates a
   fresh token $dbc$irand[j, acceptTestI].  Comparisons narrow the
   interval and accumulate the path weight exactly.
   RandomChoice[w->e] decomposes into sequential Bernoulli trials.
   Returns {nextState, pathWeight} or $OutOfBits / $dbc$outOfRange /
   $dbc$cantHandle[msg].
   ---------------------------------------------------------------- *)

$dbc$readBitsAsInt[k_Integer, readBit_] :=
  Fold[#1*2 + readBit[] &, 0, Range[k]]

RunWithBitsAT[alg_, state_, bits_List] :=
  Module[{pos = 0, weight = 1,
          nReals = 0, intervals = {},
          readBit, acceptTestI, makeRealVar, seqBernoulli, makeContToken,
          result},

    readBit[] := (
      pos++;
      If[pos > Length[bits], Throw[$OutOfBits, $dbc$tag]];
      weight *= (1/2); bits[[pos]]);

    acceptTestI[j_, p_] := Module[{lo, hi, pVal, condP, pR},
      lo   = intervals[[j,1]]; hi = intervals[[j,2]];
      pVal = If[NumericQ[p], p, PiecewiseExpand[p]];
      If[TrueQ[pVal <= lo], Return[False, Module]];
      If[TrueQ[pVal >= hi], Return[True,  Module]];
      condP = (pVal - lo) / (hi - lo);
      pos++;
      If[pos > Length[bits], Throw[$OutOfBits, $dbc$tag]];
      pR = condP /. {r_Real :> Rationalize[r]};
      If[bits[[pos]] == 1,
        weight *= pR;   intervals[[j]] = {lo, pVal}; True,
        weight *= (1-pR); intervals[[j]] = {pVal, hi}; False]];

    makeRealVar[] := (
      nReals++;
      AppendTo[intervals, {0,1}];
      $dbc$irand[nReals, acceptTestI]);

    seqBernoulli[ws_List, elems_List] :=
      Module[{n = Length[elems], remainW = Total[ws], chosen = Length[elems], var, p},
        Do[
          var = makeRealVar[];
          p   = ws[[i]] / remainW;
          If[var < p, chosen = i; Break[], remainW -= ws[[i]]],
          {i, 1, n-1}];
        elems[[chosen]]];

    makeContToken[lo_, hi_] := $dbc$contToken["Uniform", lo, hi, seqBernoulli];

    result = Catch[
      Block[{
        RandomReal = Function[Module[{args = {##}},
          Which[
            args === {},       makeRealVar[],
            MatchQ[args,{{_,_}}], makeContToken[args[[1,1]], args[[1,2]]],
            True, Throw[$dbc$cantHandle["RandomReal["<>ToString[args]<>"]: unsupported"],$dbc$tag]]]],
        Random = Function[{}, makeRealVar[]],
        RandomInteger = Function[Module[{args = {##}},
          Which[
            args==={}||args==={1}, readBit[],
            MatchQ[args,{{_Integer,_Integer}}],
              Module[{lo=args[[1,1]], hi=args[[1,2]], n, k, val},
                n=hi-lo+1; If[n==1,lo,
                  k=IntegerLength[n-1,2]; val=$dbc$readBitsAsInt[k,readBit];
                  If[val>=n, Throw[$dbc$outOfRange,$dbc$tag], lo+val]]],
            MatchQ[args,{_Integer?NonNegative}],
              Module[{n=args[[1]]+1, k, val},
                If[n==1,0,
                  k=IntegerLength[n-1,2]; val=$dbc$readBitsAsInt[k,readBit];
                  If[val>=n, Throw[$dbc$outOfRange,$dbc$tag], val]]],
            True, Throw[$dbc$cantHandle["RandomInteger["<>ToString[args]<>"]"],$dbc$tag]]]],
        RandomChoice = Function[Module[{args = {##}},
          Which[
            Length[args]==1 && !MatchQ[args[[1]],_Rule],
              Module[{list=args[[1]],n,k,idx},
                n=Length[list];
                Which[n==0, Throw[$dbc$cantHandle["RandomChoice: empty"],$dbc$tag],
                      n==1, list[[1]],
                      True, k=IntegerLength[n-1,2]; idx=$dbc$readBitsAsInt[k,readBit];
                            If[idx>=n,Throw[$dbc$outOfRange,$dbc$tag],list[[1+idx]]]]],
            Length[args]==1 && MatchQ[args[[1]],Rule[_List,_List]],
              Module[{ws=args[[1,1]], elems=args[[1,2]]},
                If[Length[ws]=!=Length[elems],
                  Throw[$dbc$cantHandle["RandomChoice w->e: length mismatch"],$dbc$tag]];
                If[Length[elems]==0,Throw[$dbc$cantHandle["RandomChoice w->e: empty"],$dbc$tag]];
                If[Length[elems]==1, elems[[1]], seqBernoulli[ws,elems]]],
            True, Throw[$dbc$cantHandle["RandomChoice: unsupported form"],$dbc$tag]]]],
        RandomVariate = Function[Module[{args = {##}},
          Which[
            MatchQ[args,{HoldPattern[UniformDistribution[{_,_}]]}],
              makeContToken[args[[1,1,1]], args[[1,1,2]]],
            MatchQ[args,{HoldPattern[NormalDistribution[_,_]]}],
              Module[{mu=args[[1,1]],sigma=args[[1,2]],nMax,vals,rawW},
                nMax=Floor[$dbcCurrentNGrid/2];
                vals=Range[Round[mu]-nMax,Round[mu]+nMax];
                rawW=Table[CDF[NormalDistribution[mu,sigma],k+1/2]-
                           CDF[NormalDistribution[mu,sigma],k-1/2],{k,vals}];
                seqBernoulli[rawW/Total[rawW],vals]],
            True, Throw[$dbc$cantHandle["RandomVariate: unsupported"],$dbc$tag]]]],
        RandomPermutation = Function[Module[{args={##},list,n,perm},
          Which[
            MatchQ[args,{_Integer?Positive}],
              n=args[[1]]; perm=Range[n];
              Do[With[{j=RandomInteger[{1,i}]},perm[[{i,j}]]=perm[[{j,i}]]],{i,n,2,-1}]; perm,
            MatchQ[args,{_List}],
              list=args[[1]]; n=Length[list]; If[n==0,Return[{},Module]];
              perm=Range[n];
              Do[With[{j=RandomInteger[{1,i}]},perm[[{i,j}]]=perm[[{j,i}]]],{i,n,2,-1}];
              list[[perm]],
            True, Throw[$dbc$cantHandle["RandomPermutation: unsupported"],$dbc$tag]]]],
        RandomSample = Function[Module[{args={##},list,k,n,perm},
          Which[
            MatchQ[args,{_List}],
              list=args[[1]]; n=Length[list]; If[n==0,Return[{},Module]];
              perm=Range[n];
              Do[With[{j=RandomInteger[{1,i}]},perm[[{i,j}]]=perm[[{j,i}]]],{i,n,2,-1}];
              list[[perm]],
            MatchQ[args,{_List,_Integer?NonNegative}],
              list=args[[1]]; k=args[[2]]; n=Length[list];
              If[k>n,Throw[$dbc$cantHandle["RandomSample: k>n"],$dbc$tag]];
              perm=Range[n];
              Table[With[{j=RandomInteger[{1,i}]},
                perm[[{i,j}]]=perm[[{j,i}]]; list[[perm[[i]]]]],
                {i,n,n-k+1,-1}],
            True, Throw[$dbc$cantHandle["RandomSample: unsupported"],$dbc$tag]]]]},
      alg[state]],
      $dbc$tag, Function[{ex}, ex]];

    Which[
      result === $OutOfBits,       $OutOfBits,
      result === $dbc$outOfRange,  $dbc$outOfRange,
      MatchQ[result,$dbc$cantHandle[_]], result,
      !FreeQ[{result,weight},$dbc$irand],
        $dbc$cantHandle["RandomReal[] was used in an unsupported way (not directly compared to a threshold)"],
      !FreeQ[{result,weight},$dbc$contToken],
        $dbc$cantHandle["RandomReal[{lo,hi}] was never discretised with Floor/Round/Ceiling"],
      True, {result, weight}]]


(* ----------------------------------------------------------------
   $dbcBuildStateLeaves
   BFS for a SINGLE state. Returns list of {bits, nextState, weight}
   leaves, or $dbc$cantHandle[msg].
   normFn: normalisation function applied to each returned state.
   ---------------------------------------------------------------- *)
$dbcBuildStateLeaves[state_, alg_, maxDepth_Integer, tlim_,
                     nGrid_Integer, normFn_] :=
  Module[{queue = {{}}, leaves = {}, t0 = AbsoluteTime[], timedOut = False,
          bits, res, ns, w},
    $dbcCurrentNGrid = nGrid;
    While[queue =!= {} && !timedOut,
      If[AbsoluteTime[]-t0 > tlim, timedOut = True; Break[]];
      bits = First[queue]; queue = Rest[queue];
      res  = RunWithBitsAT[alg, state, bits];
      Which[
        res === $OutOfBits && Length[bits] < maxDepth,
          queue = Join[queue, {Append[bits,0], Append[bits,1]}],
        res === $OutOfBits,
          Null,  (* path silently dropped at MaxBitDepth *)
        res === $dbc$outOfRange,
          Null,  (* rejection-sampled path, silently dropped *)
        MatchQ[res,$dbc$cantHandle[_]],
          Return[res, Module],
        True,
          {ns, w} = res;
          AppendTo[leaves, {bits, normFn[ns], w}]]];
    If[timedOut, Print["  WARNING: time limit for state ",state]];
    leaves]


(* ----------------------------------------------------------------
   BuildTreeAT
   State-space BFS from seedState. Discovers all reachable states.
   Returns Association[ state → {leaf,...} ] or $dbc$cantHandle[msg].
   normFn: normalisation applied to every returned state.
   ---------------------------------------------------------------- *)
Options[BuildTreeAT] = {"MaxBitDepth"->20, "TimeLimit"->120., "Verbose"->False}

BuildTreeAT[seedState_, alg_, nGrid_Integer, normFn_, OptionsPattern[]] :=
  Module[{maxDepth = OptionValue["MaxBitDepth"],
          tlim     = N @ OptionValue["TimeLimit"],
          verbose  = OptionValue["Verbose"],
          discovered, toProcess, batch, result, s, leaves, err},
    $dbcCurrentNGrid = nGrid;
    discovered = {seedState}; toProcess = {seedState}; result = <||>;
    While[toProcess =!= {},
      batch = toProcess; toProcess = {};
      If[verbose, Print["  BFS wave: ", Length[batch], " states"]];
      Do[
        s      = batch[[i]];
        leaves = $dbcBuildStateLeaves[s, alg, maxDepth, tlim, nGrid, normFn];
        If[MatchQ[leaves, $dbc$cantHandle[_]], Return[leaves, Module]];
        result[s] = leaves;
        Do[
          If[!MemberQ[discovered, leaf[[2]]],
            AppendTo[discovered, leaf[[2]]];
            AppendTo[toProcess,  leaf[[2]]]],
          {leaf, leaves}],
        {i, Length[batch]}]];
    result]

(* Build transition matrix from BuildTreeAT output. *)
TreeATToMatrix[treeData_Association] :=
  Module[{matrix = <||>},
    Do[
      Do[With[{ns=leaf[[2]],w=leaf[[3]]},
        matrix[{s,ns}] = Lookup[matrix,Key[{s,ns}],0] + w],
        {leaf, treeData[s]}],
      {s, Keys[treeData]}];
    matrix]


(* ================================================================
   SECTION 4 — TRANSLATIONAL INVARIANCE CHECK (τ BFS)
   ================================================================
   Run state-space BFS with τ-augmented positions.  τ cancels in all
   pairwise differences → transition probabilities are τ-free iff the
   algorithm is translation-invariant over all reachable states.
   Algebraically certified by a single FreeQ pass.
   ================================================================ *)

$dbcCheckTranslational[seedState_List, alg_, nGrid_Integer,
                        maxDepth_Integer, tlim_] :=
  Module[{τSeed, τNorm, τTree, τMatrix, tauFree, violations},
    τSeed = $dbcAddTau[seedState];
    τNorm = $dbcNormTauState[#, nGrid] &;
    τTree = BuildTreeAT[τSeed, alg, nGrid, τNorm,
                        "MaxBitDepth"->maxDepth, "TimeLimit"->tlim];
    If[MatchQ[τTree,$dbc$cantHandle[_]],
      Return[<|"tauFree"->$Failed,"error"->τTree[[1]]|>]];
    τMatrix   = TreeATToMatrix[τTree];
    tauFree   = AllTrue[Values[τMatrix], FreeQ[#,τr] && FreeQ[#,τc] &];
    violations = If[tauFree, {},
      Take[Select[Normal[τMatrix], !(FreeQ[#[[2]],τr]&&FreeQ[#[[2]],τc])&], UpTo[5]]];
    <|"tauFree"->tauFree, "violations"->violations|>]


(* ================================================================
   SECTION 5 — D4 GENERATOR TESTING
   ================================================================
   For each orbit rep s₀, run BFS from rot90(s₀) and reflect(s₀).
   Compare: T(rot90(s₀)→rot90(t)) === T(s₀→t)  [structural equality]
            T(reflect(s₀)→reflect(t)) === T(s₀→t)
   If all match: D4 invariance exactly verified.
   Returns <|"d4Pass"->True/False, "fails"->{...}|>.
   ================================================================ *)

$dbcVerifyD4[repLeaves_Association, alg_, nGrid_Integer,
              maxDepth_Integer, tlim_] :=
  Module[{normFn = $dbcNormState[#,nGrid]&,
          reps, d4Pass = True, fails = {},
          s0, repMatrix,
          rotState, rotLeaves, rotMatrix, rViol,
          refState, refLeaves, refMatrix, sViol},
    reps = Keys[repLeaves];
    Do[
      s0 = reps[[ri]];
      repMatrix = <||>;
      Do[With[{ns=leaf[[2]],w=leaf[[3]]},
        repMatrix[{s0,ns}]=Lookup[repMatrix,Key[{s0,ns}],0]+w],
        {leaf, repLeaves[s0]}];

      (* --- Rotation check --- *)
      rotState  = $dbcApplyRot90[s0, nGrid];
      rotLeaves = $dbcBuildStateLeaves[rotState,alg,maxDepth,tlim,nGrid,normFn];
      If[MatchQ[rotLeaves,$dbc$cantHandle[_]],
        AppendTo[fails,<|"rep"->s0,"gen"->"rot90","error"->rotLeaves[[1]]|>];
        d4Pass=False; Continue[]];
      rotMatrix = <||>;
      Do[With[{ns=leaf[[2]],w=leaf[[3]]},
        rotMatrix[{rotState,ns}]=Lookup[rotMatrix,Key[{rotState,ns}],0]+w],
        {leaf,rotLeaves}];
      rViol = None;
      KeyValueMap[Function[{pair,tVal},
        With[{dest=pair[[2]],rotDest=$dbcApplyRot90[pair[[2]],nGrid]},
          With[{tRot=Lookup[rotMatrix,Key[{rotState,rotDest}],0]},
            If[tVal=!=tRot, rViol=<|"rep"->s0,"dest"->dest,"T(rep)"->tVal,"T(rot)"->tRot|>]]]],
        repMatrix];
      If[rViol=!=None, AppendTo[fails,<|"gen"->"rot90","mismatch"->rViol|>]; d4Pass=False];

      (* --- Reflection check --- *)
      refState  = $dbcApplyReflect[s0, nGrid];
      refLeaves = $dbcBuildStateLeaves[refState,alg,maxDepth,tlim,nGrid,normFn];
      If[MatchQ[refLeaves,$dbc$cantHandle[_]],
        AppendTo[fails,<|"rep"->s0,"gen"->"reflect","error"->refLeaves[[1]]|>];
        d4Pass=False; Continue[]];
      refMatrix = <||>;
      Do[With[{ns=leaf[[2]],w=leaf[[3]]},
        refMatrix[{refState,ns}]=Lookup[refMatrix,Key[{refState,ns}],0]+w],
        {leaf,refLeaves}];
      sViol = None;
      KeyValueMap[Function[{pair,tVal},
        With[{dest=pair[[2]],refDest=$dbcApplyReflect[pair[[2]],nGrid]},
          With[{tRef=Lookup[refMatrix,Key[{refState,refDest}],0]},
            If[tVal=!=tRef, sViol=<|"rep"->s0,"dest"->dest,"T(rep)"->tVal,"T(ref)"->tRef|>]]]],
        repMatrix];
      If[sViol=!=None, AppendTo[fails,<|"gen"->"reflect","mismatch"->sViol|>]; d4Pass=False],

      {ri, Length[reps]}];
    <|"d4Pass"->d4Pass, "fails"->fails|>]


(* ================================================================
   SECTION 6 — ERGODICITY CHECK
   ================================================================ *)

CheckErgodicity[allStates_List, nGrid_Integer] :=
  Module[{s0, types, typeCounts, S, N, theoretical},
    s0          = First[allStates];
    types       = #[[2]] & /@ s0;
    typeCounts  = Values[Counts[types]];
    S           = nGrid^2;
    N           = Length[types];
    theoretical = Product[S-k,{k,0,N-1}] / Times@@(Factorial/@typeCounts);
    <|"ergodic"     -> (Length[allStates] == theoretical),
      "found"       -> Length[allStates],
      "theoretical" -> theoretical|>]


(* ================================================================
   SECTION 7 — SCHWARTZ-ZIPPEL DB CHECK
   ================================================================ *)

$szRandQ[] := RandomChoice[{-1,1}] * RandomInteger[{1,50}] / RandomInteger[{1,50}]

$dbcFS = Symbol["$dbcFallbackSentinel"];

$dbcMergeExp[expr_] := FixedPoint[Function[e,
  e /. {Power[E^x_,n_]:>E^Expand[n*x], Power[Exp[x_],n_]:>Exp[Expand[n*x]]}
    // Expand
    /. {Times[a___,E^x_,E^y_,b___]:>Times[a,E^Expand[x+y],b],
        Times[a___,Exp[x_],Exp[y_],b___]:>Times[a,Exp[Expand[x+y]],b]}],
  expr, 30]

$dbcSplitTerm[0]  := {0,0}; $dbcSplitTerm[0.] := {0,0}
$dbcSplitTerm[t_] := Module[{ep},
  ep = Cases[{t}, Exp[x_]:>x, {0,Infinity}];
  Which[Length[ep]==0,{t,0},
        Length[ep]==1,{Expand[Cancel[t/Exp[ep[[1]]]]],Expand[ep[[1]]]},
        True,$dbcFS]]

$dbcGroupByExp[splits_List] :=
  Module[{polys={},groups={},matched},
    Scan[Function[s,
      matched=SelectFirst[Range@Length[polys],
        Function[i,Expand[s[[2]]-polys[[i]]]===0],0];
      If[matched===0,
        AppendTo[polys,s[[2]]]; AppendTo[groups,{s[[1]]}],
        groups[[matched]]=Append[groups[[matched]],s[[1]]]]],
      splits];
    MapThread[{#1,#2}&,{polys,groups}]]

$dbcIsExpZero[0]    := True; $dbcIsExpZero[0.] := True
$dbcIsExpZero[expr_] :=
  Module[{merged,expanded,terms,splits,grouped},
    merged   = $dbcMergeExp[expr];
    expanded = Expand[merged];
    If[expanded===0,Return[True]];
    terms  = If[Head[expanded]===Plus, List@@expanded, {expanded}];
    splits = Map[$dbcSplitTerm,terms];
    If[MemberQ[splits,$dbcFS],Return[$dbcFS]];
    grouped = $dbcGroupByExp[splits];
    If[AllTrue[grouped,Function[g,Expand[Total[g[[2]]]]===0]],True,False]]

$dbcAbstractTrans[expr_] :=
  Module[{atoms,k=0,subs},
    atoms=DeleteDuplicates@Cases[expr,
      (Erf|Erfc|FresnelS|FresnelC|SinIntegral|CosIntegral|
       ExpIntegralEi|LogIntegral|BesselJ|BesselY|BesselI|BesselK)[__?NumericQ]|
      ExpIntegralE[_Integer,_?NumericQ],Infinity];
    If[atoms==={},Return[expr]];
    subs=Map[(#->Symbol["$dbcTr$"<>ToString[++k]])&,atoms];
    expr/.subs]

$dbcSZCheckOne[expr_, symParams_List, nReps_Integer] :=
  Module[{assign,subst,res,result=True},
    Do[
      assign=Map[#->$szRandQ[]&,symParams];
      subst =PiecewiseExpand[expr/.assign, \[Beta]>0];
      subst =$dbcAbstractTrans[subst];
      res   =$dbcIsExpZero[subst];
      Which[res===True,Null,
            res===False,result={False,subst,assign};Break[],
            True,       result=$dbcFS;Break[]],
      {i,nReps}];
    result]

Options[CheckDetailedBalanceSZ] = {
  "FullSimplify" -> False,
  "SZRepeats"   -> 30,
  "FailFast"    -> False}

CheckDetailedBalanceSZ[matrix_Association, allStates_List, symEnergy_,
                        extraAssumptions_List:{}, OptionsPattern[]] :=
  Module[{useFS, nReps, failFast, assm, symParams,
          stateToIdx, energyCache, pairs, violations,
          si, sj, tij, tji, ei, ej, expr, res},
    useFS    = TrueQ@OptionValue["FullSimplify"];
    nReps    = OptionValue["SZRepeats"];
    failFast = TrueQ@OptionValue["FailFast"];
    assm     = Join[{\[Beta]>0}, extraAssumptions];
    symParams= Cases[extraAssumptions, Element[x_,Reals]:>x, Infinity];

    stateToIdx  = AssociationThread[allStates->Range[Length[allStates]]];
    energyCache = AssociationThread[allStates->
      Map[symEnergy[#]/.r_Real:>Rationalize[r]&,allStates]];
    pairs = DeleteDuplicates[Sort/@Select[
      Map[{stateToIdx[#[[1]]],stateToIdx[#[[2]]]}&,Keys[matrix]],
      #[[1]]=!=#[[2]]&]];
    If[Length[pairs]==0, Return[{}]];

    violations={};
    Do[
      si=allStates[[pairs[[k,1]]]]; sj=allStates[[pairs[[k,2]]]];
      tij=Lookup[matrix,Key[{si,sj}],0]; tji=Lookup[matrix,Key[{sj,si}],0];
      ei=energyCache[si]; ej=energyCache[sj];
      expr=tij*Exp[-\[Beta]*ei]-tji*Exp[-\[Beta]*ej];
      If[expr=!=0,
        res=If[useFS,
          With[{fs=FullSimplify[PiecewiseExpand[expr],Assumptions->assm]},
            If[fs===0,True,False]],
          $dbcSZCheckOne[expr,symParams,nReps]];
        If[res=!=True,
          AppendTo[violations,<|"pair"->{si,sj},"tij"->tij,"tji"->tji,
                                 "ei"->ei,"ej"->ej,"expr"->expr|>];
          If[failFast, Return[violations,Module]]]],
      {k,Length[pairs]}];
    violations]


(* ================================================================
   SECTION 7b — DIRECT SZ CHECK FROM LEAVES (fast path)
   ================================================================
   Key insight: the leaf weight leaf[[3]] is identical for all 72
   G-orbit members of an orbit rep — only positions change, not
   probabilities.  Evaluate each weight ONCE, then scatter to all
   orbit members via a precomputed integer G-action table.

   Pre-computation (done once before the SZ loop):
     gActionTable[g, i] = j  (integer state indices)
   Inner SZ loop: only numeric adds — no applyGElem, no pattern match.
   ================================================================ *)

$dbcSZCheckLeaves[repLeaves_Association, repToOrbitMap_Association,
                   energy_, allStates_List, nGrid_Integer,
                   symParams_List, nReps_Integer,
                   tol_:1*^-7] :=
  Module[{violations = {}, betaVal, assign,
          nStates, stateToIdx, energyExprs,
          (* G-action integer table: gActI[[gIdx, stateIdx]] = newStateIdx *)
          allGKeys, gIdxOf, gActI,
          (* Precomputed per-rep orbit expansion as integer index pairs *)
          repList, repOrbitPairs,
          (* SZ loop temporaries *)
          wVals, rowT, pairs, pairsArr,
          betaArr, lhs, rhs, ei, ej},

    nStates    = Length[allStates];
    stateToIdx = AssociationThread[allStates -> Range[nStates]];

    (* Collect every G element referenced in any orbit map *)
    allGKeys = DeleteDuplicates @ Flatten[Values /@ Values[repToOrbitMap], 1];
    gIdxOf   = AssociationThread[allGKeys -> Range[Length[allGKeys]]];

    (* Build integer G-action table: gActI[[gi, si]] = target state index *)
    gActI = Table[0, {Length[allGKeys]}, {nStates}];
    Do[
      With[{gi = gIdxOf[g]},
        Do[
          gActI[[gi, si]] = stateToIdx[$dbcApplyGElem[g, allStates[[si]], nGrid]],
          {si, nStates}]],
      {g, allGKeys}];

    (* Pre-expand: for each rep, for each leaf, list of {srcIdx,tgtIdx} pairs *)
    repList = Keys[repLeaves];
    repOrbitPairs = Table[
      With[{repKey = repList[[ri]], orbitMap = repToOrbitMap[repList[[ri]]]},
        Table[
          With[{tgtIdx = stateToIdx[repLeaves[repKey][[li, 2]]]},
            Map[Function[sKey,
                {stateToIdx[sKey],
                 gActI[[gIdxOf[orbitMap[sKey]], tgtIdx]]}],
              Keys[orbitMap]]],
          {li, Length[repLeaves[repKey]]}]],
      {ri, Length[repList]}];

    (* Pre-compute symbolic energy for each state (for fast substitution) *)
    energyExprs = Map[energy, allStates];

    (* Precompute unique unordered pairs that can appear *)
    pairsArr = DeleteDuplicates @ Sort @ Map[Sort,
      Flatten[repOrbitPairs, 2]];
    pairsArr = Select[pairsArr, #[[1]] =!= #[[2]] &];

    (* ---- SZ loop ---- *)
    Do[
      betaVal = RandomReal[{0.3, 5.0}];
      assign  = Join[
        Map[# -> N[$szRandQ[]] &, symParams],
        {\[Beta] -> betaVal, numBeta -> betaVal}];

      (* Evaluate all leaf weights once per leaf (not per orbit member) *)
      wVals = Table[
        N[repLeaves[repList[[ri]]][[All, 3]] /. assign],
        {ri, Length[repList]}];

      (* Accumulate numerical T matrix using integer table *)
      rowT = SparseArray[{}, {nStates, nStates}];
      Do[
        Do[
          With[{w = wVals[[ri, li]]},
            If[NumericQ[w] && w != 0.,
              Scan[Function[st, rowT[[st[[1]], st[[2]]]] += w],
                   repOrbitPairs[[ri, li]]]]],
          {li, Length[repOrbitPairs[[ri]]]}],
        {ri, Length[repList]}];

      (* Energy vector for this coupling assignment *)
      betaArr = N[energyExprs /. assign];

      (* Detailed balance check *)
      Do[
        With[{i = pair[[1]], j = pair[[2]]},
          tij = rowT[[i, j]]; tji = rowT[[j, i]];
          ei  = betaArr[[i]]; ej  = betaArr[[j]];
          lhs = tij * Exp[-betaVal * ei];
          rhs = tji * Exp[-betaVal * ej];
          If[Abs[lhs - rhs] > tol * Max[Abs[lhs], Abs[rhs], 1*^-30],
            AppendTo[violations,
              <|"pair" -> {allStates[[i]], allStates[[j]]},
                "tij"  -> tij, "tji"  -> tji,
                "ei"   -> ei,  "ej"   -> ej,
                "beta" -> betaVal,
                "lhs"  -> lhs, "rhs"  -> rhs|>]]],
        {pair, pairsArr}],
    {nReps}];
    violations]


(* ================================================================
   SECTION 8 — NUMERICAL MCMC CHECK (optional)
   ================================================================ *)

Options[RunNumericalMCMCAT] = {"NSteps"->100000,"WarmupFrac"->0.1}

RunNumericalMCMCAT[allStates_List, alg_, numBeta_, nGrid_Integer, OptionsPattern[]] :=
  Module[{nSteps=OptionValue["NSteps"],
          nWarmup=Round[OptionValue["NSteps"]*OptionValue["WarmupFrac"]],
          state,counts},
    state=RandomChoice[allStates];
    counts=AssociationThread[allStates->0];
    Do[Block[{\[Beta]=numBeta,$nGrid=nGrid,$dbcCurrentNGrid=nGrid},
      state=alg[state]],{nWarmup}];
    Do[Block[{\[Beta]=numBeta,$nGrid=nGrid,$dbcCurrentNGrid=nGrid},
      state=alg[state]];
      If[KeyExistsQ[counts,state],counts[state]++],
      {nSteps-nWarmup}];
    counts]

BoltzmannWeightsAT[allStates_List,energy_,numBeta_] :=
  Module[{ws,Z},
    ws=N[Exp[-numBeta*energy[#]]&/@allStates]; Z=Total[ws];
    AssociationThread[allStates->ws/Z]]
