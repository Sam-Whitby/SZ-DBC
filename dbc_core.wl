(* ================================================================
   SZ-DBC — Detailed Balance Checker (core library)
   ================================================================
   State format: sorted list of {{row,col}, type} pairs.
   Positions are integers in {1,..,nGrid}.

   Key features:
   - COUPLING SYMMETRY: automatically canonicalises coupling atoms
     so that couplingJ[b,a,d] = couplingJ[a,b,d] for b>a.
   - TRANSLATIONAL INVARIANCE verified algebraically via τ-BFS:
     runs from orbit representatives (not full state space) so
     that the τ check costs ~2 s instead of ~135 s.
   - D4 ROTATIONAL INVARIANCE verified numerically (SZ-style) at
     random parameter points — more robust than structural equality.
   - ERGODICITY verified via reachability BFS on the orbit-expanded
     transition graph derived from the orbit-rep BFS output.
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
   SECTION 0b — COUPLING SYMMETRY
   ================================================================
   Adds a canonicalisation rule for every coupling head found in
   symCouplings: coupling[b,a,...] with b>a → coupling[a,b,...].
   This ensures that all leaf weights and energy expressions use
   the canonical (a ≤ b) form, so the SZ substitution covers all
   atoms that actually appear.
   Call AFTER ClearAll on the coupling head and BEFORE any BFS.
   ================================================================ *)

$dbcAddCouplingSymmetry[symCouplings_List] :=
  Module[{heads = DeleteDuplicates[Head /@ symCouplings]},
    Do[
      With[{h = heads[[i]]},
        h[b_Integer, a_Integer, rest__Integer] /; b > a := h[a, b, rest]],
      {i, Length[heads]}]]


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
    rotIdxs = Which[
      MemberQ[symGroup, "D4"], Range[0, 7],
      MemberQ[symGroup, "C4"], Range[0, 3],
      True, {0}];
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
      (* Build orbit map using group elements relative to the CANONICAL
         REP (minImg), not relative to s.  When s ≠ minImg (always true
         for D4 since the canonical rep may differ from the first
         unvisited state), storing g_s where g_s·s = img would cause the
         orbit expansion to apply g_s to leaf states from the BFS of
         minImg, producing wrong target positions.  We need h where
         h·minImg = img.  Recompute imgPairs from minImg to get these. *)
      With[{imgPairsFromRep =
              {$dbcApplyGElem[#, minImg, nGrid], #} & /@ allGElems},
        orbitMap = <||>;
        Scan[Function[ip,
            With[{img=ip[[1]], g=ip[[2]]},
              If[KeyExistsQ[stateSet,img] && !KeyExistsQ[orbitMap,img],
                orbitMap[img] = g]]],
          imgPairsFromRep]];
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
                  If[val>=n, Throw[$dbc$outOfRange,$dbc$tag], weight*=2^k/n; lo+val]]],
            MatchQ[args,{_Integer?NonNegative}],
              Module[{n=args[[1]]+1, k, val},
                If[n==1,0,
                  k=IntegerLength[n-1,2]; val=$dbc$readBitsAsInt[k,readBit];
                  If[val>=n, Throw[$dbc$outOfRange,$dbc$tag], weight*=2^k/n; val]]],
            True, Throw[$dbc$cantHandle["RandomInteger["<>ToString[args]<>"]"],$dbc$tag]]]],
        RandomChoice = Function[Module[{args = {##}},
          Which[
            Length[args]==1 && !MatchQ[args[[1]],_Rule],
              Module[{list=args[[1]],n,k,idx},
                n=Length[list];
                Which[n==0, Throw[$dbc$cantHandle["RandomChoice: empty"],$dbc$tag],
                      n==1, list[[1]],
                      True, k=IntegerLength[n-1,2]; idx=$dbc$readBitsAsInt[k,readBit];
                            If[idx>=n,Throw[$dbc$outOfRange,$dbc$tag],weight*=2^k/n; list[[1+idx]]]]],
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
   Warns if any paths are silently dropped at maxDepth.
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
          (* A BFS path reached maxDepth without the algorithm terminating.
             This means transitions exist that are not captured — the DB check
             would silently operate on an incomplete transition matrix.
             This is a hard error: increase -maxDepth to cover all paths. *)
          Return[$dbc$cantHandle[
            "BFS incomplete: a random-number path reached maxDepth=" <>
            ToString[maxDepth] <> " before the algorithm returned a state. " <>
            "Increase -maxDepth (current: " <> ToString[maxDepth] <>
            ") to ensure all transitions are covered."], Module],
        res === $dbc$outOfRange,
          Null,
        MatchQ[res,$dbc$cantHandle[_]],
          Return[res, Module],
        True,
          {ns, w} = res;
          AppendTo[leaves, {bits, normFn[ns], w}]]];
    If[timedOut, Print["  WARNING: time limit reached for state ", state]];
    leaves]




(* ================================================================
   SECTION 4 — TRANSLATIONAL INVARIANCE CHECK (τ BFS)
   ================================================================
   Run $dbcBuildStateLeaves with τ-augmented positions for each
   supplied state (typically the orbit representatives).
   τ cancels in all pairwise differences → transition probabilities
   are τ-free iff the algorithm is translation-invariant.
   Checking from orbit reps is ~63× faster than the full state-space
   BFS while covering all distinct states (given G-symmetry).
   ================================================================ *)

$dbcCheckTranslational[states_List, alg_, nGrid_Integer,
                        maxDepth_Integer, tlim_] :=
  Module[{τNorm, tauFree = True, violations = {}},
    τNorm = $dbcNormTauState[#, nGrid] &;
    Do[
      With[{τLeaves = $dbcBuildStateLeaves[
              $dbcAddTau[states[[si]]], alg, maxDepth, tlim, nGrid, τNorm]},
        If[MatchQ[τLeaves, $dbc$cantHandle[_]],
          Return[<|"tauFree"->$Failed, "error"->τLeaves[[1]]|>, Module]];
        Scan[Function[leaf,
          If[!FreeQ[leaf[[3]], τr] || !FreeQ[leaf[[3]], τc],
            tauFree = False;
            AppendTo[violations, leaf[[3]]]]],
          τLeaves]],
      {si, Length[states]}];
    <|"tauFree" -> tauFree, "violations" -> Take[violations, UpTo[5]]|>]






(* ================================================================
   SECTION 6 — ERGODICITY CHECKS
   ================================================================ *)

(* Combinatorial state-count check: verifies that $nGrid and
   $particleTypes are consistent (always passes for valid inputs
   but catches misconfiguration).  Separate from reachability. *)
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

(* Reachability ergodicity check using orbit-BFS data.
   Expands orbit-rep leaves via G-action to get all transitions,
   then BFS from state 1 to check whether every state is reachable.
   This tests the algorithm's actual connectivity, not just the
   state count.  Requires orbit BFS to have been completed. *)
$dbcCheckErgodicityFromLeaves[repLeaves_Association, repToOrbitMap_Association,
                               allStates_List, nGrid_Integer, seedState_List] :=
  Module[{nStates, stateToIdx, adjFwd, repList, repKey, orbitMap,
          srcIdx, dstIdx, reachableSet, queue, curIdx, seedIdx},
    nStates    = Length[allStates];
    stateToIdx = AssociationThread[allStates -> Range[nStates]];
    repList    = Keys[repLeaves];

    (* Build forward adjacency via explicit loops — avoids nested Join/Apply *)
    adjFwd = Association[Table[i -> {}, {i, nStates}]];
    Do[
      repKey   = repList[[ri]];
      orbitMap = repToOrbitMap[repKey];
      KeyValueMap[
        Function[{sKey, g},
          srcIdx = stateToIdx[sKey];
          Do[
            dstIdx = Lookup[stateToIdx,
                       Key[$dbcApplyGElem[g, repLeaves[repKey][[li, 2]], nGrid]], 0];
            If[dstIdx =!= 0 && dstIdx =!= srcIdx,
              adjFwd[srcIdx] = Union[adjFwd[srcIdx], {dstIdx}]],
            {li, Length[repLeaves[repKey]]}]],
        orbitMap],
      {ri, Length[repList]}];

    (* BFS from the actual seed state *)
    seedIdx = Lookup[stateToIdx, Key[seedState], 1];
    reachableSet = <|seedIdx -> True|>;
    queue = {seedIdx};
    While[queue =!= {},
      curIdx = First[queue]; queue = Rest[queue];
      Scan[Function[n,
        If[!KeyExistsQ[reachableSet, n],
          reachableSet[n] = True;
          AppendTo[queue, n]]],
        adjFwd[curIdx]]];

    <|"ergodic" -> (Length[reachableSet] == nStates),
      "reached" -> Length[reachableSet],
      "total"   -> nStates|>]


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
   Key insight: the leaf weight leaf[[3]] is identical for all |G|
   orbit members of an orbit rep — only positions change, not
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
  Module[{violations = <||>, done = False, betaVal, assign,
          nStates, stateToIdx, energyExprs,
          allGKeys, gIdxOf, gActI,
          repList, repOrbitPairs,
          wVals, rowT, pairsArr,
          betaArr, tij, tji, lhs, rhs, ei, ej},

    nStates    = Length[allStates];
    stateToIdx = AssociationThread[allStates -> Range[nStates]];

    allGKeys = DeleteDuplicates @ Flatten[Values /@ Values[repToOrbitMap], 1];
    gIdxOf   = AssociationThread[allGKeys -> Range[Length[allGKeys]]];

    gActI = Table[0, {Length[allGKeys]}, {nStates}];
    Do[
      With[{gi = gIdxOf[g]},
        Do[
          gActI[[gi, si]] = stateToIdx[$dbcApplyGElem[g, allStates[[si]], nGrid]],
          {si, nStates}]],
      {g, allGKeys}];

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

    energyExprs = Map[energy, allStates];

    pairsArr = DeleteDuplicates @ Sort @ Map[Sort,
      Flatten[repOrbitPairs, 2]];
    pairsArr = Select[pairsArr, #[[1]] =!= #[[2]] &];

    Do[
      If[done, Break[]];
      betaVal = RandomReal[{0.3, 5.0}];
      assign  = Join[
        Map[# -> N[$szRandQ[]] &, symParams],
        {\[Beta] -> betaVal, numBeta -> betaVal}];

      wVals = Table[
        N[repLeaves[repList[[ri]]][[All, 3]] /. assign],
        {ri, Length[repList]}];

      rowT = SparseArray[{}, {nStates, nStates}];
      Do[
        Do[
          With[{w = wVals[[ri, li]]},
            If[NumericQ[w] && w != 0.,
              Scan[Function[st, rowT[[st[[1]], st[[2]]]] += w],
                   repOrbitPairs[[ri, li]]]]],
          {li, Length[repOrbitPairs[[ri]]]}],
        {ri, Length[repList]}];

      betaArr = N[energyExprs /. assign];

      Do[
        With[{i = pair[[1]], j = pair[[2]]},
          tij = rowT[[i, j]]; tji = rowT[[j, i]];
          ei  = betaArr[[i]]; ej  = betaArr[[j]];
          If[NumericQ[ei] && NumericQ[ej],
            lhs = tij * Exp[-betaVal * ei];
            rhs = tji * Exp[-betaVal * ej];
            If[Abs[lhs - rhs] > tol * Max[Abs[lhs], Abs[rhs], 1*^-30],
              If[!KeyExistsQ[violations, pair],
                violations[pair] = <|"pair" -> {allStates[[i]], allStates[[j]]},
                                     "tij"  -> tij, "tji"  -> tji,
                                     "ei"   -> ei,  "ej"   -> ej,
                                     "beta" -> betaVal|>];
              If[Length[violations] >= 10, done = True; Break[]]]]],
        {pair, pairsArr}],
    {nReps}];
    Values[violations]]


(* ================================================================
   SECTION 7c — EBE (EXHAUSTIVE BRANCH ENUMERATION) DB CHECK
   ================================================================
   Exact DB verification by covering every coupling-constant branch
   region in the symbolic leaf weights.

   Phase 1 — Condition extraction (once):
     Scan all symbolic leaf weights and collect the k distinct
     Piecewise branch conditions.  Each condition is a comparison
     involving coupling-constant atoms (e.g. J_{12,1}-J_{12,2} < 0).

   Phase 2 — Region enumeration (once):
     Enumerate all 2^k sign patterns σ ∈ {0,1}^k.  For each, solve
     a small LP (FindInstance over Q) to find a rational coupling
     point J*(σ) in that region, or declare infeasible.

   Phase 3 — Exact DB check per feasible region:
     For each feasible region, substitute the rational J* into every
     leaf weight.  Piecewise conditions are now concrete True/False
     rationals, so every weight collapses to 0 or a concrete
     rational x Exp[-β x rational] expression.  β stays symbolic.
     Assemble T, form the DB expression for each communicating pair,
     and check exactly with $dbcIsExpZero — no floating point.

   Falls back to $dbcSZCheckLeaves when k > ebeMaxK.
   ================================================================ *)

$dbcEBECheckLeaves[repLeaves_Association, repToOrbitMap_Association,
                   energy_, allStates_List, nGrid_Integer,
                   allSymParams_List,
                   ebeMaxK_Integer : 50,
                   szFallbackReps_Integer : 30,
                   tol_ : 1*^-7] :=
  Module[{nStates, stateToIdx, energyExprs,
          allGKeys, gIdxOf, gActI, repList, repOrbitPairs, pairsArr,
          pairToLeafSrc,
          allLeaves, allConds, k,
          feasibleRegions, sigma, constraints, jStar,
          bfsQueue, visitedSigmas, sigma2, constraints2, jStar2,
          violations, done,
          jStarAssign, wVals, leafSrcs, tij, tji, ei, ej, expr, res},

    (* ---- Pre-compute G-action integer table (same as $dbcSZCheckLeaves) ---- *)
    nStates    = Length[allStates];
    stateToIdx = AssociationThread[allStates -> Range[nStates]];

    allGKeys = DeleteDuplicates @ Flatten[Values /@ Values[repToOrbitMap], 1];
    gIdxOf   = AssociationThread[allGKeys -> Range[Length[allGKeys]]];

    gActI = Table[0, {Length[allGKeys]}, {nStates}];
    Do[
      With[{gi = gIdxOf[g]},
        Do[gActI[[gi, si]] =
             stateToIdx[$dbcApplyGElem[g, allStates[[si]], nGrid]],
           {si, nStates}]],
      {g, allGKeys}];

    repList = Keys[repLeaves];
    repOrbitPairs = Table[
      With[{repKey = repList[[ri]],
            orbitMap = repToOrbitMap[repList[[ri]]]},
        Table[
          With[{tgtIdx = stateToIdx[repLeaves[repKey][[li, 2]]]},
            Map[Function[sKey,
                {stateToIdx[sKey],
                 gActI[[gIdxOf[orbitMap[sKey]], tgtIdx]]}],
              Keys[orbitMap]]],
          {li, Length[repLeaves[repKey]]}]],
      {ri, Length[repList]}];

    energyExprs = Map[energy, allStates];

    pairsArr = DeleteDuplicates @ Sort @ Map[Sort,
      Flatten[repOrbitPairs, 2]];
    pairsArr = Select[pairsArr, #[[1]] =!= #[[2]] &];

    (* Pre-compute: for each directed pair {src,tgt}, which (ri,li)
       orbit-expanded leaves contribute to T(src->tgt)?
       Stored as pairToLeafSrc[{src,tgt}] = {{ri,li}, ...}          *)
    pairToLeafSrc = <||>;
    Do[
      Scan[Function[st,
        With[{key = {st[[1]], st[[2]]}},
          pairToLeafSrc[key] =
            Append[Lookup[pairToLeafSrc, Key[key], {}], {ri, li}]]],
        repOrbitPairs[[ri, li]]],
      {ri, Length[repList]},
      {li, Length[repOrbitPairs[[ri]]]}];

    (* ============================================================
       PHASE 1 — Extract distinct Piecewise conditions
       ============================================================ *)
    allLeaves = Flatten[Values[repLeaves], 1];
    allConds = DeleteDuplicates @ Select[
      Flatten @ Map[Function[lf,
        (* HoldPattern prevents Mathematica evaluating the Piecewise LHS
           during matching (which would trigger Piecewise::pairs warnings
           for non-standard clause structures).  We only extract the
           second element of each 2-element clause.                      *)
        Cases[lf[[3]],
              HoldPattern[Piecewise[cl_List, _]] :>
                (#[[2]] & /@ Select[cl, Length[#] == 2 &]),
              Infinity]],
        allLeaves],
      Function[c,
        (* Keep only coupling-constant comparisons; discard True, False,
           === Infinity sentinels, and any non-comparison expressions.  *)
        MatchQ[c, _Less | _LessEqual | _Greater | _GreaterEqual] &&
        !FreeQ[c, Alternatives @@ allSymParams]]];
    k = Length[allConds];

    Print["  EBE: k=", k, " branch condition(s)"];

    (* ---- Fallback when k is too large ---- *)
    If[k > ebeMaxK,
      Print["  EBE: k=", k, " exceeds ebeMaxK=", ebeMaxK,
            " — falling back to probabilistic SZ (",
            szFallbackReps, " random points)"];
      Return[$dbcSZCheckLeaves[
               repLeaves, repToOrbitMap, energy,
               allStates, nGrid, allSymParams,
               szFallbackReps, tol]]];

    (* ============================================================
       PHASE 2 — Region enumeration via arrangement BFS
       ============================================================
       Replaces exhaustive 2^k enumeration with O(T×k) BFS over the
       hyperplane arrangement adjacency graph, where T is the number
       of real cells (feasible regions).

       Key guarantee: the d-dimensional cells of any hyperplane
       arrangement are connected under facet-adjacency (cells sharing a
       codimension-1 face differ in exactly one sign bit).  Therefore
       BFS from any starting cell reaches ALL cells, visiting each
       exactly once.  Total FindInstance calls <= T x k.
       ============================================================ *)
    feasibleRegions = {};
    If[k == 0 || Length[allSymParams] == 0,
      (* No Piecewise conditions: single trivial region *)
      feasibleRegions = {{Table[1, {k}], {}}};
      Print["  EBE: 1 feasible region (no branch conditions)"],
      (* BFS over the hyperplane arrangement *)
      jStar = Quiet[FindInstance[True, allSymParams, Rationals, 1]];
      If[jStar === {} || jStar === {{}},
        (* Fallback: shouldn't happen for real coupling constants *)
        feasibleRegions = {{Table[1, {k}], {}}};
        Print["  EBE: 1 feasible region (fallback — no rational start found)"],
        jStar = First[jStar];
        sigma = Table[If[TrueQ[allConds[[i]] /. jStar], 1, 0], {i, k}];
        visitedSigmas = <|sigma -> jStar|>;
        feasibleRegions = {{sigma, jStar}};
        bfsQueue = {sigma};
        While[bfsQueue =!= {},
          sigma = First[bfsQueue]; bfsQueue = Rest[bfsQueue];
          Do[
            sigma2 = ReplacePart[sigma, i -> 1 - sigma[[i]]];
            If[!KeyExistsQ[visitedSigmas, sigma2],
              constraints2 = And @@ Table[
                If[sigma2[[j]] == 1, allConds[[j]], !allConds[[j]]], {j, k}];
              jStar2 = Quiet[FindInstance[constraints2, allSymParams, Rationals, 1]];
              If[jStar2 =!= {} && jStar2 =!= {{}},
                visitedSigmas[sigma2] = First[jStar2];
                AppendTo[feasibleRegions, {sigma2, First[jStar2]}];
                AppendTo[bfsQueue, sigma2],
                visitedSigmas[sigma2] = None]],
            {i, k}]];
        Print["  EBE: ", Length[feasibleRegions], " feasible region(s)  (",
              Length[visitedSigmas], " sign patterns visited via BFS)"]]];

    (* ============================================================
       PHASE 3 — Exact DB check for each feasible region
       ============================================================ *)
    violations = <||>;
    done       = False;

    Do[
      If[done, Break[]];
      {sigma, jStarAssign} = feasibleRegions[[ri]];

      (* Substitute rational J* into every leaf weight.
         Piecewise conditions (now rational comparisons) auto-resolve:
           Piecewise[{{val, True}}, 0]  →  val
           Piecewise[{{val, False}}, 0] →  0
         β remains symbolic throughout.                              *)
      wVals = Table[
        repLeaves[repList[[rj]]][[All, 3]] /. jStarAssign,
        {rj, Length[repList]}];

      (* Energy values at J* — rational, β-free *)
      With[{eArr = energyExprs /. jStarAssign},

        Do[
          With[{i = pair[[1]], j = pair[[2]]},

            leafSrcs = Lookup[pairToLeafSrc, Key[{i, j}], {}];
            tij = If[leafSrcs === {}, 0,
                     Total[wVals[[#[[1]], #[[2]]]] & /@ leafSrcs]];

            leafSrcs = Lookup[pairToLeafSrc, Key[{j, i}], {}];
            tji = If[leafSrcs === {}, 0,
                     Total[wVals[[#[[1]], #[[2]]]] & /@ leafSrcs]];

            (* Skip pairs with no transitions in this region *)
            If[tij === 0 && tji === 0, Continue[]];

            ei = eArr[[i]]; ej = eArr[[j]];

            (* DB expression: T(i→j)·exp(-βE(i)) - T(j→i)·exp(-βE(j))
               For a correct algorithm this is algebraically zero.    *)
            expr = Expand[
              tij * Exp[-\[Beta] * ei] - tji * Exp[-\[Beta] * ej]];
            If[expr === 0, Continue[]];

            res = $dbcIsExpZero[expr];
            If[res =!= True,
              If[!KeyExistsQ[violations, pair],
                violations[pair] =
                  <|"pair"   -> {allStates[[i]], allStates[[j]]},
                    "tij"    -> tij,    "tji"   -> tji,
                    "ei"     -> ei,     "ej"    -> ej,
                    "region" -> sigma,
                    "jStar"  -> jStarAssign|>];
              If[Length[violations] >= 10,
                done = True; Break[]]]],

          {pair, pairsArr}]],

      {ri, Length[feasibleRegions]}];

    Values[violations]]


(* ================================================================
   SECTION 7d — ROTATIONAL INVARIANCE CHECK (numerical)
   ================================================================
   Checks T(s→t) = T(R·s → R·t) for all communicating pairs and
   each D4 element R in rotIdxsToCheck, using nReps random coupling
   points (SZ-style).

   rotIdxsToCheck: subset of {1..7}
     C4 rotations:  1 (rot90CW), 2 (rot180), 3 (rot270CW)
     Reflections:   4 (LR), 5 (UD), 6 (main diag), 7 (anti-diag)

   repLeaves/repToOrbitMap should be the TRANSLATION-ONLY BFS output
   (not yet reduced by any rotational symmetry).  The function builds
   the full T matrix by orbit expansion and checks the symmetry.

   Returns <|"pass" -> True/False, "violations" -> {...}|>.
   ================================================================ *)

$dbcCheckRotational[repLeaves_Association, repToOrbitMap_Association,
                     allStates_List, nGrid_Integer, symParams_List,
                     rotIdxsToCheck_List, nReps_Integer : 5,
                     tol_ : 1*^-6] :=
  Module[{nStates, stateToIdx, allGKeys, gIdxOf, gActI,
          repList, repOrbitPairs, pairsArr,
          rotAct, betaVal, assign, wVals, rowT,
          violations = {}, done = False,
          rsi, rsj, diffFwd, diffRev},

    nStates    = Length[allStates];
    stateToIdx = AssociationThread[allStates -> Range[nStates]];

    (* G-action table from the translation-only orbit map *)
    allGKeys = DeleteDuplicates @ Flatten[Values /@ Values[repToOrbitMap], 1];
    gIdxOf   = AssociationThread[allGKeys -> Range[Length[allGKeys]]];
    gActI    = Table[0, {Length[allGKeys]}, {nStates}];
    Do[
      With[{gi = gIdxOf[g]},
        Do[gActI[[gi, si]] =
             stateToIdx[$dbcApplyGElem[g, allStates[[si]], nGrid]],
           {si, nStates}]],
      {g, allGKeys}];

    repList = Keys[repLeaves];
    repOrbitPairs = Table[
      With[{repKey = repList[[ri]],
            orbitMap = repToOrbitMap[repList[[ri]]]},
        Table[
          With[{tgtIdx = stateToIdx[repLeaves[repKey][[li, 2]]]},
            Map[Function[sKey,
                {stateToIdx[sKey],
                 gActI[[gIdxOf[orbitMap[sKey]], tgtIdx]]}],
              Keys[orbitMap]]],
          {li, Length[repLeaves[repKey]]}]],
      {ri, Length[repList]}];

    pairsArr = DeleteDuplicates @ Sort @ Map[Sort,
      Flatten[repOrbitPairs, 2]];
    pairsArr = Select[pairsArr, #[[1]] =!= #[[2]] &];

    (* Rotation action table: rotAct[[ri, si]] = index of R_ri(allStates[[si]]) *)
    rotAct = Table[
      With[{g = <|"rot" -> rotIdxsToCheck[[ri]], "dr" -> 0, "dc" -> 0|>},
        Table[stateToIdx[$dbcApplyGElem[g, allStates[[si]], nGrid]],
              {si, nStates}]],
      {ri, Length[rotIdxsToCheck]}];

    Do[
      If[done, Break[]];
      betaVal = RandomReal[{0.3, 5.0}];
      assign  = Join[
        Map[# -> N[$szRandQ[]] &, symParams],
        {\[Beta] -> betaVal, numBeta -> betaVal}];

      wVals = Table[
        N[repLeaves[repList[[ri]]][[All, 3]] /. assign],
        {ri, Length[repList]}];

      rowT = SparseArray[{}, {nStates, nStates}];
      Do[
        Do[
          With[{w = wVals[[ri, li]]},
            If[NumericQ[w] && w != 0.,
              Scan[Function[st, rowT[[st[[1]], st[[2]]]] += w],
                   repOrbitPairs[[ri, li]]]]],
          {li, Length[repOrbitPairs[[ri]]]}],
        {ri, Length[repList]}];

      (* For each communicating pair and each rotation, check symmetry *)
      Do[
        If[done, Break[]];
        With[{i = pair[[1]], j = pair[[2]]},
          Do[
            rsi = rotAct[[ri, i]]; rsj = rotAct[[ri, j]];
            diffFwd = Abs[rowT[[i, j]] - rowT[[rsi, rsj]]];
            diffRev = Abs[rowT[[j, i]] - rowT[[rsj, rsi]]];
            If[diffFwd > tol * Max[Abs[rowT[[i, j]]], Abs[rowT[[rsi, rsj]]], 1*^-30] ||
               diffRev > tol * Max[Abs[rowT[[j, i]]], Abs[rowT[[rsj, rsi]]], 1*^-30],
              AppendTo[violations,
                <|"rot"   -> rotIdxsToCheck[[ri]],
                  "s"     -> allStates[[i]],
                  "t"     -> allStates[[j]],
                  "Tst"   -> rowT[[i, j]],
                  "TRsRt" -> rowT[[rsi, rsj]]|>];
              If[Length[violations] >= 5, done = True; Break[]]],
            {ri, Length[rotIdxsToCheck]}]],
        {pair, pairsArr}],
    {nReps}];

    <|"pass" -> (violations === {}), "violations" -> violations|>]


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
