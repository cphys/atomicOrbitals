(* ::Package:: *)

ClearAll["Global`*"]
<<"hydrogenAtom`"
ParallelNeeds["hydrogenAtom`"]

numb=40;
(* Import NIST data *)
nistDataNR=Import[FileNameJoin[{NotebookDirectory[],"inputData","nistDataNR.csv"}],"CSV"];
nistDataNR=Take[nistDataNR,numb-1];
nistData=Import[FileNameJoin[{NotebookDirectory[],"inputData","nistData.csv"}],"CSV"];

(* Here all oscillator strengths are calculated 
   (the rest is just formatting) *)
valuesNR = ParallelTable[hydrogenOscillatorsNR[n], {n, 2, numb}];
values = Flatten[
   ParallelTable[
    Prepend[#, ToString[n]] & /@ 
     MapThread[
      Prepend, {hydrogenOscillators[
        n], {"\!\(\*FractionBox[\(1\), \(2\)]\)", 
        "\!\(\*FractionBox[\(3\), \(2\)]\)"}}], {n, 2, 6}], 1];
        
(* Rounding and formatting *)
valuesFNR = Table[{
    ToString[i + 1],
    NumberForm[valuesNR[[i]][[1]], {6, 6}],
    NumberForm[valuesNR[[i]][[2]], {6, 6}]},
   {i, Length[valuesNR]}];

valuesF = Table[{
    values[[i]][[1]],
    values[[i]][[2]],
    NumberForm[values[[i]][[3]], {6, 6}],
    NumberForm[values[[i]][[4]], {5, 5}]},
   {i, Length[values]}];
   
pDiff[f1_, f2_] := NumberForm[Abs[f2 - f1]/f1*100, {7, 5}]

(* Calculate % differences between our calculations and NIST data *)
pdifsNR =
  Table[
   {pDiff[valuesNR[[i]][[1]], nistDataNR[[i]][[1]]], 
    pDiff[valuesNR[[i]][[2]], nistDataNR[[i]][[2]]]}, {i, 
    Length[valuesNR]}];
pdifs = Table[
   {pDiff[values[[i]][[-2]], nistData[[i]][[1]]], 
    pDiff[values[[i]][[-1]], nistData[[i]][[2]]]},
   {i, Length[values]}];
   
appendF[m1_,m2_]:=FlattenAt[#,-1]&/@MapThread[Append,{m1,m2}]

Export[
FileNameJoin[{NotebookDirectory[],"gridNR.png"}],
Grid[Prepend[
  Fold[appendF, valuesFNR, {nistDataNR, pdifsNR}],
  {"n", "\!\(\*SubscriptBox[\(E\), \(n0\)]\)", 
   "\!\(\*SubscriptBox[\(f\), \(n0\)]\)", 
   "\!\(\*SubscriptBox[\(E\), \(n0\)]\) (NIST)", 
   "\!\(\*SubscriptBox[\(f\), \(n0\)]\)(NIST)", 
   "% diff \!\(\*SubscriptBox[\(E\), \(n0\)]\)", 
   "% diff \!\(\*SubscriptBox[\(f\), \(n0\)]\)"}],
 Dividers -> All]];

Export[
FileNameJoin[{NotebookDirectory[],"grid.png"}],
Grid[Prepend[
  Fold[appendF, valuesF, {nistData, pdifs}],
  {"n", "j", "\!\(\*SubscriptBox[\(E\), \(n0\)]\)", 
   "\!\(\*SubscriptBox[\(f\), \(n0\)]\)", 
   "\!\(\*SubscriptBox[\(E\), \(n0\)]\) (NIST)", 
   "\!\(\*SubscriptBox[\(f\), \(n0\)]\)(NIST)", 
   "% diff \!\(\*SubscriptBox[\(E\), \(n0\)]\)", 
   "% diff \!\(\*SubscriptBox[\(f\), \(n0\)]\)"}],
 Dividers -> All]];
 
CloseKernels[];
Exit[];


(* ::Input:: *)
(**)
