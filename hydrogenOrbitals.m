(* ::Package:: *)

ClearAll["Global`*"]
<<"hydrogenAtom`"
ParallelNeeds["hydrogenAtom`"]


(* Change fastPlot to True to produce low resolution plots
Here there may be room for a medium setting
depending on the system runnning the process
currently running at highest setting (pp=50)
requires ~40 gigs of Ram, changing to a pp to a
lower number while keeping perf = "Quality" will
allow for lower ram usage *)
fastPlot = True;
If[fastPlot,
  {perf = "Speed", (* performance of density plot *)
   pp = 15}, (* number of plot points used in density plot *)
  {perf = "Quality",
   pp = 50}];
   
(* Animation will run 0 \[Rule] 2pi animationIncrements
controles the radians per frame for the output *)   
animationIncrements = \[Pi]/8;

(* True will save animated table *)
saveTableAnimation = True;
(* True will save a single non animated table *)
saveSingleTable = True;
(* True will save single orbital *)
saveSingleOrbital = True;

(* Name of output file *)
tableAnimationFName = "orbitalTableAnimation.gif";

(* viewpoint about the z-axis for output file *)
tableViewPoint = \[Pi]/4;
(* Name of table output file *)
tableFName = "orbitalTable_"<>ToString[N@tableViewPoint]<>".png";

(* viewpoint about the z-axis for output file *)
orbitalViewPoint = \[Pi]/4;
(* parameters for output orbital *)
nO=5;
lO=3;
mO=0;

(* Name of single orbital file *)
orbitalFName =
 "atomicOrb_"<>
 "n"<>ToString[nO]<>
 "_l"<>ToString[lO]<>
 "_m"<>ToString[mO]<>
 "_v"<>ToString[N@orbitalViewPoint]<>".png";
borRad = 1; (* value of the Bohr radius *)
isize = 500; (* Image size *)
R[n_, l_, r_, ao_] :=
 R[n, l, r, ao] =
    (* A function giving the radial wave function for the
    Hydrogen atom.
    inputs:
    n:    (int) Primary atomic number
    l:    (int) Azimuthal quantum number
    r:    (float) radius
    ao:   (float) Bohr Radius *)
  Sqrt[(2/(n*ao))^3*(n - l - 1)!/(2*n ((n + l)!))]*
   Exp[-r/n*ao]*((2*r)/(n*ao))^l*
   LaguerreL[n - l - 1, 2 l + 1, (2*r)/(n*ao)]
   
   
plotWaveFunction[n_, l_, m_, rot_] := (
    (* Here we write the orbitals in {0,+,-} basis
     where pz orbital is the same as the p0 orbital
     and the p+ and p- orbitals are formed by taking
     linear combinations of the m=+-1 states.
    inputs:
    n:    (int) Primary atomic number
    l:    (int) Azimuthal quantum number
    m:    (int) Magnetic quantum number
    rot:  (float) an angle controlling the viewpoint of
           the camera about the z-axis rot=0 returns a view
           of the x-y plane *)

  (* this controles the scale of the area being to be plotted
  currently region function is set to be a sphere about the origin
  with a radius .9*zScale *)
  zScale = 7*n;
  xScale = zScale;
  yScale = zScale;
  Legended[DensityPlot3D[
    Chop[Conjugate[
       hydrogenWaveFunction[n, l, m, Sqrt[x^2 + y^2 + z^2], 
        ArcTan[z, Sqrt[x^2 + y^2]], ArcTan[x, y]]]*
      hydrogenWaveFunction[n, l, m, Sqrt[x^2 + y^2 + z^2], 
       ArcTan[z, Sqrt[x^2 + y^2]], ArcTan[x, y]]],
    {x, -xScale*borRad, xScale*borRad},
    {y, -yScale*borRad, yScale*borRad},
    {z, -zScale*borRad, zScale*borRad},
    ColorFunction -> "SunsetColors",
    ColorFunctionScaling -> True,
    ViewVector -> {2*xScale*Cos[rot], 2*yScale*Sin[rot], 0},
    PerformanceGoal -> perf,
    PlotPoints -> pp,
    BoxRatios -> {xScale, yScale, zScale},
    PlotTheme -> "Scientific",
    Boxed -> False,
    Axes -> False,
    ImageSize -> isize,
    ImagePadding -> None,
    RegionFunction -> 
     Function[{x, y, z, f}, x^2 + y^2 + z^2 <= (.9*zScale)^2]],     
   If[MemberQ[Table[in,{in,0,n-1}],l],
   Placed[
    Framed[
     Style[
      "n=" <> ToString[n] <> "  l=" <> ToString[l] <> "  m=" <> 
       ToString[m],
      Black,
      25,
      DigitBlock -> 3],
     FrameStyle -> Black,
     RoundingRadius -> 10,
     FrameMargins -> 2,
     Background -> CMYKColor[0, 0, 0, 0, .7]],
    {.5, .125}], None]])

  
rowPlot[n_, rot_] :=
  Flatten[
  Table[
  If[And[
  MemberQ[Table[iin,{iin,0,n-1}],il],
  MemberQ[Table[iil,{iil,0,il}],Abs[im]]],
  plotWaveFunction[n, il,im, rot],
  Nothing],
   {il, 0, 3},
   {im,-3, 3}]]
  
graphGrid[rot_]:=GraphicsGrid[
	Table[
		rowPlot[in, rot], {in, 1, 6}],
	ImageSize -> {16*isize/2.25, 5*isize/1.8},
	Alignment -> Left]
	  
parGraphGrid[rot_]:=GraphicsGrid[
	ParallelTable[
		rowPlot[in, rot], {in, 1, 6}],
	ImageSize -> {16*isize/2.25, 5*isize/1.8},
	Alignment -> Left]
	
createAnimatedTable[inc_,name_]:= Export[
		FileNameJoin[{NotebookDirectory[],name}],
		ParallelTable[graphGrid[t], {t, 0, 2*\[Pi], inc}],
		AnimationRepetitions->\[Infinity]]
	
If[saveSingleOrbital, Export[Evaluate[FileNameJoin[{NotebookDirectory[], orbitalFName}]], plotWaveFunction[nO, lO, mO, orbitalViewPoint]]];
If[saveSingleTable, Export[FileNameJoin[{NotebookDirectory[], tableFName}], parGraphGrid[tableViewPoint]]];
If[saveTableAnimation, createAnimatedTable[animationIncrements, tableAnimationFName]]; 

