(* ::Package:: *)

BeginPackage["hydrogenAtom`"];
 
radialWaveFunction::usage
  "radialWaveFunction[n, l, r]
      A function giving the radial wave function for the
      Hydrogen atom.
      inputs:
      n:    (int) Primary atomic number
      l:    (int) Azimuthal quantum number
      r:    (float) radius";
      
angularWaveFunction::usage
  "angularWaveFunction[l, m, \[Theta], \[Phi]]
      A function giving the radial wave function for the
      Hydrogen atom
     inputs :
      l :    (int) Azimuthal quantum number
      m :    (int) Magnetic quantum number
      \[Theta] :    (float) polar angle
      \[Phi] :    (float) azimuthal angle";
      
hydrogenWaveFunction::usage
 "hydrogenWaveFunction[n, l, m, r, \[Theta], \[Phi]]
    A function combining the radial and
    angular parts of the wavefunction
    inputs:
    n:    (int) Primary atomic number
    l:    (int) Azimuthal quantum number
    m:    (int) Magnetic quantum number
    r:    (float) radial distance
    \[Theta]:    (float) polar angle
    \[Phi]:    (float) azimuthal angle";
    
hydrogenWaveFunctionXYZ::usage
 "hydrogenWaveFunctionXYZ[n,l,m,x,y,z]
    A function combining the radial and
    angular parts of the wavefunction
    inputs:
    n:    (int) Primary atomic number
    l:    (int) Azimuthal quantum number
    m:    (int) Magnetic quantum number
    x:    (float) x distance
    y:    (float) y distance
    z:    (float) z distance";
    
hydrogenEnergy::usage
 "hydrogenEnergy[n]
    energy levels of hydrogen
    inputs:
    n:    (int) Primary atomic number";    
    
hydrogenOscillators::usage
 "hydrogenOscillators[n]
    energy differences and oscillator strengths
	for hydrogen (relativistic correction)
    inputs:
    n:    (int) Primary atomic number";
    
hydrogenOscillatorsNR::usage
 "hydrogenOscillatorsNR[n]
    energy differences and oscillator strengths
	for hydrogen (nonrelativistic)
    inputs:
    n:    (int) Primary atomic number";
    
hydrogenFineStructure::usage
 "energyFineStructure[n,m]
    energy levels of hydrogen including 
    the fine structure
    inputs:
    n:    (int) Primary atomic number
    m:    (int) Magnetic quantum number";

energyDifference::usage
 "energyDifference[n,m]
    differences in energy levels including
	the fine structure
    inputs:
    n:    (int) Primary atomic number
    m:    (int) Magnetic quantum number";

energyDifferenceNR::usage
 "energyDifferenceNR[n]
    differences in energy levels including
    inputs:
    n:    (int) Primary atomic number";
    


Begin["`Private`"]

wper=32;
ehRy=SetPrecision[-13.5984859254535,wper];
catom=SetPrecision[137.036,wper];
(* Hartree energy in [eV] *)
EheV=SetPrecision[27.211386245988,wper];

radialWaveFunction[nNumb_, lNumb_, radius_,
	  OptionsPattern[{bohrRadius->1}]] :=
      Module[{     
      n = nNumb,
      l = lNumb,
      r = radius,
      ao = OptionValue[bohrRadius]},
      
      Sqrt[(2/(n*ao))^3*(n - l - 1)!/(2*n ((n + l)!))]
      *Exp[-r/n*ao]*((2*r)/(n*ao))^l
      *LaguerreL[n - l - 1, 2 l + 1, (2*r)/(n*ao)]]

angularWaveFunction[lNumb_, mNumb_, theta_, phi_] :=
      Module[{     
      l = lNumb,
      m = mNumb,
      \[Theta] = theta,
      \[Phi] = phi},      
      SphericalHarmonicY[l, m, \[Theta], \[Phi]]] 

hydrogenWaveFunction[nNumb_, lNumb_, mNumb_, radius_, theta_, phi_,
      OptionsPattern[{
      bohrRadius->1,
	  basis->"nlm"}]] :=
      Module[{ 
      wfp,wfm,angp,angm,rad,    
      n = nNumb,
      l = lNumb,
      m = mNumb,
      r = radius,
      \[Theta] = theta,
      \[Phi] = phi,
      ao = OptionValue[bohrRadius],
      bs = OptionValue[basis]},
      
      bs=ToString[bs];
      angp = SphericalHarmonicY[l, m, \[Theta], \[Phi]];
      rad = radialWaveFunction[n, l, r, bohrRadius -> ao];
      wfp=angp*rad; 
      
      If[bs!="xyz",
        wfp,     
        angm = SphericalHarmonicY[l, -m, \[Theta], \[Phi]];
        wfm=angm*rad;
        If[m == 0,
          wfp,
          If[m > 0,
            I/Sqrt[2] * (wfp + wfm),
            1/Sqrt[2] * (wfm - wfp)]]]]
            
hydrogenWaveFunctionXYZ[nNumb_, lNumb_, mNumb_, xVal_, yVal_, zVal_,
      OptionsPattern[{
      bohrRadius->1,
	  basis->"nlm"}]] :=
      Module[{ 
      r, \[Theta], \[Phi],    
      n = nNumb,
      l = lNumb,
      m = mNumb,
      x = xVal,
      y = yVal,
      z = zVal,
      ao = OptionValue[bohrRadius],
      bs = OptionValue[basis]},
           
      bs= ToString[bs];
      r = Sqrt[x^2+y^2+z^2];
      \[Theta] = ArcTan[z, Sqrt[x^2+y^2]];
      \[Phi] = ArcTan[x,y];          
      hydrogenWaveFunction[n,l,m, r,\[Theta],\[Phi],bohrRadius->ao,basis->bs]]
      
hydrogenEnergy[nNumb_] :=
      Module[{     
      n = nNumb},      
      ehRy*(1/n^2)/EheV]
      
hydrogenFineStructure[nNumb_,mNumb_]:=
      Module[{     
      n = nNumb,
      m = mNumb},
      hydrogenEnergy[n]*(1+1/(catom^2*n^2)*(n/(If[m==0,3/2,1/2]+1/2)-3/4))]

energyDifference[nNumb_,mNumb_]:=
      Module[{     
      n = nNumb,
      m = mNumb},
      (hydrogenFineStructure[n,m]-hydrogenEnergy[1])]

energyDifferenceNR[nNumb_]:=
      Module[{     
      n = nNumb},
      (ehRy*(1/n^2-1)/EheV)]
      
hydrogenOscillators[
  nNumb_,
  OptionsPattern[{
    bohrRadius -> 1,
    workingPrecision->MachinePrecision,
    precisionGoal-> Automatic}]] :=
 Module[{
   fosc, assums,components,
   n = nNumb,
   ao = OptionValue[bohrRadius],
   wp = OptionValue[workingPrecision],
   pg = OptionValue[precisionGoal]},
  
  assums = {ao > 0, x \[Element] Reals, y \[Element] Reals, z \[Element] Reals};
  components={x,y,z};
  
  fosc = DeleteDuplicatesBy[SortBy[
     Table[{energyDifference[n, im],
       2/3*Sum[Norm[
          Assuming[assums,
           NIntegrate[
             FullSimplify[              
              Conjugate[
                hydrogenWaveFunctionXYZ[1, 0, 0, x, y, z, 
                 bohrRadius -> ao]]*components[[ii]]*
               hydrogenWaveFunctionXYZ[n, 1, im, x, y, z, 
                bohrRadius -> ao], assums],
             {x, -\[Infinity], \[Infinity]},
             {y, -\[Infinity], \[Infinity]},
             {z, -\[Infinity], \[Infinity]},
             WorkingPrecision->wp,
             PrecisionGoal->pg] // Quiet]^2],
         {ii, Length[components]}]*energyDifference[n, im]*
        If[im == 1, 1, 2]},
      {im, -1, 1}], First], First]]
      
hydrogenOscillatorsNR[
  nNumb_,
  OptionsPattern[{
    bohrRadius -> 1,
    workingPrecision->MachinePrecision,
    precisionGoal-> Automatic}]] :=
 Module[{
   fosc, assums,components,
   n = nNumb,
   ao = OptionValue[bohrRadius],
   wp = OptionValue[workingPrecision],
   pg = OptionValue[precisionGoal]},
  
  assums = {ao > 0, x \[Element] Reals, y \[Element] Reals, z \[Element] Reals};
  components={x,y,z};
  
  {energyDifferenceNR[n],
       2/3*Sum[Sum[Norm[
Assuming[assums,
NIntegrate[
Conjugate[hydrogenWaveFunctionXYZ[1,0,0,x,y,z,bohrRadius -> ao]]*components[[ii]]*hydrogenWaveFunctionXYZ[n,1,im,x,y,z,
bohrRadius -> ao],
{x,-\[Infinity],\[Infinity]},
{y,-\[Infinity],\[Infinity]},
{z,-\[Infinity],\[Infinity]},
WorkingPrecision->wp,
PrecisionGoal->pg]//Quiet]^2],
{ii,Length[components]}]*energyDifference[n,im],{im,-1,1}]}]
     


End[]
EndPackage[]
