# Program to create density plots of hydrogen atomic orbitals
* Program has the option of creating a single atomic orbial, a table of all atomic orbitals (up to n=6, l=3, m=+-3 by default), and an animated gif of atomic orbitals (up to n=6, l=3, m=+-3 by default) rotating around the z axis.

## Instructions
### Running code will create plots to be saved in the same directory as the directory which this script is saved
* saveTableAnimation, saveSingleTable, and saveSingleOrbital are bools allowing user to select which type of file to save
* fastPlot is a bool allowing user to run low resolution code for testing or higher resolution for final plots.
* Currently code is very memory intensive adjusting pp variable will lower resolution as well as lowring run time and memory used.

## Example outputs
* A plot of a single orbital n=5, l=3, m=0, rotated pi/4 about the z axis run at pp=100 (changing the filename within the code from foo.png to foo.pdf will increase the resolution as well as the file size)
![](atomicOrb_n5_l3_m0_v0.785398.png)

* A table of orbitals up to n=6, l=3, m=+-3 taken at pp=100  (changing the filename within the code from foo.png to foo.pdf will increase the resolution as well as the file size)
![](orbitalTable_0.785398.png)

* An animated gif for a table orbitals up to n=6, l=3, m=+-3 taken at pp=50 rotating the camera view abour the z axis.
![](orbitalTableAnimation.gif)

## Potential Improvements
* Adjusting the radius under region function so that it doesn't cut of density plots for high n and l.

* There are several places to improve current code involving memory usage and efficiency
- By adding more read/writes animation can be made to be higher resolution while using less memory. Writing files to hard drive and later uploading when creating gif is probably more efficient.

* Adjusting the opacity may make plots more clear.

