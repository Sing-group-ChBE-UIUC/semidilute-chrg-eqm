#ifndef main_h
#define main_h

#include "Parameters.h"
#include "Initialization.c"

void resetForce(); // Resets forces on beads after each time step

float gasdev(long *idum); // Box-Mueller Transform output from ran1 (uniform) to a Gaussian Dist

void bondforce(); // Bonded forces. Normally stretching and bending given by kappas, kappab.

void LJforce(); // LJ EV forces given by epsilon

void verletlist(); // Function for creating neighbor lists. See Frenkel & Smit p. 545 for details

void celllist();

void LJcell();

int binmod(int i,int j);

void ewaldrpy(); // Ewald sum of the RPY tensor to account for effects of periodic lattice images
// See Saadat, A.; Khomami, B. Phys. Rev. E 2015, 92 (3). I use very similar notation here.

void ewaldself(int ii);

void gwdecompmod(); // Geyer-Winter's Truncated Expansion Ensatz (TEA) for the decomposition of the MM matrix.
// For details on this see Geyer & Winter, JCP 130, 114905 (2009). I use very similar notation.

void printTEA();

void getNoise(); // Get Brownian random velocity from the Gaussian dist function

void updateChains(); // Apply forces, HI, and noise using the TEA decomposition.
// Note this program does not implement deformed BC's and so cannot do flowing systems, only equilibrium systems.

void updateBins();

void updateFD();

void checkVerlet(); // Check if neighbor list needs to be updated. Again see See Frenkel & Smit p. 545 for details

void printTrajectory(); // Save trajectory coordinates to the xyz file.

void printMM(); //Saves mobility matrix

void resetAverage();

void calcRgCM(); //Calculates RG + C.O.M of chains

// void printRgCM();

#endif // main_h
