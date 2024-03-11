#ifndef Initialization_h
#define Initialization_h

#include "Parameters.h"
#include <getopt.h>

void Initialization(int argc, char * argv[]); // Perform initialization operations

void ParseInput(int argc, char * argv[]); // Parse command line inputs for variables

void readInput(); // Read constant inputs from file, eg Input.txt

void initBox(); // Define the box dimensions given number of beads and concentration

void initVerlet(); // Define the verlet neighbor list cutoffs

void initCell(); // Initial linked list cell variables

void initEwald(); // Calculates the constant Ewald sum k-space term M2
// See Saadat, A.; Khomami, B. Phys. Rev. E 2015, 92 (3). I use very similar notation here.

void initGrid();

void allocate(); // Allocate memory for the position, force, mobility matrix etc arrays.
// Also other miscellaneous constant definitions here.

void printOutput(); // Prints simulation parameters to the output file

void initChains(); // Initialize chain positions or read from xyz file if restarting

pid_t getpid(void); // Just here for the initRan function, not sure exactly what it does

long initRan(); // Initialize the RNG.

float ran1(long *idum); // Uniform RNG [0,1) from numerical recipes textbook

#endif // Initialization_h
