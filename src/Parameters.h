#ifndef Parameters_h
#define Parameters_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
// #include <complex.h>
// #include <stdbool.h>
// #include <stdint.h>
// #include <assert.h>
#include <string.h>
// #include <limits.h>
// #include <unistd.h>
// #include <ctype.h>
// #include <sys/time.h>
// #include <sys/types.h>
#include <omp.h>

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define NDIV (1+IMM1/NTAB)
#define max(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b ? _a : _b; })
#define min(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a < _b ? _a : _b; })

/////////////////////////////
// Global Variables
/////////////////////////////
///Charge Variables///
void coulombForce(); //applies coulomb interaction
void chargeHop() ; //charges hop with coulomb interactions
void chargeHop_nonint(); //charges hop without coloumb interactions
void calcMSD(); //calculates MSD
void calcDis(); //calculates displacement of charges
void calcCH(); //adds up and prints displacement of charges
void initPos();//gets the initial position of charges
int *Charge;  //0 if uncharged //1 if charged
int *Charge_track ;  //tracks charge
int *Charge_old ; //tracks charge
int *Charge_indices ; //tracks charge
int *check_ch ; //checks if charge hopped thru PBC in dx,dy, and dz | 0 - regular displacement 1 - PBC displacement
int Ncharges ; //Number of charged beads
int i,cnumber,counter1,counter2,counter3, counter4, counter5, counter6;
double lambda_d ; //debye screening length
double lambda_b ; //bjerrum length
double kappa_d ; //inverse debye screening length(1/lambda_d)
int chargehop_type ; //0 - no charge hop, 1 - interacting charge hop, 2 - non-interacting charge hop
int MStep ; // do MC every MSteps
int *i_values ;
double barrier ; //barrier for adjacent hopping
double barrier2 ; //barrier for non-adjacent hopping
char* str ;
char* str2 ;
FILE *outputfile1, *chfile;
double Lbar[9]; // cell vectors {L1,L2,L3}
char *outp1;
FILE *datafile,*outputfile2;
int printprops;
double MSD_ch;
double MSD_com;
double conc; // concentration normalized by c_star
double *dxt, *dyt, *dzt;
double *rbi; //initial position of charge at timestep
int MSDstart; //startpoint of MSD calculation in tmax/printprops
int tcount;
double E_initial; //Initial Coulomb Energy used for kMC
//////////////////////
double c_norm; // concentration normalized by c_star
int Nb; // number of beads per chain
int Nc; // number of chains
int N; // total number of beads
int Ncm; // Number of centers of mass per chain
unsigned long tmax,t,ttemp,tstart; // Max number of time steps
double epsilon; // LJ interaction strength
double kappas; // Spring constant
double kappab; // Bending potential constant
double Rg; // Root mean squared radius of gyration
double dt; // Time step
int restart; // Restart condition. 0 -> new run (overwrite), 1 -> continue iteration
double c_star; // Overlap concentration
double box_volume; // Volume of the box given N and c_norm
double box_length; // Length of the box from the volume
double box_side; // Location of the side of the box
double sigma; // LJ hard-core overlap distance = 2a
double rc; // LJ cutoff radius
double rv,rv2; // Verlet list cutoff radius
double r2cut;
double rnew; // Distance which determines when to update neighbor list
int rn; // Number of cells in each direction
double rl; // Size of a cell
int *nlist; // Number of neighbors in neighbor list for each bead
int *list; // Neighbor list in format (bead,neighbor indices)
double *xo,*yo,*zo;
// long int ***hoc;
int *hoc;
int *ll;
double kcoeff; // 2*PI/box_length
double rtpi; // Square root PI
// alpha controls the weighting of real space and k-space sums. Lower alpha -> higher real space weight. Higher alpha -> higher k-space weight.
double alpha,alpha2,alpha3,alpha4,alpha5,alpha7; // alphan = alpha to power n
double selfconst; // Self mobility correction to Ewald sum
int kmax; // Number of wavevectors in each direction, typically converges with 3-4
int numk; // Total number of k-vectors given kmax
double *M2; // Constant coefficient of the k-space sum
// int kx,ky,kz; // Wavevector indices in x,y,z directions
// double rkx,rky,rkz; // Wavevectors in x,y,z directions
// double rkk; // Wavevector magnitude
// double m2; // Coefficient for M2
int sampf; // frequency of HI sampling, eg once every sampf time steps
double *rx,*ry,*rz; // Position arrays, eg rx[0] = x-pos of bead 0, ry[4] = y-pos of bead 4
double *px,*py,*pz;
double *rgx,*rgy,*rgz,*rg;
double *comx,*comy,*comz;
double *reex,*reey,*reez;
double *fx,*fy,*fz; // Force arrays
double *R; // Gaussian dist random velocities, R[0:2] = x,y,z random velocities of bead 0
double *C; // TEA decomposition parameters, see Geyer, T.; Winter, U. J. Chem. Phys 2009, 130.
double *Cavg;
double *Crun;
double *Dijavg;
double *Dijselfavg; // Self mobility matrix as a 2d array, eg Dij[0][0] = u_00_xx, Dij[0][2] = u_01_xx, Dij[1][3] = u_01_yy
double *Dij;
double *Dijself;
// float *Dij;
// float *Dijself;
double **D;
double *Mij;
// float *Mij;
long int *start_ij;
double p,pc; // Width of the Gaussian dist random velocity to satisfy fluctuation-dissipation
long *idum; // Random number generator seed
char *xyz; // xyz filename
char *outp; // output filename
char *matrix;
char *decomposition;
char *rgt;
char *cmt;
char *reet;
char *betavt;
char *clustero;
FILE *outputfile, *xyzfile, *MMfile, *DCfile, *Bfile, *rgtfile, *cmtfile, *reetfile,*clusterout;
double betaii; // TEA decomposition parameter
double betaij,bavg,brun; // TEA decomposition parameter
int gwcount;
// double *COMx,*COMy,*COMz; // Chain center of mass arrays
double G1,G2; // Size of grids in units of bead radii eg G1 = 1 -> first grid pts separated by 1
double g1l,g2l; // Real space position of the edge of the first and second grids
int g1i,g2i; // Bin of the edge of the first and second grids
int num_g; // number of grid points in each direction
int num_bin; // number of bins in each direction from origin
int tot_bins; // Total number of grid pts
int bin_dir; // Bins in each direction from 0,0,0
double bin_size;
double gs1,gs2,gl1;
int gd1,gd2,gn1,gn2,gnt,gh1,gh2;
long int *count,selfcount;
double eps,eps2;
double Dtr;
double *rowsum;
int itercount,itermax,iterstart;
int num_threads;
int trace;
int printperiod;
int spp;
unsigned long *tprint;
double ljtimetest,ljtime;

#endif // Parameters.h
