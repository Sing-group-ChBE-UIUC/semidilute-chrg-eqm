#include "Initialization.h"

void Initialization(int argc, char * argv[]){
    // clustero = malloc(sizeof(char)*100);
    // sprintf(clustero,"qsub/out.txt");
    ParseInput(argc,argv);
    readInput();
    initBox();
    // initVerlet();
    initCell();
    initGrid();
    allocate();
    initEwald();
    printOutput();
    initChains();
    for(i=0;i<N;i++){
      Charge[i] = 0 ;
    }
    for(i=0;i<Ncharges;i++){
      cnumber=ran1(idum)*N;
      if(Charge[cnumber]==0){
        Charge[cnumber]=1;
        Charge_indices[i]=cnumber;
        Charge_track[i]=cnumber;
        Charge_old[i]=cnumber;
      }
      else
      {
        i=i-1;
      }
    }

    // for(i=0;i<N;i++){
    //   printf("Charge: %d\n",Charge[i]);
    // }

}

void ParseInput(int argc, char * argv[]){

    // Default values if command line argument not given
    double m = 0.56701369, b = 0.184561671; // Normal LJ, attractive Scaling parameters for log(Rg) = m*log(N) + b
    // double m = 0.58942199, b = 0.155783717; // Repulsive only semiflexible eps = 1.0, Scaling parameters for log(Rg) = m*log(N) + b
    // double m = 0.606974309, b = 0.135486626; // Repulsive only flexible? Scaling parameters for log(Rg) = m*log(N) + b
    c_norm = 1.0;
    Nb = 100;
    Nc = 13;
    tmax = 10000000;
    num_threads = 1;
    trace = 1;
    Ncharges = 1;
    lambda_d = 5.0;
    restart = 0;
    barrier = 3.0;
    barrier2 = 3.0;

    int option = 0;

    if(argc < 8) printf("Warning: Using defaults for unspecified command line arguments.\n");
    while((option = getopt(argc, argv, "c:b:a:t:p:n:r:h:l:x:y:")) != -1){
        switch(option){
            case 'c':
                sscanf(optarg, "%lf", &c_norm);
                break;
            case 'b':
                sscanf(optarg, "%d", &Nb);
                break;
            case 'a':
                sscanf(optarg, "%d", &Nc);
                break;
            case 't':
                sscanf(optarg, "%lu", &tmax);
                break;
            case 'p':
                sscanf(optarg, "%d", &trace);
                break;
            case 'n':
                sscanf(optarg, "%d", &num_threads);
                break;
            case 'r':
                sscanf(optarg, "%d", &restart);
                break;

            //Charge Stuff//  //also add variable to option list in line 60
            case 'h':
                sscanf(optarg, "%d", &Ncharges);
                break;
            case 'l':
                sscanf(optarg, "%lf", &lambda_d);
                break;
            case 'x':
                sscanf(optarg, "%lf", &barrier);
                break;
            case 'y':
                sscanf(optarg, "%lf", &barrier2);
                break;
            ////////////////
            case '?':
                printf("Unknown option -%c.\n Execution abort.", optopt);
                exit(EXIT_FAILURE);
        }
    }
    N = Nb*Nc; // Total number of beads

    if(Ncharges>N){
      printf("Error - Ncharges exceeds N, exiting\n");
      exit(1);
    }

    //double b = 0.7924, m = 0.6392; // Kremer-Grest, Rg = b*N_{b}^m to match old PEF
    //double b = 0.4082, m = 3/5; //scaling coefficient in good solvent
    //Rg = b*pow(Nb,m);
    Rg = pow(10,m*log10(Nb)+b);
    // printf("%lf\n",Rg);

}

void readInput(){

    FILE *inputfile;
    inputfile = fopen("Input.txt","r");
    fscanf(inputfile, "epsilon = %lf\n", &epsilon);
    fscanf(inputfile, "kappas = %lf\n", &kappas);
    fscanf(inputfile, "kappab = %lf\n", &kappab);
    // fscanf(inputfile, "Rg = %lf\n", &Rg); // Mean Squared Radius of Gyration
    fscanf(inputfile, "dt = %lf\n", &dt);
    fscanf(inputfile, "itermax = %d\n", &itermax);
    fscanf(inputfile, "iterstart = %d\n", &iterstart);
    //Charge Stuff//
    fscanf(inputfile, "lambda_b = %lf\n", &lambda_b);
    fscanf(inputfile, "chargehop_type = %d\n", &chargehop_type);
    fscanf(inputfile, "MStep = %d\n", &MStep);
    fscanf(inputfile, "MSDstart = %d\n", &MSDstart);
    /////////////////
    fclose(inputfile);
}

void initBox(){

    int itertemp,tempthreads;
    unsigned long ttemp;
    // Set N, find box length for c = c_star*c_norm
    c_star = Nb/(M_PI*4/3*Rg*Rg*Rg); // Overlap concentration
    conc = c_star*c_norm;
    // c_star = Nb/(4/3*M_PI*Rg*Rg*Rg); // *** THIS IS WRONG!!! *** For some reason, 4/3* ignored in this form, use M_PI*4/3
    box_volume = N/(c_norm*c_star); // Box volume at the normalized concentration
    box_length = pow(box_volume,(double)1/3); // For cubic lattice, L = V**(1/3)
    box_side = box_length/2; // Box centered at origin, so sides at +- L/2

    if(box_length<3*Nc*Rg){
      printf("Warning - Box Too Small (Make larger than %lf)\n",3*Nc*Rg);
      //exit(1);
    }

    outp = malloc(sizeof(char)*100);
    // sprintf(outp, "txt/P%.3lf_%d_%d_%d.txt", c_norm, Nb, Nc, trace);
    sprintf(outp, "txt/P%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%d_%d.txt",Nb,Nc,c_norm,box_length,Ncharges,barrier,barrier2,lambda_d,trace,iterstart);
    //printf("%s\n",outp);

    if(restart==1){
        outputfile = fopen(outp, "r");
        if(!outputfile){
            printf("Error: trying to restart %s with no existing files\n",outp);
            printf("Try checking input or change to restart = 0\n");
            exit(1);
        }
        fclose(outputfile);
    }

    //Charge Stuff//
    Lbar[0] = box_length;
		Lbar[1] = 0.0;
		Lbar[2] = 0.0;
		Lbar[3] = 0.0;
		Lbar[4] = box_length;
		Lbar[5] = 0.0;
		Lbar[6] = 0.0;
		Lbar[7] = 0.0;
		Lbar[8] = box_length;
    ///////////////
}

void initVerlet(){

    int i;
    sigma = 2.0; // 1 particle diameter = 2 radii, HS overlap, get LJ cutoff from this
    rc = 2.5*sigma; // LJ cutoff range (assumption: goes to zero outside this range)
    r2cut = rc*rc;
    rv =  5.0*sigma; // Verlet list radius. Particles inside this radius are added to the varlet list
    rv2 = rv*rv;
    rnew = 0.5*(rv-rc);
    // printf("%lf %lf %lf\n",rc,rv,rnew);
    nlist = calloc(N, sizeof(long int)); // Number of neighbors in neighbor list for each bead
    // List for neighbors in format list(bead,neighbor indices)
    // Size NxN because there are N beads requiring neighbor lists and a maximum of N possible neighbors for each bead
    list = calloc(N*N, sizeof(long int));
    xo = calloc(N,sizeof(double));
    yo = calloc(N,sizeof(double));
    zo = calloc(N,sizeof(double));
}

void initCell(){
    int i,j;
    sigma = 2.0;
    rc = 2.5*sigma;
    r2cut = rc*rc;
    rn = floor(box_length/rc);
    rl = box_length/rn;
    hoc = calloc(rn*rn*rn,sizeof(int));
    ll = calloc(N,sizeof(int));
    // printf("%lf %d %lf\n",rc,rn,rl);
}

void initEwald(){

    int i,j,k,l,m,n,o,kk,tempbx,tempby,tempbz,tempiter,tempcount,a,b,start,kx,ky,kz;
    double rkx,rky,rkz,rkk,m2;
    gwcount = 0;
    sampf = 10000;
    kmax = 3;
	kcoeff = 2*M_PI/box_length;
    // printf("%lf %lf\n",kcoeff,box_length);
    // exit(1);
	alpha = 6/box_length;
    alpha2 = alpha*alpha; alpha3 = alpha2*alpha; alpha4 = alpha2*alpha2; alpha5 = alpha4*alpha; alpha7 = alpha4*alpha2*alpha;
	rtpi = sqrt(M_PI);
	selfconst = 1-6/rtpi*alpha + 40/(3*rtpi)*alpha3;
    numk = kmax*2 + 1;
    M2 = calloc(numk*numk*numk*9, sizeof(double));
    for(kx = -kmax; kx < kmax+1; kx++){
        rkx = kx*kcoeff;
        for(ky = -kmax; ky < kmax+1; ky++){
            rky = ky*kcoeff;
            for(kz = -kmax; kz < kmax +1; kz++){
                rkz = kz*kcoeff;
                rkk = rkx*rkx + rky*rky + rkz*rkz;
                m2 = (1.0-rkk/3.0)*(1.0+rkk/(4.0*alpha2)+rkk*rkk/(8.0*alpha4))*6.0*M_PI/rkk*exp(-rkk/(4.0*alpha2));
                start = 9*((kx+kmax)*numk*numk+(ky+kmax)*numk+(kz+kmax));
                M2[start] = m2*(1-rkx*rkx/rkk)/box_volume;
                M2[start+1] = m2*(-rkx*rky/rkk)/box_volume;
                M2[start+2] = m2*(-rkx*rkz/rkk)/box_volume;
                M2[start+3] = m2*(-rky*rkx/rkk)/box_volume;
                M2[start+4] = m2*(1-rky*rky/rkk)/box_volume;
                M2[start+5] = m2*(-rky*rkz/rkk)/box_volume;
                M2[start+6] = m2*(-rkz*rkx/rkk)/box_volume;
                M2[start+7] = m2*(-rkz*rky/rkk)/box_volume;
                M2[start+8] = m2*(1-rkz*rkz/rkk)/box_volume;
            }
        }
    }
    Dtr = selfconst;
    for(kx = -kmax; kx < kmax+1; kx++){
        for(ky = -kmax; ky < kmax+1; ky++){
            for(kz = -kmax; kz < kmax +1; kz++){
                kk = kx*kx + ky*ky + kz*kz;
                if(kk != 0){
                    Dtr += M2[9*((kx+kmax)*numk*numk+(ky+kmax)*numk+(kz+kmax))];
                }
            }
        }
    }
    eps = 0.0; eps2 = 0.0;
    if(restart==1){
        MMfile = fopen(matrix,"r");
        fscanf(MMfile,"%lu %ld\n",&tstart,&selfcount);
        for(i=0;i<Nb;i++){
            for(j=0;j<Nb;j++){
                start = 9*Nb*i + 9*j;
                fscanf(MMfile,"%d %d\n",&tempbx,&tempby);
                for(a=0;a<3;a++){
                    for(b=0;b<3;b++){
                        fscanf(MMfile,"%le ",&Dijselfavg[start+3*a+b]);
                    }
                    fscanf(MMfile,"\n");
                }
            }
        }
        for(i=0;i<gd1;i++){
            for(j=0;j<gd1;j++){
                for(k=0;k<gd1;k++){
                    start = i*gd1*gd1+j*gd1+k;
                    fscanf(MMfile,"%d %d %d %ld\n",&tempbx,&tempby,&tempbz,&count[start]);
                    for(l=0;l<3;l++){
                        for(m=0;m<3;m++){
                            fscanf(MMfile,"%le ",&Dijavg[9*start+3*l+m]);
                        }
                        fscanf(MMfile,"\n");
                    }
                }
            }
        }
        for(i=0;i<gd2;i++){
            for(j=0;j<gd2;j++){
                for(k=0;k<gd2;k++){
                    start = gn1 + i*gd2*gd2 + j*gd2 + k;
                    fscanf(MMfile,"%d %d %d %ld\n",&tempbx,&tempby,&tempbz,&count[start]);
                    for(l=0;l<3;l++){
                        for(m=0;m<3;m++){
                            fscanf(MMfile,"%le ",&Dijavg[9*start+3*l+m]);
                        }
                        fscanf(MMfile,"\n");
                    }
                }
            }
        }
        fclose(MMfile);
        for(i=0;i<9*Nb*Nb;i++){
            Dijselfavg[i] *= selfcount;
        }
        for(i=0;i<gnt;i++){
            for(j=0;j<9;j++){
                Dijavg[9*i+j] *= count[i];
            }
        }
        DCfile = fopen(decomposition,"r");
        fscanf(DCfile,"%lu %le\n",&ttemp,&brun);
        for(i=0;i<N;i++){
            fscanf(DCfile,"%d %le %le %le\n",&tempbx,&Crun[3*i],&Crun[3*i+1],&Crun[3*i+2]);
        }
        gwcount = ttemp/sampf + 1;
        brun *= gwcount;
        for(i=0;i<3*N;i++){
            Crun[i] *= gwcount;
        }
    }
    if(iterstart>0){
        sprintf(matrix, "mm/M%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%d_%d.txt",Nb,Nc,c_norm,box_length,Ncharges,barrier,barrier2,lambda_d,trace,iterstart-1);
        MMfile = fopen(matrix,"r");
        if(!MMfile){
            printf("Error: No M file from previous iteration\n");
            exit(1);
        }
            fscanf(MMfile,"%lu %d\n",&ttemp,&tempcount);
            for(i=0;i<Nb;i++){
                for(j=0;j<Nb;j++){
                    start = 9*Nb*i + 9*j;
                    fscanf(MMfile,"%d %d\n",&tempbx,&tempby);
                    for(a=0;a<3;a++){
                        for(b=0;b<3;b++){
                            // fscanf(MMfile,"%lf ",&Dijself[start+3*a+b]);
                            fscanf(MMfile,"%le ",&Dijself[start+3*a+b]);
                        }
                        fscanf(MMfile,"\n");
                    }
                }
            }
            for(i=0;i<gd1;i++){
                for(j=0;j<gd1;j++){
                    for(k=0;k<gd1;k++){
                        start = i*gd1*gd1+j*gd1+k;
                        fscanf(MMfile,"%d %d %d %d\n",&tempbx,&tempby,&tempbz,&tempcount);
                        for(l=0;l<3;l++){
                            for(m=0;m<3;m++){
                                // fscanf(MMfile,"%lf ",&Dij[9*start+3*l+m]);
                                fscanf(MMfile,"%le ",&Dij[9*start+3*l+m]);
                            }
                            fscanf(MMfile,"\n");
                        }
                    }
                }
            }
            for(i=0;i<gd2;i++){
                for(j=0;j<gd2;j++){
                    for(k=0;k<gd2;k++){
                        start = gn1 + i*gd2*gd2 + j*gd2 + k;
                        fscanf(MMfile,"%d %d %d %d\n",&tempbx,&tempby,&tempbz,&tempcount);
                        for(l=0;l<3;l++){
                            for(m=0;m<3;m++){
                                // fscanf(MMfile,"%lf ",&Dij[9*start+3*l+m]);
                                fscanf(MMfile,"%le ",&Dij[9*start+3*l+m]);
                            }
                            fscanf(MMfile,"\n");
                        }
                    }
                }
            }
        fclose(MMfile);
        sprintf(matrix, "mm/M%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%d_%d.txt",Nb,Nc,c_norm,box_length,Ncharges,barrier,barrier2,lambda_d,trace,iterstart);
        sprintf(decomposition, "decomp/DC%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%d_%d.txt",Nb,Nc,c_norm,box_length,Ncharges,barrier,barrier2,lambda_d,trace,iterstart-1);
        DCfile = fopen(decomposition,"r");
        if(!DCfile){
            printf("Error: No DC file from previous iteration\n");
            exit(1);
        }
        while(!feof(DCfile)){
            fscanf(DCfile,"%lu %le",&ttemp,&bavg);
            for(i=0;i<N;i++){
                fscanf(DCfile,"%d %le %le %le\n",&tempbx,&Cavg[3*i],&Cavg[3*i+1],&Cavg[3*i+2]);
            }
        }
        fclose(DCfile);
        sprintf(decomposition, "decomp/DC%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%d_%d.txt",Nb,Nc,c_norm,box_length,Ncharges,barrier,barrier2,lambda_d,trace,iterstart);
    }
    if(restart==2){
        DCfile = fopen(decomposition,"w");
        fprintf(DCfile,"");
        fclose(DCfile);
        Bfile = fopen(betavt,"w");
        fprintf(Bfile,"");
        fclose(Bfile);
        // Cfile = fopen(cvt,"w");
        // fprintf(Cfile,"");
        // fclose(Cfile);
        MMfile = fopen(matrix,"w");
        fprintf(MMfile,"");
        fclose(MMfile);
    }
}

void allocate(){

    int i,j,k,l;
    printperiod = 1000;
    printprops = 1000;
    spp = printperiod*dt; // Samples per period
    xyz = malloc(sizeof(char)*100);
    // outp = malloc(sizeof(char)*100);
    matrix = malloc(sizeof(char)*100);
    betavt = malloc(sizeof(char)*100);
    // cvt = malloc(sizeof(char)*100);
    decomposition = malloc(sizeof(char)*100);
    rgt = malloc(sizeof(char)*100);
    cmt = malloc(sizeof(char)*100);
    reet = malloc(sizeof(char)*100);
    // sprintf(xyz, "xyz/R%.3lf_%d_%d_%d_%d.xyz", c_norm, Nb, Nc, trace, iterstart);
    // sprintf(decomposition, "decomp/DC%.3lf_%d_%d_%d_%d.txt", c_norm, Nb, Nc, trace, iterstart);
    // sprintf(matrix, "mm/M%.3lf_%d_%d_%d_%d.txt", c_norm, Nb, Nc, trace, iterstart);
    // sprintf(betavt, "decomp/B%.3lf_%d_%d_%d_%d.txt", c_norm, Nb, Nc, trace, iterstart);
    // sprintf(rgt, "prop/RG%.3lf_%d_%d_%d_%d.txt", c_norm, Nb, Nc, trace, iterstart);
    // sprintf(cmt, "prop/CM%.3lf_%d_%d_%d_%d.txt", c_norm, Nb, Nc, trace, iterstart);
    sprintf(xyz, "xyz/R%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%d_%d.xyz",Nb,Nc,c_norm,box_length,Ncharges,barrier,barrier2,lambda_d,trace,iterstart);
    sprintf(decomposition, "decomp/DC%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%d_%d.txt",Nb,Nc,c_norm,box_length,Ncharges,barrier,barrier2,lambda_d,trace,iterstart);
    sprintf(matrix, "mm/M%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%d_%d.txt",Nb,Nc,c_norm,box_length,Ncharges,barrier,barrier2,lambda_d,trace,iterstart);
    sprintf(betavt, "decomp/B%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%d_%d.txt",Nb,Nc,c_norm,box_length,Ncharges,barrier,barrier2,lambda_d,trace,iterstart);
    sprintf(rgt, "prop/RG%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%d_%d.txt",Nb,Nc,c_norm,box_length,Ncharges,barrier,barrier2,lambda_d,trace,iterstart);
    sprintf(cmt, "prop/CM%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%d_%d.txt",Nb,Nc,c_norm,box_length,Ncharges,barrier,barrier2,lambda_d,trace,iterstart);
    sprintf(reet, "prop/REE%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%d_%d.txt",Nb,Nc,c_norm,box_length,Ncharges,barrier,barrier2,lambda_d,trace,iterstart);

    //Charge Stuff//
    str2 = malloc(sizeof(char)*300);
    outp1 = malloc(sizeof(char)*300);
    sprintf(str2,"prop/CH%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%d_%d.txt",Nb,Nc,c_norm,box_length,Ncharges,barrier,barrier2,lambda_d,trace,iterstart);
    sprintf(outp1,"prop/MSD%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%d_%d.txt",Nb,Nc,c_norm,box_length,Ncharges,barrier,barrier2,lambda_d,trace,iterstart);
    Charge  = calloc(N, sizeof(int));
    Charge_track = calloc(Ncharges, sizeof(int));
    Charge_old = calloc(Ncharges, sizeof(int));
    Charge_indices = calloc(Ncharges, sizeof(int));
    check_ch = calloc(3*Ncharges, sizeof(int));
    i_values = calloc(N, sizeof(int));
    dxt = calloc(Ncharges*printprops, sizeof(double));
  	dyt = calloc(Ncharges*printprops, sizeof(double));
  	dzt = calloc(Ncharges*printprops, sizeof(double));
  	rbi = calloc(3*Ncharges, sizeof(double));
    ////////////////

    // sprintf(cvt, "decomp/C%.3lf_%d_%d.txt", c_norm, N,trace);
    rx = calloc(N, sizeof(double));
    ry = calloc(N, sizeof(double));
    rz = calloc(N, sizeof(double));
    px = calloc(N, sizeof(double));
    py = calloc(N, sizeof(double));
    pz = calloc(N, sizeof(double));
    // rgx = calloc(Nc,sizeof(double));
    // rgy = calloc(Nc,sizeof(double));
    // rgz = calloc(Nc,sizeof(double));
    rg = calloc(Nc,sizeof(double));
    comx = calloc(Nc,sizeof(double));
    comy = calloc(Nc,sizeof(double));
    comz = calloc(Nc,sizeof(double));
    reex = calloc(Nc,sizeof(double));
    reey = calloc(Nc,sizeof(double));
    reez = calloc(Nc,sizeof(double));
    // for(i=0;i<Nc;i++){
    //     rgx[i] = calloc(spp,sizeof(double));
    //     rgy[i] = calloc(spp,sizeof(double));
    //     rgz[i] = calloc(spp,sizeof(double));
    //     rg[i] = calloc(spp,sizeof(double));
    //     comx[i] = calloc(spp,sizeof(double));
    //     comy[i] = calloc(spp,sizeof(double));
    //     comz[i] = calloc(spp,sizeof(double));
    //     reex[i] = calloc(spp,sizeof(double));
    //     reey[i] = calloc(spp,sizeof(double));
    //     reez[i] = calloc(spp,sizeof(double));
    // }
    // tprint = calloc(spp,sizeof(unsigned long));
    fx = calloc(N, sizeof(double));
    fy = calloc(N, sizeof(double));
    fz = calloc(N, sizeof(double));
    R = calloc(N*3, sizeof(double));
    C = calloc(3*N, sizeof(double));
    Cavg = calloc(3*N,sizeof(double));
    Crun = calloc(3*N,sizeof(double));
    Dijselfavg = calloc(9*Nb*Nb, sizeof(double));
    Dijself = calloc(9*Nb*Nb, sizeof(double));
    Dijavg = calloc(9*gnt, sizeof(double));
    Dij = calloc(9*gnt, sizeof(double));
    // D = calloc(3*N,sizeof(double));
    // for(i=0;i<3*N;i++){
    //     D[i] = calloc(3*N,sizeof(double));
    // }
    // Mij = calloc(9*N*N,sizeof(double));
    // Mij = calloc(9*(N*(N-1)/2+N),sizeof(double));
    start_ij = calloc(N*N,sizeof(long int));
    p = sqrt(2.0/dt);
    pc = sqrt(2.0*dt);
    idum = malloc(sizeof(long)); // = initRan(); // for testing, set seed
    *idum = initRan();
    //*idum = 2;
    // *idum = 541207347;
    if(*idum==2){
        printf("Warning! Seed is fixed!\n");
    }
    betaii = 1.0;
    brun = 0.0;
    count = calloc(gnt, sizeof(long int));
    rowsum = calloc(3*N, sizeof(double));

}

void printOutput(){

    omp_set_num_threads(num_threads);
    // #pragma omp parallel
    // {
    //     #pragma omp single
    //         num_threads = omp_get_num_threads();
    // }

	outputfile = fopen(outp, "w");
	fprintf(outputfile, "epsilon = %lf\n", epsilon);
	fprintf(outputfile, "kappas = %lf\n", kappas);
	fprintf(outputfile, "kappab = %lf\n", kappab);
  fprintf(outputfile, "c_norm = %lf\n", c_norm);
  fprintf(outputfile, "conc = %lf\n", conc);
	fprintf(outputfile, "N = %d\n", N);
	fprintf(outputfile, "Nb = %d\n", Nb);
	fprintf(outputfile, "Nc = %d\n", Nc);
  ///charge stuff
  fprintf(outputfile,"Ncharges = %d\n",Ncharges);
  fprintf(outputfile,"lambda_d = %lf\n",lambda_d);
  fprintf(outputfile,"lambda_b = %lf\n",lambda_b);
  fprintf(outputfile,"MStep = %d\n",MStep);
  fprintf(outputfile,"chargehop_type = %d\n",chargehop_type);
  fprintf(outputfile,"barrier = %lf\n",barrier);
  fprintf(outputfile,"barrier2 = %lf\n",barrier2);
  //
	fprintf(outputfile, "Rg = %lf\n", Rg);
    fprintf(outputfile, "dt = %lf\n", dt);
    fprintf(outputfile, "tmax = %lu\n", tmax);
    fprintf(outputfile, "itermax = %d\n", itermax);
	fprintf(outputfile, "box_length = %lf\n", box_length);
	fprintf(outputfile, "rv = %lf\n", rv);
	fprintf(outputfile, "kmax = %d\n", kmax);
	fprintf(outputfile, "alpha = %lf\n", alpha);
    fprintf(outputfile, "gd1 = %d\n", gd1);
    fprintf(outputfile, "gd2 = %d\n", gd2);
    fprintf(outputfile, "gnt = %d\n", gnt);
    fprintf(outputfile, "gs1 = %lf\n", gs1);
    fprintf(outputfile, "gs2 = %lf\n", gs2);
    fprintf(outputfile, "num_threads = %d\n",num_threads);
    fprintf(outputfile, "trace = %d\n", trace);
    fprintf(outputfile, "restart = %d\n",restart);
	fprintf(outputfile, "\n");
	fprintf(outputfile, "\n");
    fprintf(outputfile, "SEED %ld\n", *idum);
	fclose(outputfile);

}

void initChains(){

    int i,j,k,index,monindex,test,tempN,initcount,inittest;
    double phi,theta,dx,dy,dz;
    char tempname[1024];
    // If starting a new run, clear the xyz file and randomly distribute non-overlapping chains in the box
    if(restart==0){
        //xyzfile = fopen(xyz, "r");
        // if(xyzfile){
        //     printf("Error: trying to write over an existing simulation\n");
        //     printf("Try checking input or change to restart = 2\n");
        //     exit(1);
        // }
        xyzfile = fopen(xyz,"w");
        fprintf(xyzfile,"");
        fclose(xyzfile);
        rgtfile = fopen(rgt,"w");
        fprintf(rgtfile,"");
        fclose(rgtfile);
        cmtfile = fopen(cmt,"w");
        fprintf(cmtfile,"");
        fclose(cmtfile);
        reetfile = fopen(reet,"w");
        fprintf(reetfile,"");
        fclose(reetfile);
        //Charge Stuff/////
        chfile = fopen(str2,"w");
    		fprintf(chfile, "");
    		fclose(chfile);
        outputfile=fopen(outp1,"w");
    		fprintf(outputfile, "");
    		fclose(outputfile);
        //outputfile2=fopen("disp1.csv","w");
    		//fprintf(outputfile2, "");
    		//fclose(outputfile2);
        //////////////////
        tstart = 0;
        inittest = 0;
        while(inittest==0){
            inittest = 1;
            for(i = 0; i<Nc; ++i){
                initcount = 0;
                test = 0;
                index = Nb*i; // Index of first monomer in  i
                while(test==0){
                    test = 1;
                    rx[index] = ran1(idum)*box_length - box_side;
                    ry[index] = ran1(idum)*box_length - box_side;
                    rz[index] = ran1(idum)*box_length - box_side;
                    for(k = 0; k<index; k++){
                        dx = rx[index]-rx[k]; dy = ry[index]-ry[k]; dz = rz[index]-rz[k];
                        dx -= box_length*round(dx/box_length); dy -= box_length*round(dy/box_length); dz -= box_length*round(dz/box_length);
                        if(dx*dx+dy*dy+dz*dz<5.0)
                        {
                            test = 0;
                        }
                    }
                    initcount++;
                    if(initcount>1e3){
                        // printf("max attempts reached\n");
                        // exit(1);
                        inittest = 0;
                        break;
                    }
                }
                if(inittest==0){
                    break;
                }
                for(j=1;j<Nb;j++){
                    monindex = index+j; // Index of monomer j in  i
                    test = 0;
                    while(test==0){
                        test = 1;
                        theta = ran1(idum)*2.0*M_PI;
                        phi = acos(2.0*ran1(idum)-1.0);
                        rx[monindex] = rx[monindex-1] + 2.05*sin(phi)*cos(theta);
                        ry[monindex] = ry[monindex-1] + 2.05*sin(phi)*sin(theta);
                        rz[monindex] = rz[monindex-1] + 2.05*cos(phi);
                        rx[monindex] -= box_length*round(rx[monindex]/box_length);
                        ry[monindex] -= box_length*round(ry[monindex]/box_length);
                        rz[monindex] -= box_length*round(rz[monindex]/box_length);
                        for(k = 0; k<monindex; k++){
                            dx = rx[monindex]-rx[k]; dy = ry[monindex]-ry[k]; dz = rz[monindex]-rz[k];
                            dx -= box_length*round(dx/box_length); dy -= box_length*round(dy/box_length); dz -= box_length*round(dz/box_length);
                            if(dx*dx+dy*dy+dz*dz<4.0)
                            {
                                test = 0;
                            }
                        }
                        initcount++;
                        if(initcount>1e3){
                            // printf("max attempts reached\n");
                            // exit(1);
                            inittest = 0;
                            break;
                        }
                    }
                    if(inittest==0){
                        break;
                    }
                }
                if(inittest==0){
                    break;
                }
            }
        }
    }
    else if(restart==1) {
        // sprintf(xyz, "xyz/readin_%.3lf.xyz",c_norm);
        xyzfile = fopen(xyz,"r");
        if(!xyzfile){
            printf("Error: trying to continue %s (does not exist)\n",xyz);
            printf("Try checking input or change to restart = 0\n");
            exit(1);
        }
        while(!feof(xyzfile)){
            fscanf(xyzfile, "%d\n", &tempN);
            fscanf(xyzfile, "%lu\n", &tstart);
            for(j=0;j<tempN;j++){
                fscanf(xyzfile, "A %lf %lf %lf\n",&rx[j],&ry[j],&rz[j]);
            }
            // if(ttemp==tstart){
            //     break;
            // }
        }
        fclose(xyzfile);
    }
}

void initGrid(){

    // double bin_sep = 2.0;
    // // if(restart==0){
    // //     num_g = (int)(box_length/bin_sep);
    // //     bin_size = box_length/num_g;
    // //     if(num_g%2==0){
    // //         num_g += 1;
    // //     }
    // // }
    // num_g = (int)(box_length/2.0/bin_sep);
    // num_g = 2*num_g+1;
    // bin_size = box_length/num_g;
    // tot_bins = num_g*num_g*num_g;
    // double gs1,gs2,gl1;
    // int gd1,gd2,gn1,gn2,gnt;
    gs2 = 4.0;
    gd2 = 2*ceil(box_length/2.0/gs2) + 1;
    gs2 = box_length/2.0/((gd2-1)/2);
    gs1 = gs2/2.0;
    gl1 = 30.0;
    gd1 = 2*ceil(gl1/gs1) + 1;
    gl1 = gs1*(gd1-1)/2;
    gn1 = gd1*gd1*gd1;
    gn2 = gd2*gd2*gd2;
    gnt = gn1+gn2;
    gh1 = (gd1-1)/2;
    gh2 = (gd2-1)/2;
    // int gridtest = round(-1.5/1.0);
    // printf("%d\n",gridtest);
    // exit(1);
    // printf("gs1 %lf gs2 %lf gd1 %d gd2 %d gn1 %d gn2 %d gh1 %d gh2 %d gl1 %lf gnt %d\n",gs1,gs2,gd1,gd2,gn1,gn2,gh1,gh2,gl1,gnt);
    // exit(1);
}

float ran1(long *idum){
	int j;
	long k;
	static long idum2 = 123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if(*idum <= 0)
	{
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2 = (*idum);
		for(j=NTAB+7;j>=0;--j)
		{
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if(*idum<0) *idum+=IM1;
			if(j<NTAB) iv[j] = *idum;
		}
		iy = iv[0];
	}
	k = (*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if(*idum<0) *idum += IM1;
	k=idum2/IQ2;
	if(*idum<0) idum2+= IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if(iy<1) iy += IMM1;
	if((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

long initRan(){
    //time_t seconds;
    //time(&seconds);
    //return -1*(unsigned long)(seconds/12345); This is bad.  :(

    //This will hopefully allow us to have a unique seed even if executed multiple times a second-Got from Mike
    //http://stackoverflow.com/questions/322938/recommended-way-to-initialize-srand
    unsigned long a = clock();
    unsigned long b = time(NULL);
    unsigned long c = getpid();
    a=a-b;  a=a-c;  a=a^(c >> 13);
    b=b-c;  b=b-a;  b=b^(a << 8);
    c=c-a;  c=c-b;  c=c^(b >> 13);
    a=a-b;  a=a-c;  a=a^(c >> 12);
    b=b-c;  b=b-a;  b=b^(a << 16);
    c=c-a;  c=c-b;  c=c^(b >> 5);
    a=a-b;  a=a-c;  a=a^(c >> 3);
    b=b-c;  b=b-a;  b=b^(a << 10);
    c=c-a;  c=c-b;  c=c^(b >> 15);
    return c%1000000000; //careful here.  Another 0 might break the ran1 (long long instead of just long)
}
