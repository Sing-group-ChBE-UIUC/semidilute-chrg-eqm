#include "main.h"

/////////////////////////
// Structure of the Code
////////////////////////
// Calculation of forces and updates are all performed using functions in the main.c file
// Some setup operations which only need to be performed at the start of the simulation are
// seperated into another file, Initialization.c.
// Variables and constants not defined within functions can be found in the Parameters.h file
// along with a brief description and the indexing format.
// For declaration and descriptions of various functions, see the accompanying .h file

int main (int argc,char * argv[]){

	clock_t tic,toc;
	double start_time,time_step,ksteps,dx,dy,dz,r;
	int tdelay = 0,i;
	ljtimetest = 0.0;

	Initialization(argc,argv); // Some initialization steps. See Initialization.c
	// verletlist(); // Create the initial neighbor lists
	for(itercount=iterstart;itercount<itermax;itercount++){
		sprintf(xyz, "xyz/R%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%d_%d.xyz",Nb,Nc,c_norm,box_length,Ncharges,barrier,barrier2,lambda_d,trace,itercount);
	    sprintf(decomposition, "decomp/DC%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%d_%d.txt",Nb,Nc,c_norm,box_length,Ncharges,barrier,barrier2,lambda_d,trace,itercount);
	    sprintf(matrix, "mm/M%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%d_%d.txt",Nb,Nc,c_norm,box_length,Ncharges,barrier,barrier2,lambda_d,trace,itercount);
	    sprintf(betavt, "decomp/B%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%d_%d.txt",Nb,Nc,c_norm,box_length,Ncharges,barrier,barrier2,lambda_d,trace,itercount);
	    sprintf(rgt, "prop/RG%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%d_%d.txt",Nb,Nc,c_norm,box_length,Ncharges,barrier,barrier2,lambda_d,trace,itercount);
	    sprintf(cmt, "prop/CM%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%d_%d.txt",Nb,Nc,c_norm,box_length,Ncharges,barrier,barrier2,lambda_d,trace,itercount);
		sprintf(reet, "prop/REE%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%d_%d.txt",Nb,Nc,c_norm,box_length,Ncharges,barrier,barrier2,lambda_d,trace,itercount);
		sprintf(str2,"prop/CH%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%d_%d.txt",Nb,Nc,c_norm,box_length,Ncharges,barrier,barrier2,lambda_d,trace,itercount);
	  sprintf(outp1,"prop/MSD%d_%d_%.6lf_%.6lf_%d_%.2lf_%.2lf_%.2lf_%d_%d.txt",Nb,Nc,c_norm,box_length,Ncharges,barrier,barrier2,lambda_d,trace,itercount);
		if(restart==0){
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
			//////////////////
		}


		start_time = omp_get_wtime();
		ksteps = omp_get_wtime();
		tcount=tstart%printprops;
	    for(t = tstart; t<(tmax+1); ++t){
			initPos();
			resetForce(); // Reset forces after each time step
			celllist(); // Renew cell lists
			bondforce(); // Calculate bonded forces from  connectivity (eg stretching, bending)
			// LJforce(); // Calculate non-bonded forces using neighbor lists (eg excluded volume aka LJ)
			LJcell();
			if(chargehop_type==0){
				coulombForce();
			}
			else if(chargehop_type==1){
				coulombForce();
				chargeHop();
			}
			else if(chargehop_type==2){
				coulombForce();
				chargeHop_nonint();
			}
			else{
				printf("Error: chargehop_type %d not recognized - exiting\n",chargehop_type);
				exit(1);
			}
			// if(t%sampf==0 && t%100000!=0){
			// 	for(i=0;i<Nc;i++){
			// 		ewaldself(i);
			// 	}
			// 	// ewaldrpy();
			// 	// gwdecompmod();
			// }
			// if(t%100000==0){
			// 	eps = 0.0;
			// 	for(i=0;i<3*N;i++){
			// 		rowsum[i] = 0.0;
			// 	}
			// 	ewaldrpy();
			// 	gwdecompmod();
			// 	printMM();
			// 	printTEA();
			// }
			if(t%sampf==0){
				ewaldrpy();
				gwdecompmod();
			}
			getNoise(); // Get new random velocities

			// updateChains(); // Update Positions GW Decomp
			if(itercount==0){
				updateFD();
			}
			else{
				// updateBins(); // Use bi HI for all bead pairs
				updateChains();
			}

			calcDis(); //
			if(t%printprops==0){
				calcCH();
				for(i=0;i<Ncharges*printprops;i++){
					dxt[i] = 0;
					dyt[i] = 0;
					dzt[i] = 0;
				}
				tcount=0;
			}
			// checkVerlet();

			//printf("t: %lu rx: %lf\n",t,rx[Charge_indices[0]]);
			//printperiod=1;
			if(t%printperiod==0){
				// printf("%lu\n",t);
				if(restart==0 || restart==2 || (restart==1 && t>tstart)){
					printTrajectory(); // Save trajectory to xyzfile every tau = 1/dt
				}
			}
			if(t%1000000==0){
				printMM();
				printTEA();
			}
			//printprops=1;
			if(t%printprops==0){
				if(restart == 0 || restart == 2 || (restart==1 && t>tstart)){
					calcRgCM(); // prints the radius of gyration and center of mass of each chain
				}
				// printf("%lu 1k ts in %lf seconds, Mult time %lf seconds, rx[100] = %lf\n",t,omp_get_wtime()-ksteps,ljtimetest,rx[100]);
				// ksteps = omp_get_wtime();
				// ljtimetest = 0.0;
			}
			tstart = 0;
	    }
		printf("%lu %.3e seconds \n",t,omp_get_wtime()-start_time);
		//printf("c1: %d c2: %d c3: %d c4: %d c5: %d c6: %d\n",counter1,counter2,counter3,counter4,counter5,counter6);
		calcMSD();
		printMM();
		printTEA();
		resetAverage(); // Reset the average on HI and TEA parameters and calculate the averages to be used in the next iteration
	}
    return 0;
}

void resetForce(){

	int i;
	for(i = 0; i<N; ++i){
		fx[i] = 0.0; fy[i] = 0.0; fz[i] = 0.0;
	}

}
float gasdev(long *idum){
	static int iset=0;
	static float gset;
	float fac,rsq,v1,v2;

	if (*idum < 0) iset = 0;
	if (iset == 0)
	{
		do
		{
			v1 = 2.0*ran1(idum)-1.0;
			v2 = 2.0*ran1(idum)-1.0;
			rsq = v1*v1+v2*v2;
		}
		while(rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	}
	else
	{
		iset = 0;
		return gset;
	}
}
void bondforce(){
	double dx1, dy1, dz1, rr, r, Fs, dx2, dy2, dz2, amp1, amp2, var, theta, coeff;
	int i, j;
	for(i = 0; i<Nc; ++i){
		for(j = 1; j < Nb; j++){
			dx1 = rx[Nb*i+j] - rx[Nb*i+j-1]; dy1 = ry[Nb*i+j] - ry[Nb*i+j-1]; dz1 = rz[Nb*i+j] - rz[Nb*i+j-1];
			dx1 -= box_length*round(dx1/box_length); dy1 -= box_length*round(dy1/box_length); dz1 -= box_length*round(dz1/box_length);
			rr = dx1*dx1+dy1*dy1+dz1*dz1;
			r = sqrt(rr);
			// if(r>4.0){
			// 	printf("%d %d\n",i,j);
			// 	exit(1);
			// }
			Fs = -kappas*(r-2.0)/r;
			fx[Nb*i+j] += Fs*dx1; fy[Nb*i+j] += Fs*dy1; fz[Nb*i+j] += Fs*dz1;
			fx[Nb*i+j-1] += -Fs*dx1; fy[Nb*i+j-1] += -Fs*dy1; fz[Nb*i+j-1] += -Fs*dz1;
		}
	}
	if(kappab>0.0){
		// printf("%lu\n",t);
		#pragma omp parallel
		{
			int i,j;
			double dx1, dy1, dz1, rr, r, Fs, dx2, dy2, dz2, amp1, amp2, var, theta, coeff;
			#pragma omp for schedule(static)
				for(i = 0; i<Nc; ++i){
					for(j = 1; j < Nb-1; j++){
						dx1 = rx[Nb*i+j] - rx[Nb*i+j-1]; dy1 = ry[Nb*i+j] - ry[Nb*i+j-1]; dz1 = rz[Nb*i+j] - rz[Nb*i+j-1];
						dx1 -= box_length*round(dx1/box_length); dy1 -= box_length*round(dy1/box_length); dz1 -= box_length*round(dz1/box_length);
						dx2 = rx[Nb*i+j] - rx[Nb*i+j+1]; dy2 = ry[Nb*i+j] - ry[Nb*i+j+1]; dz2 = rz[Nb*i+j] - rz[Nb*i+j+1];
						dx2 -= box_length*round(dx2/box_length); dy2 -= box_length*round(dy2/box_length); dz2 -= box_length*round(dz2/box_length);
						amp1 = sqrt(dx1*dx1+dy1*dy1+dz1*dz1);
						amp2 = sqrt(dx2*dx2+dy2*dy2+dz2*dz2);
						var = (dx1*dx2 + dy1*dy2 + dz1*dz2)/(amp1*amp2);
						theta = M_PI - acos(var); // Bond angle, Flory coordNb*i+jnates conventNb*i+jon
						coeff = -kappab*theta/sqrt(1-var*var); // theta_0 = 0 Nb*i+js the relaxed state
						#pragma omp atomic
							fx[Nb*i+j] += coeff*((dx1+dx2)/(amp1*amp2)-var*(dx1/(amp1*amp1)+dx2/(amp2*amp2)));
							fx[Nb*i+j-1] += coeff*(var*dx1/(amp1*amp1)-dx2/(amp1*amp2));
							fx[Nb*i+j+1] += coeff*(var*dx2/(amp2*amp2)-dx1/(amp1*amp2));
							fy[Nb*i+j] += coeff*((dy1+dy2)/(amp1*amp2)-var*(dy1/(amp1*amp1)+dy2/(amp2*amp2)));
							fy[Nb*i+j-1] += coeff*(var*dy1/(amp1*amp1)-dy2/(amp1*amp2));
							fy[Nb*i+j+1] += coeff*(var*dy2/(amp2*amp2)-dy1/(amp1*amp2));
							fz[Nb*i+j] += coeff*((dz1+dz2)/(amp1*amp2)-var*(dz1/(amp1*amp1)+dz2/(amp2*amp2)));
							fz[Nb*i+j-1] += coeff*(var*dz1/(amp1*amp1)-dz2/(amp1*amp2));
							fz[Nb*i+j+1] += coeff*(var*dz2/(amp2*amp2)-dz1/(amp1*amp2));
					}
				}
		}
	}
}
void LJforce(){
	double dx, dy, dz, rr, coeff, r6, ratio;
	int i, j, k;
	// #pragma omp parallel for schedule(static,1) private(j,k,dx,dy,dz,rr,coeff,r6,ratio) // not thread safe, need critical section
	#pragma omp parallel
	{
		int j,k;
		double dx,dy,dz,rr,coeff,r6,ratio;
		double *fxt = calloc(N,sizeof(double));
		double *fyt = calloc(N,sizeof(double));
		double *fzt = calloc(N,sizeof(double));
		#pragma omp for schedule(static,1)
			for(i = 0; i<N; ++i){
				for(j = 0; j<nlist[i]; ++j){
					k = list[i*N+j];
					dx = rx[i] - rx[k]; dy = ry[i] - ry[k]; dz = rz[i] - rz[k];
					dx -= box_length*round(dx/box_length); dy -= box_length*round(dy/box_length); dz -= box_length*round(dz/box_length);
					rr = dx*dx+dy*dy+dz*dz;
					if(rr<r2cut){
						ratio = 4.00/rr;
						r6 = ratio*ratio*ratio;
						if(r6>3) r6 = 3;
						// coeff = (48*epsilon/rr)*(r6*r6-0.5*r6);
						coeff = (12*epsilon/rr)*(r6*r6-r6);
						// #pragma omp atomic
							fxt[i] += coeff*dx;
							fyt[i] += coeff*dy;
							fzt[i] += coeff*dz;
							fxt[k] -= coeff*dx;
							fyt[k] -= coeff*dy;
							fzt[k] -= coeff*dz;
					}
				}
			}
		#pragma omp critical
		{
			for(i=0;i<N;i++){
				fx[i] += fxt[i];
				fy[i] += fyt[i];
				fz[i] += fzt[i];
			}
			free(fxt); free(fyt); free(fzt);
		}
	}
}
void coulombForce(){
	double dx,dy,dz,rr,r,coeff,r6,ratio,Fc;
	int i,j,k,l;

	// for(i=0;i<N;i++){
	// 	printf("Charge[i]: %d \n", Charge[i]) ;
	// }
	// printf("lambda_d: %lf \n", lambda_d);
	kappa_d = 1/lambda_d ;

	for (i=0; i<N; ++i){
		for (j=i+1; j<N; ++j){
			if (Charge[i]!=0 && Charge[j]!=0)/// If Bead is charged and testing the next bead neighboring bead
			{
				dx = rx[i] - rx[j]; dy = ry[i] - ry[j]; dz = rz[i] - rz[j];
				dx -= box_length*round(dx/box_length); dy -= box_length*round(dy/box_length); dz -= box_length*round(dz/box_length);
				rr = dx*dx+dy*dy+dz*dz;
				r = sqrt(rr) ;


				Fc = lambda_b*Charge[i]*Charge[j]*exp(-kappa_d*r)*(1+r*kappa_d)/rr;

				E_initial += (lambda_b*Charge[i]*Charge[j]/r)*exp(-kappa_d*r); // Calculating for Initial Energy
				//printf("%lf \n",Fc);

				fx[i] += Fc*dx/r;
				fy[i] += Fc*dy/r;
				fz[i] += Fc*dz/r;
				fx[j] -= Fc*dx/r;
				fy[j] -= Fc*dy/r;
				fz[j] -= Fc*dz/r;


			}
		}
	}
}
void chargeHop(){

	int i, j, k,c,l, test, test2, cnumber, cnumber2, count, NCharges,num_of_particles;
	double E_final, dx, dy, dz, r, rr , mm, dice_roll;
	double rxtempi, rxtempf,rytempi,rytempf,rztempi,rztempf;
	rxtempi = 0, rxtempf = 0, rytempi = 0, rytempf = 0, rztempi = 0, rztempf = 0;

	NCharges = Ncharges;
	num_of_particles = Nb*Nc ;
	int num_of_p_particles = Nb*Nc ;

	if (t%MStep==0) // Do Monte Carlo Move every MSteps
	{
		for (j=0;j<NCharges;++j)
		{
			for (i=0; i<NCharges;++i)
			{
				if (Charge_track[j]==Charge_indices[i] && i==j)
				{i_values[i]=j;}
			}
		}

		for(i=0;i<NCharges;++i)
		{	Charge_indices[i]=Charge_track[i];}

	for (k=0; k<NCharges; ++k) // Do MC_Move for each charge
	{
		for(i=0;i<num_of_p_particles;++i)
		{
			Charge[i]=0;
		}

		for(i=0;i<num_of_p_particles;++i)
		{
			for(j=0;j<NCharges;++j)
			{
				if(i==Charge_track[j])
				Charge[i]=1;

			}
		}

		//E_initial calculated in coulombforce()//

		mm=ran1(idum);
		cnumber=Charge_indices[k]; //find the next charge on the shuffled list
		test=0;// test for am I close?
		test2=0;

		for(j = 0; j<num_of_p_particles; ++j) //for loop to check distance between cnumber and every other particle
					{
						if (j!=cnumber && j!=cnumber+1 && j!=cnumber-1)
						{
							dx = rx[cnumber] - rx[j]; dy = ry[cnumber] - ry[j]; dz = rz[cnumber] - rz[j];
							dx -= box_length*round(dx/box_length); dy -= box_length*round(dy/box_length); dz -= box_length*round(dz/box_length);

								if(dx*dx+dy*dy+dz*dz<4.41)
								{
					if (mm<.33) //Move to the right
					{
						cnumber2=cnumber+1; //right hand neighbor
						//printf("t=%d RR cnumber=%d cnumber2=%d\n",t,cnumber,cnumber2);
						for(i=0;i<NCharges;++i)
						{

							if((cnumber==(num_of_particles-1) || cnumber%(Nb)==Nb-1 || cnumber==Nb-1 || Charge_indices[i]==cnumber2)) // don't check if right hand neighbor is charged
							{
								test2=1;
							}
						}
						if (test2!=0)
						{
							//printf("RR Skip\n");
							continue;
						}
						for(i=0;i<NCharges;++i)
						{
							if (Charge_indices[i]==cnumber )
							{
								Charge_indices[i]=cnumber2;
								//printf("RR Swap\n");
								test2=0;
							}
						}

						test=1;//possibility of hopping to NOT neighbor occurred

					}
					if(mm>0.66) //move to the left - if random # is less than .5
					{
						cnumber2=cnumber-1; //left hand neighbor
						for (i=0;i<NCharges;++i)
						{
							if(cnumber2<0 || cnumber==0 || cnumber%Nb==0 || Charge_indices[i]==cnumber2)// don't attempt to go left at bead 0
							{
								test2=1;

							}
						}
						if (test2!=0)
						{
							//printf("LL Skip\n");
								continue;
						}
						else{
							for(i=0;i<NCharges;++i)
							{
								if (Charge_indices[i]==cnumber )
								{	Charge_indices[i]=cnumber2;
									//printf("LL Swap\n");
									test2=0;
								}
							}
							test=1;//possibility of hopping to NOT neighbor occurred
						}

					}
					else if(mm<=0.66 && mm >=0.33) ///swap with j particle that's close to cnumber
					{
						cnumber2=j;
						//printf("t=%d P cnumber=%d cnumber2=%d\n",t,cnumber,cnumber2);
						for(i=0;i<NCharges;++i)
						{
							if (Charge_indices[i]==cnumber2)
							{
								test2=1;
							}
						}
						if (test2!=0)
						{
							//printf("Pop Skip\n");
								continue;
						}
						else{
							for(i=0;i<NCharges;++i)
							{
								if (Charge_indices[i]==cnumber)
								{
									//printf("Pop\n");  //Pop=non-adjacent hop
									Charge_indices[i]=cnumber2;
									test2=0;
								}
							}
							test=2;
						}
					}
					break;
				}
			}
		}

	//////////// Perform swap if NOT close enough, only to a neighbor bead///////////////
		if(test==0){
				if (mm>.5) //Move to the right
				{

				cnumber2=cnumber+1; //right hand neighbor
				//printf("t=%d R cnumber=%d cnumber2=%d\n",t,cnumber,cnumber2);

					for(i=0;i<NCharges;++i)
					{

						if(cnumber==(num_of_particles-1) || cnumber%(Nb)==Nb-1 || cnumber==Nb-1 || Charge_indices[i]==cnumber2) // don't check if right hand neighbor is charged
						{
							test2=1;
						}
					}
					if (test2!=0)
					{
						//printf("R Skip\n");
							continue;
					}
					else{
						for (i=0;i<NCharges;++i)
						{
							if (Charge_indices[i]==cnumber )
								{
									Charge_indices[i]=cnumber2;
									//printf("R Swap\n");
									test2=0;
								}
						}
					}

				}
				else //move to the left - if random # is less than .5
				{

					cnumber2=cnumber-1; //left hand neighbor
					//printf("t=%d L cnumber=%d cnumber2=%d\n",t,cnumber,cnumber2);
					for (i=0;i<NCharges;++i)
					{
						if(cnumber2<0 || cnumber==0 || cnumber%Nb==0 || Charge_indices[i]==cnumber2)// don't attempt to go left at bead 0
						{
							test2=1;
						}
					}
					if (test2!=0)
					{
						//printf("L Skip\n");
							continue;
					}
					else{
						for (i=0;i<NCharges;++i)
						{
							if (Charge_indices[i]==cnumber)
							{
								Charge_indices[i]=cnumber2;
								//printf("L Swap\n");
								test2=0;
								test=0;

							}
						}
					}
				}
			}

			// Update the charge positions to test new energy
			for(i=0;i<num_of_p_particles;++i)
			{
				Charge[i]=0;
			}

			for(i=0;i<num_of_p_particles;++i)
			{
				for(j=0;j<NCharges;++j)
				{
					if(i==Charge_indices[j])
					Charge[i]=1;

				}
			}

			//Calculate New Energy
			E_final=0;
			for (i=0; i<num_of_particles; ++i)
			{
				for (j=i+1; j<num_of_particles; ++j)
						{
							if (Charge[i]!=0 && Charge[j]!=0)/// If Bead is charged and testing the next bead neighboring bead
							{
								dx = rx[i] - rx[j]; dy = ry[i] - ry[j]; dz = rz[i] - rz[j];
								dx -= box_length*round(dx/box_length); dy -= box_length*round(dy/box_length); dz -= box_length*round(dz/box_length);
								rr = dx*dx+dy*dy+dz*dz;
								r = sqrt(rr) ;
								E_final += (lambda_b*Charge[i]*Charge[j]/r)*exp(-kappa_d*r); // Calculating for Final Energy
							}
						}
			}
			// printf("this is Efinal=");
			// printf("%lf\n",E_final);
			 //printf("this delta E=(%lf)\n", (E_final-E_initial));

		 if (test==0 || test==1) // Perform MC but under neighbor conditions - barrier
		 {
		 //Actual Monte Carlo Test
			 if (1<=exp((-barrier-0.5*(E_final - E_initial))))
			 {
				 count=0;
								 for(i = 0; i<NCharges; ++i)
								 {
										 if(Charge_track[i]==cnumber && i==i_values[i] && count==0 )
										 { 	count+=1;
											 Charge_old[i]=cnumber;
											 Charge_indices[i]=cnumber2;
												 Charge_track[i]=cnumber2;
								 // 				Output2 = fopen(str2, "a");
								 // 				fprintf(Output2,"%d %d %d %d %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]]);
								 // // // printf("%d %d %d %d %lf %lf %lf %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]], Chainrx[Charge_old[i]], Chainry[Charge_old[i]],Chainrz[Charge_old[i]]);
								 // fclose(Output2);
										 }
								 }
								 counter1++ ;
			 }
			 else
			 {
					 dice_roll=ran1(idum);
					 if (dice_roll<=exp((-barrier-0.5*(E_final - E_initial)))) // if greater, accept swap
					 {	count=0;
												 for(i = 0; i<NCharges; ++i)
												 {
														 if(Charge_track[i]==cnumber && i==i_values[i] && count==0)
														 {
															 count+=1;
															 Charge_old[i]=cnumber;
															 Charge_indices[i]=cnumber2;
																 Charge_track[i]=cnumber2;
											 // 				Output2 = fopen(str2, "a");
											 // 				fprintf(Output2,"%d %d %d %d %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]]);
											 // fclose(Output2);
											 // printf("%d %d %d %d %lf %lf %lf %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]], Chainrx[Charge_old[i]], Chainry[Charge_old[i]],Chainrz[Charge_old[i]]);
														 }
												 }
												 counter2++ ;
					 }
					 else
					 {	count=0;
						 for(i=0;i<NCharges; ++i)
						 {
							 if (Charge_indices[i]==cnumber2 && i==i_values[i] && count==0)
							 {	count+=1;
								 // printf("count=%d tt=%d Charge_indices=%d\n",count,t,cnumber2);
								 Charge_indices[i]=cnumber;
								 Charge_old[i]=cnumber;
								 }
						 }
						 counter3++ ;
					 }
				 }
			 }
			 ///////// If there's a possibility of a bead (not neighbor) close by to hop a charge
			 else if(test==2)
			 {
				 if (1<=exp((-barrier2-0.5*(E_final - E_initial))))
				 {
					 count=0;
					 for(i = 0; i<NCharges; ++i)
					 {
							 if(Charge_track[i]==cnumber && i==i_values[i] && count==0 )
							 { 	count+=1;
								 Charge_old[i]=cnumber;
								 Charge_indices[i]=cnumber2;
									 Charge_track[i]=cnumber2;
					 // 				Output2 = fopen(str2, "a");
					 // 				fprintf(Output2,"%d %d %d %d %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]]);
					 // // // printf("%d %d %d %d %lf %lf %lf %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]], Chainrx[Charge_old[i]], Chainry[Charge_old[i]],Chainrz[Charge_old[i]]);
					 // fclose(Output2);
							 }
					 }
					 counter4++ ;
				 }
				 else{
					 dice_roll=ran1(idum);
					 if (dice_roll<=exp((-barrier2-0.5*(E_final - E_initial)))) // if greater, accept swap
					 {	count=0;
												 for(i = 0; i<NCharges; ++i)
												 {
														 if(Charge_track[i]==cnumber && i==i_values[i] && count==0)
														 {
															 count+=1;
															 Charge_old[i]=cnumber;
															 Charge_indices[i]=cnumber2;
																 Charge_track[i]=cnumber2;
											 // 				Output2 = fopen(str2, "a");
											 // 				fprintf(Output2,"%d %d %d %d %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]]);
											 // fclose(Output2);
											 // printf("%d %d %d %d %lf %lf %lf %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]], Chainrx[Charge_old[i]], Chainry[Charge_old[i]],Chainrz[Charge_old[i]]);
														 }
												 }
												 counter5++ ;
					 }
					 else
					 {	count=0;
						 for(i=0;i<NCharges; ++i)
						 {
							 if (Charge_indices[i]==cnumber2 && i==i_values[i] && count==0)
							 {	count+=1;
								 // printf("count=%d tt=%d Charge_indices=%d\n",count,t,cnumber2);
								 Charge_indices[i]=cnumber;
								 Charge_old[i]=cnumber;
								 }
						 }
						 counter6++ ;
					 }
				 }
			 }

		}//end of k loop
	}//end of MC step
	for(i=0;i<num_of_p_particles;++i)
	{
		Charge[i]=0;
	}

	for(i=0;i<num_of_p_particles;++i)
	{
		for(j=0;j<NCharges;++j)
		{
			if(i==Charge_track[j])
			Charge[i]=1;
		}
	}

}//end of function
void chargeHop_nonint(){

	int i, j, k,c,l, test, test2, cnumber, cnumber2, count, NCharges,num_of_particles;
	double E_initial, E_final, dx, dy, dz, r, rr , mm, dice_roll;
	double rxtempi, rxtempf,rytempi,rytempf,rztempi,rztempf;
	rxtempi = 0, rxtempf = 0, rytempi = 0, rytempf = 0, rztempi = 0, rztempf = 0;
	for(i=0;i<Ncharges;i++){
		check_ch[3*i]=0;
		check_ch[3*i+1]=0;
		check_ch[3*i+2]=0;
	}

	NCharges = Ncharges;
	num_of_particles = N ;

	if (t%MStep==0) // Do Monte Carlo Move every MSteps
	{
		for (j=0;j<NCharges;++j)
		{
			for (i=0; i<NCharges;++i)
			{
				if (Charge_track[j]==Charge_indices[i] && i==j)
				{i_values[i]=j;}
			}
		}

		for(i=0;i<NCharges;++i)
		{	Charge_indices[i]=Charge_track[i];}

	for (k=0; k<NCharges; ++k) // Do MC_Move for each charge
	{		mm=ran1(idum);
			//printf("mm=%f\n",mm);

			cnumber=Charge_indices[k]; //find the next charge on the shuffled list
			test=0;// test for am I close?
			test2=0;
			for(j = 0; j<num_of_particles; ++j) //for loop to check distance between cnumber and every other particle
						{
							if (j!=cnumber && j!=cnumber+1 && j!=cnumber-1)
							{
								dx = rx[cnumber] - rx[j]; dy = ry[cnumber] - ry[j]; dz = rz[cnumber] - rz[j];
								dx -= box_length*round(dx/box_length); dy -= box_length*round(dy/box_length); dz -= box_length*round(dz/box_length);


									if(dx*dx+dy*dy+dz*dz<4.41)
									{
						if (mm<.33) //Move to the right
						{
							cnumber2=cnumber+1; //right hand neighbor
							//printf("t=%d RR cnumber=%d cnumber2=%d\n",t,cnumber,cnumber2);
							if(cnumber==(num_of_particles-1) || cnumber%(Nb)==Nb-1 || cnumber==Nb-1) // don't check if right hand neighbor is charged
							{
								continue;
							}

							for(i=0;i<NCharges;++i)
							{
								if (Charge_indices[i]==cnumber )
								{
									Charge_indices[i]=cnumber2;
									//printf("RR Swap\n");
								}
							}

							test=1;//possibility of hopping to NOT neighbor occurred

						}
						if(mm>0.66) //move to the left - if random # is less than .5
						{
							cnumber2=cnumber-1; //left hand neighbor

							if(cnumber==0 || cnumber%Nb==0)// don't attempt to go left at bead 0
							{
								continue;
							}

							for(i=0;i<NCharges;++i)
							{
								if (Charge_indices[i]==cnumber )
								{	Charge_indices[i]=cnumber2;
									//printf("LL Swap\n");
								}
							}
							test=1;//possibility of hopping to NOT neighbor occurred
						}
						else if(mm<=0.66 && mm >=0.33) ///swap with j particle that's close to cnumber
						{
							cnumber2=j;
							//printf("t=%d P cnumber=%d cnumber2=%d\n",t,cnumber,cnumber2);

							for(i=0;i<NCharges;++i)
							{
								if (Charge_indices[i]==cnumber)
								{
									//printf("Pop\n");  //Pop=non-adjacent hop
									Charge_indices[i]=cnumber2;
								}
							}
							test=2;
						}
						break;
					}
				}
			}


//////////// Perform swap if NOT close enough, only to a neighbor bead///////////////
		if(test==0)
		{
			if (mm>.5) //Move to the right
			{

			cnumber2=cnumber+1; //right hand neighbor
			//printf("t=%d R cnumber=%d cnumber2=%d\n",t,cnumber,cnumber2);

				if(cnumber==(num_of_particles-1) || cnumber%(Nb)==Nb-1 || cnumber==Nb-1) // don't check if right hand neighbor is charged
				{
					continue;
				}

				for (i=0;i<NCharges;++i)
				{
					if (Charge_indices[i]==cnumber )
						{
							Charge_indices[i]=cnumber2;
							//printf("R Swap\n");
						}
				}
			}
			else //move to the left - if random # is less than .5
			{

				cnumber2=cnumber-1; //left hand neighbor
				//printf("t=%d L cnumber=%d cnumber2=%d\n",t,cnumber,cnumber2);

				if(cnumber==0 || cnumber%Nb==0)// don't attempt to go left at bead 0
				{
					continue;
				}

				for (i=0;i<NCharges;++i)
				{
					if (Charge_indices[i]==cnumber)
					{
						Charge_indices[i]=cnumber2;
						//printf("L Swap\n");
					}
				}
			}
		}

		if (test==0 || test==1) // Perform MC but under neighbor conditions - barrier
		{


			//Actual Monte Carlo Test
			E_final = 0;
			E_initial = 0;
			if (1<=exp((-barrier-0.5*(E_final - E_initial)))) //if energy change is greater than 1 accept change
			{	count=0;
									for(i = 0; i<NCharges; ++i)
									{
											if(Charge_track[i]==cnumber && i==i_values[i] && count==0 )
											{ 	count+=1;
												Charge_old[i]=cnumber;
												Charge_indices[i]=cnumber2;
													Charge_track[i]=cnumber2;
									// 				Output2 = fopen(str2, "a");
									// 				fprintf(Output2,"%d %d %d %d %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]]);
									// // // printf("%d %d %d %d %lf %lf %lf %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]], Chainrx[Charge_old[i]], Chainry[Charge_old[i]],Chainrz[Charge_old[i]]);
									// fclose(Output2);
											}
									}
									counter1++ ;
			}

			else
			{
				dice_roll=ran1(idum);
				if (dice_roll<=exp((-barrier-0.5*(E_final - E_initial)))) // if greater, accept swap
				{	count=0;
											for(i = 0; i<NCharges; ++i)
											{
													if(Charge_track[i]==cnumber && i==i_values[i] && count==0)
													{
														count+=1;
														Charge_old[i]=cnumber;
														Charge_indices[i]=cnumber2;
														Charge_track[i]=cnumber2;
										// 				Output2 = fopen(str2, "a");
										// 				fprintf(Output2,"%d %d %d %d %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]]);
										// fclose(Output2);
										// printf("%d %d %d %d %lf %lf %lf %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]], Chainrx[Charge_old[i]], Chainry[Charge_old[i]],Chainrz[Charge_old[i]]);
													}
											}
											counter2++ ;
				}
				else
				{	count=0;
					for(i=0;i<NCharges; ++i)
					{
						if (Charge_indices[i]==cnumber2 && i==i_values[i] && count==0)
						{	count+=1;
							// printf("count=%d tt=%d Charge_indices=%d\n",count,t,cnumber2);
							Charge_indices[i]=cnumber;
							}
					}
					counter3++ ;
				}
			}
		}

		///////// If there's a possibility of a bead (not neighbor) close by to hop a charge


		else if(test==2) ///Perform MC under NOT neighbor conditions - barrier2
		{
			//Actual Monte Carlo Test
			if (1<=exp((-barrier2-0.5*(E_final - E_initial)))) //if energy change is greater than 1 accept change
			{
				count=0;
									for(i = 0; i<NCharges; ++i)
									{
											if(Charge_track[i]==cnumber && i==i_values[i] && count==0)
											{	count+=1;

												Charge_old[i]=cnumber;
													Charge_track[i]=cnumber2;
													Charge_indices[i]=cnumber2;
									// 				Output2 = fopen(str2, "a");
									// 	fprintf(Output2,"%d %d %d %d %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]]);
									// fclose(Output2);
									// printf("%d %d %d %d %lf %lf %lf %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]], Chainrx[Charge_old[i]], Chainry[Charge_old[i]],Chainrz[Charge_old[i]]);

											}
									}
									counter4++ ;
			}

			else
			{
				dice_roll=ran1(idum);
				if (dice_roll<=exp((-barrier2-0.5*(E_final - E_initial)))) // if greater, accept swap
				{
					count=0;
											for(i = 0; i<NCharges; ++i)
											{
													if(Charge_track[i]==cnumber && i==i_values[i] && count==0)
													{count+=1;
														Charge_old[i]=cnumber;
															Charge_track[i]=cnumber2;
															Charge_indices[i]=cnumber2;
										// 					Output2 = fopen(str2, "a");
										// 					fprintf(Output2,"%d %d %d %d %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]]);
										// fclose(Output2);
										// printf("%d %d %d %d %lf %lf %lf %lf %lf %lf\n",t, i,Charge_track[i],Charge_old[i], Chainrx[Charge_track[i]], Chainry[Charge_track[i]],Chainrz[Charge_track[i]], Chainrx[Charge_old[i]], Chainry[Charge_old[i]],Chainrz[Charge_old[i]]);
													}
											}
											counter5++ ;
				}
				else
				{counter6++ ; //reject and swap back

					count=0;
					for(i=0;i<NCharges; ++i)
					{
						if (Charge_indices[i]==cnumber2 && i==i_values[i] && count==0)
						{count+=1;
							// printf("t= %d Charge_indices=%d\n",t,Charge_indices[i]);
						Charge_indices[i]=cnumber;
						}
					}

				}
			}
		}


	}//end of k loop
}//end of MC step


}//end of function
void calcDis()
{
	double dx, dy, dz,L;
	int i,j;

	L = box_length;
//int c = 0;
	//outputfile2=fopen("disp1.csv","a");
	tcount++ ;
	for(i=0;i<Ncharges;i++){
		j = Charge_track[i];
		dx = rx[j] - rbi[3*i]  ;
		dy = ry[j] - rbi[3*i+1]  ;
		dz = rz[j] - rbi[3*i+2]  ;

		if(dx>L/2){dx-=L;} if(dx<-L/2){dx+=L;}
		if(dy>L/2){dy-=L;} if(dy<-L/2){dy+=L;}
		if(dz>L/2){dz-=L;} if(dz<-L/2){dz+=L;}

		dxt[i*printprops+tcount] = dx  ;
		dyt[i*printprops+tcount] = dy  ;
		dzt[i*printprops+tcount] = dz  ;
	}

}
void calcCH(){
	int i,j,ind2;
	double dx,dy,dz;

	chfile = fopen(str2,"a");
	//printf("t: %lu\n",t);
	if(t==0){
		dx=0;
		dy=0;
		dz=0;
		for(j=0;j<Ncharges;j++){
			fprintf(chfile,"%lu %d %d %d %lf %lf %lf\n",t,j,Charge_old[j],Charge_track[j],dx,dy,dz);
		}
	}
	else{
		for(j=0;j<Ncharges;j++){
			dx=0;
			dy=0;
			dz=0;
			for(i=1;i<printprops+1;i++){
					ind2 = j*printprops+i ;
					//printf("j: %d t: %d dxt[%d]: %lf\n",j,i,ind2,dxt[ind2]);
					dx += dxt[ind2] ;
					dy += dyt[ind2] ;
					dz += dzt[ind2] ;
				}
			//printf("t: %lu j: %d dx[0]: %lf\n",t,j,dx);
			fprintf(chfile,"%lu %d %d %d %lf %lf %lf\n",t,j,Charge_old[j],Charge_track[j],dx,dy,dz);
		}
	}
	fclose(chfile);
}
void initPos(){
	int i,j;
	for(i=0;i<Ncharges;i++){
		j = Charge_track[i];
		rbi[3*i] = rx[j];
		rbi[3*i+1] = ry[j];
		rbi[3*i+2] = rz[j];
	}
}
void verletlist(){
	int i, j;
	double dx, dy, dz, r,start_time;
	for(i=0;i<N;i++){
		nlist[i] = 0; // Reset lists on update
		xo[i] = rx[i]; yo[i] = ry[i]; zo[i] = rz[i];
		// printf("%d %lf %lf %lf\n",i,xo[i],yo[i],zo[i]);
	}
	// Only loop over each pair once. Each pair of neighbors only needs to be in one table or the other, not both
	// start_time = omp_get_wtime();
	// omp_set_num_threads(4);
	#pragma omp parallel
	{
		// int id = omp_get_thread_num();
		int i,j;
		double dx,dy,dz,r;
		#pragma omp for schedule(static,1)
			for(i=0;i<N-1;i++){
				for(j=i+1;j<N;j++){
					dx = rx[i] - rx[j]; dy = ry[i] - ry[j]; dz = rz[i] - rz[j];
					dx -= box_length*round(dx/box_length); dy -= box_length*round(dy/box_length); dz -= box_length*round(dz/box_length);
					r = sqrt(dx*dx+dy*dy+dz*dz);
					// printf("%d %d %d\n",id,i,j);
					// Check - if particle j within rv or particle i, add to NL
					if(r < rv){
						list[i*N+nlist[i]] = j;
						// list[j*N+nlist[j]] = i;
						#pragma omp atomic
							nlist[i]++;
							// nlist[j]++;
					}
				}
			}
	}
	// printf("Verlet list time %lf seconds\n",omp_get_wtime()-start_time);
	// for(i=0;i<nlist[100];i++){
	// 	printf("%d\n",list[100*N+i]);
	// }
	// exit(1);
}
int binmod(int i, int j){
	int bin,grid,gridx,gridy,gridz,bin1x,bin1y,bin1z,bin2x,bin2y,bin2z;
	double dx,dy,dz,dxabs,dyabs,dzabs;
	dx = rx[j] - rx[i]; dy = ry[j] - ry[i]; dz = rz[j] - rz[i];
	dx -= box_length*round(dx/box_length); dy -= box_length*round(dy/box_length); dz -= box_length*round(dz/box_length);
	// res.x = round(dx/bin_size)+num_g/2;
	// res.y = round(dy/bin_size)+num_g/2;
	// res.z = round(dz/bin_size)+num_g/2;
	// grid = (1 - round((dx+0.5*gs2)/gl1))*(1 - round((dy+0.5*gs2)/gl1))*(1 - round((dz+0.5*gs2)/gl1));
	if(fabs(dx)<gl1 && fabs(dy)<gl1 && fabs(dz)<gl1){
		bin1x = round(dx/gs1) + gh1;
		bin1y = round(dy/gs1) + gh1;
		bin1z = round(dz/gs1) + gh1;
		bin = bin1x*gd1*gd1 + bin1y*gd1 + bin1z;
		// printf("i %d j %d Grid 1 dx %lf dy %lf dz %lf bin1x %d bin1y %d bin1z %d bin %d\n",i,j,fabs(dx),fabs(dy),fabs(dz),bin1x,bin1y,bin1z,bin);
	}
	else{
		bin2x = round(dx/gs2) + gh2;
		bin2y = round(dy/gs2) + gh2;
		bin2z = round(dz/gs2) + gh2;
		bin = gn1 + bin2x*gd2*gd2 + bin2y*gd2 + bin2z;
		// printf("i %d j %d Grid 2 dx %lf dy %lf dz %lf bin2x %d bin2y %d bin2z %d bin %d\n",i,j,fabs(dx),fabs(dy),fabs(dz),bin2x,bin2y,bin2z,bin);
	}
	// bin = grid*(bin1x*gd1*gd1 + bin1y*gd1 + bin1z) + (1-grid)*(gn1 + bin2x*gd2*gd2 + bin2y*gd2 + bin2z);
	// printf("i %d j %d dx %lf dy %lf dz %lf grid %d bin1x %d bin1y %d bin1z %d bin2x %d bin2y %d bin2z %d bin %d\n",i,j,dx,dy,dz,grid,bin1x,bin1y,bin1z,bin2x,bin2y,bin2z,bin);
	return bin;
}
void ewaldrpy(){
	int i,j,kx,ky,kz,kk;
	double rkx,rky,rkz,rkk;
	double start_time = omp_get_wtime();
	omp_set_num_threads(num_threads);
	#pragma omp parallel
	{
		int id = omp_get_thread_num();
		double epsp,eps2p;
		// double *Dijp = calloc(9*tot_bins,sizeof(double));
		// double *rowsump = calloc(3*N,sizeof(double));
		int i,j,k,l,m,n,a,b,startij,start,startji,binij,binji;
		int kx,ky,kz,kk;
		double dx, dy, dz, r, rr, C1, C2, C3, C4, coeff1, coeff2, term1, term2, coskr,incr,rkx,rky,rkz;
		double Dijtemp[3][3];
		double M1[3][3];
		// struct bin tbij; struct bin tbji;
		epsp = 0.0; eps2p = 0.0;
		#pragma omp for schedule(static,1)
			for(i=0;i<N;i++){
				for(j=((int)(i/Nb)+1)*Nb;j<N;j++){
					// tbij = binij(i,j);
					// tbji = binij(j,i);
					// printf("%d %d %d %d %d\n",i,j,tbij.x,tbij.y,tbij.z);
					binij = binmod(i,j);
					binji = binmod(j,i);
					// startij = tbij.x*num_g*num_g+tbij.y*num_g+tbij.z;
					// startji = tbji.x*num_g*num_g+tbji.y*num_g+tbji.z;
					count[binij]++;
					count[binji]++;
					for(a=0;a<3;a++){
						for(b=0;b<3;b++){
							Dijtemp[a][b] = 0.0;
						}
					}
					dx = rx[i] - rx[j]; dy = ry[i] - ry[j]; dz = rz[i] - rz[j];
					dx -= box_length*round(dx/box_length); dy -= box_length*round(dy/box_length); dz -= box_length*round(dz/box_length);
					rr = dx*dx+dy*dy+dz*dz; r = sqrt(rr);
					// Real Part
					if(r>2.0){
						C1 = 0.75/r+0.5/(r*rr);
						C2 = 4.0*alpha7*rr*rr + 3.0*alpha3*rr - 20.0*alpha5*rr - 4.5*alpha + 14.0*alpha3 + alpha/rr;
						C3 = 0.75/r-1.5/(r*rr);
						C4 = -4.0*alpha7*rr*rr - 3.0*alpha3*rr + 16.0*alpha5*rr + 1.5*alpha - 2.0*alpha3 - 3.0*alpha/rr;
					}
					else{
						C1 = 1.0 - 0.28125*r;
						C2 = 64.0*alpha7 + 12.0*alpha3 - 80.0*alpha5 - 4.5*alpha + 14.0*alpha3 + 0.25*alpha;
						C3 = 0.09375*r;
						C4 = -64.0*alpha7 - 12.0*alpha3 + 64.0*alpha5 + 1.5*alpha - 2.0*alpha3 - 0.75*alpha;
					}
					term1 = C1*erfc(alpha*r) + C2*exp(-alpha2*rr)/rtpi;
					term2 = C3*erfc(alpha*r) + C4*exp(-alpha2*rr)/rtpi;
					M1[0][0] = term1 + term2*dx*dx/rr;
					M1[0][1] = term2*dx*dy/rr;
					M1[0][2] = term2*dx*dz/rr;
					M1[1][0] = M1[0][1];
					M1[1][1] = term1 + term2*dy*dy/rr;
					M1[1][2] = term2*dy*dz/rr;
					M1[2][0] = M1[0][2];
					M1[2][1] = M1[1][2];
					M1[2][2] = term1 + term2*dz*dz/rr;
					Dijtemp[0][0] += M1[0][0];
					Dijtemp[0][1] += M1[0][1];
					Dijtemp[0][2] += M1[0][2];
					Dijtemp[1][0] += M1[1][0];
					Dijtemp[1][1] += M1[1][1];
					Dijtemp[1][2] += M1[1][2];
					Dijtemp[2][0] += M1[2][0];
					Dijtemp[2][1] += M1[2][1];
					Dijtemp[2][2] += M1[2][2];
					// Loop over wavevectors for Drecip
					for(kx = -kmax; kx < kmax+1; kx++){
						rkx = kx*kcoeff;
						for(ky = -kmax; ky < kmax+1; ky++){
							rky = ky*kcoeff;
							for(kz = -kmax; kz < kmax+1; kz++){
								rkz = kz*kcoeff;
								kk = kx*kx + ky*ky + kz*kz;
								if(kk != 0){
									coskr = cos(rkx*dx + rky*dy + rkz*dz);
									start = 9*((kx+kmax)*numk*numk+(ky+kmax)*numk+(kz+kmax));
									Dijtemp[0][0] += M2[start]*coskr;
									Dijtemp[0][1] += M2[start+1]*coskr;
									Dijtemp[0][2] += M2[start+2]*coskr;
									Dijtemp[1][0] += M2[start+3]*coskr;
									Dijtemp[1][1] += M2[start+4]*coskr;
									Dijtemp[1][2] += M2[start+5]*coskr;
									Dijtemp[2][0] += M2[start+6]*coskr;
									Dijtemp[2][1] += M2[start+7]*coskr;
									Dijtemp[2][2] += M2[start+8]*coskr;
								}
							}
						}
					}
					for(a=0;a<3;a++){
						for(b=0;b<3;b++){
							incr = Dijtemp[a][b]/Dtr;
							epsp += 2.0*incr;
							eps2p += 2.0*incr*incr;
							#pragma omp atomic
								Dijavg[9*binij+3*a+b] += Dijtemp[a][b];
								Dijavg[9*binji+3*a+b] += Dijtemp[a][b];
								rowsum[3*i+a] += incr*incr;
								rowsum[3*j+a] += incr*incr;
							// Dijp[9*startij+3*a+b] += Dijtemp[a][b];
							// Dijp[9*startji+3*a+b] += Dijtemp[a][b];
							// rowsump[3*i+a] += incr*incr;
							// rowsump[3*j+a] += incr*incr;
						}
					}
				}
			}
		#pragma omp critical
		{
			// printf("%d %lf\n",id,Dijp[9*(num_bin+1)*(num_bin+1)*(num_bin+1)]);
			eps += epsp; eps2 += eps2p;
			// for(i=0;i<9*tot_bins;i++){
			// 	Dijavg[i] += Dijp[i];
			// }
			// for(i=0;i<3*N;i++){
			// 	rowsum[i] += rowsump[i];
			// }
			// free(Dijp); free(rowsump);
		}
	}
	// printf("%lu %lf\n",t,rowsum[0]);
	for(i=0;i<Nc;i++){
		ewaldself(i);
	}
	// printf("Time Step %lu Ewald sum time %lf eps = %lf Dij[0] = %lf Rowsum[0] = %lf Num threads = %d\n",t,omp_get_wtime()-start_time,eps,Dijavg[0],rowsum[0],omp_get_num_threads());
	// exit(1);
}
void ewaldself(int ii){
	int i,j,k,l,m,n,kk,a,b,xbin,ybin,zbin,startij,startji,startii,start,kx,ky,kz;
	double dx, dy, dz, r, rr, C1, C2, C3, C4, coeff1, coeff2, term1, term2, coskr,incr,rkx,rky,rkz,rkk;
	double Dijtemp[3][3];
	double M1[3][3];
	selfcount++;
	for(i=0;i<Nb;i++){
		for(j=i+1;j<Nb;j++){
			for(a=0;a<3;a++){
				for(b=0;b<3;b++){
					Dijtemp[a][b] = 0.0;
				}
			}
			dx = rx[ii*Nb+i] - rx[ii*Nb+j]; dy = ry[ii*Nb+i] - ry[ii*Nb+j]; dz = rz[ii*Nb+i] - rz[ii*Nb+j];
			dx -= box_length*round(dx/box_length); dy -= box_length*round(dy/box_length); dz -= box_length*round(dz/box_length);
			rr = dx*dx+dy*dy+dz*dz; r = sqrt(rr);
			// Real Part
			if(r>2.0){
				C1 = 0.75/r+0.5/(r*rr);
				C2 = 4.0*alpha7*rr*rr + 3.0*alpha3*rr - 20.0*alpha5*rr - 4.5*alpha + 14.0*alpha3 + alpha/rr;
				C3 = 0.75/r-1.5/(r*rr);
				C4 = -4.0*alpha7*rr*rr - 3.0*alpha3*rr + 16.0*alpha5*rr + 1.5*alpha - 2.0*alpha3 - 3.0*alpha/rr;
			}
			else{
				C1 = 1.0 - 0.28125*r;
				C2 = 64.0*alpha7 + 12.0*alpha3 - 80.0*alpha5 - 4.5*alpha + 14.0*alpha3 + 0.25*alpha;
				C3 = 0.09375*r;
				C4 = -64.0*alpha7 - 12.0*alpha3 + 64.0*alpha5 + 1.5*alpha - 2.0*alpha3 - 0.75*alpha;
			}
			term1 = C1*erfc(alpha*r) + C2*exp(-alpha2*rr)/rtpi;
			term2 = C3*erfc(alpha*r) + C4*exp(-alpha2*rr)/rtpi;
			M1[0][0] = term1 + term2*dx*dx/rr;
			M1[0][1] = term2*dx*dy/rr;
			M1[0][2] = term2*dx*dz/rr;
			M1[1][0] = M1[0][1];
			M1[1][1] = term1 + term2*dy*dy/rr;
			M1[1][2] = term2*dy*dz/rr;
			M1[2][0] = M1[0][2];
			M1[2][1] = M1[1][2];
			M1[2][2] = term1 + term2*dz*dz/rr;
			Dijtemp[0][0] += M1[0][0];
			Dijtemp[0][1] += M1[0][1];
			Dijtemp[0][2] += M1[0][2];
			Dijtemp[1][0] += M1[1][0];
			Dijtemp[1][1] += M1[1][1];
			Dijtemp[1][2] += M1[1][2];
			Dijtemp[2][0] += M1[2][0];
			Dijtemp[2][1] += M1[2][1];
			Dijtemp[2][2] += M1[2][2];
			// Loop over wavevectors for Drecip
			for(kx = -kmax; kx < kmax+1; kx++){
				rkx = kx*kcoeff;
				for(ky = -kmax; ky < kmax+1; ky++){
					rky = ky*kcoeff;
					for(kz = -kmax; kz < kmax+1; kz++){
						rkz = kz*kcoeff;
						kk = kx*kx + ky*ky + kz*kz;
						if(kk != 0){
							coskr = cos(rkx*dx + rky*dy + rkz*dz);
							start = 9*((kx+kmax)*numk*numk+(ky+kmax)*numk+(kz+kmax));
							Dijtemp[0][0] += M2[start]*coskr;
							Dijtemp[0][1] += M2[start+1]*coskr;
							Dijtemp[0][2] += M2[start+2]*coskr;
							Dijtemp[1][0] += M2[start+3]*coskr;
							Dijtemp[1][1] += M2[start+4]*coskr;
							Dijtemp[1][2] += M2[start+5]*coskr;
							Dijtemp[2][0] += M2[start+6]*coskr;
							Dijtemp[2][1] += M2[start+7]*coskr;
							Dijtemp[2][2] += M2[start+8]*coskr;
						}
					}
				}
			}
			startij = 9*Nb*i+9*j;
			startji = 9*Nb*j+9*i;
			for(a=0;a<3;a++){
				for(b=0;b<3;b++){
					Dijselfavg[startij+3*a+b] += Dijtemp[a][b];
					Dijselfavg[startji+3*a+b] += Dijtemp[a][b];
					incr = Dijtemp[a][b]/Dtr;
					eps += 2.0*incr;
					eps2 += 2.0*incr*incr;
					rowsum[3*Nb*ii+3*i+a] += incr*incr;
					rowsum[3*Nb*ii+3*j+a] += incr*incr;
				}
			}
		}
	}
	for(i=0;i<Nb;i++){
		startii = 9*Nb*i + 9*i;
		for(kx = -kmax; kx < kmax+1; kx++){
			for(ky = -kmax; ky < kmax+1; ky++){
				for(kz = -kmax; kz < kmax +1; kz++){
					kk = kx*kx + ky*ky + kz*kz;
					if(kk != 0){
						start = 9*((kx+kmax)*numk*numk+(ky+kmax)*numk+(kz+kmax));
						Dijselfavg[startii] += M2[start];
						Dijselfavg[startii+4] += M2[start+4];
						Dijselfavg[startii+8] += M2[start+8];
					}
				}
			}
		}
		Dijselfavg[startii] += selfconst;
		Dijselfavg[startii+4] += selfconst;
		Dijselfavg[startii+8] += selfconst;
	}
}
void gwdecompmod(){
	double tempd;
	int i,j,k,l,m;
	int ROW = 3*N;
	// eps /= (3*N);
	// eps2 /= (3*N);
	eps /= (ROW*ROW);
	eps2 = eps*eps;
    betaij = (1.0-sqrt(1.0-(ROW*eps2-ROW*eps)))/(ROW*eps2-ROW*eps);
	// printf("%lu %lf %lf %lf\n",t,betaij,eps,eps2);
	if(betaij<1.0){
		// printf("opt1\n");
		// printf("%lu %lf %lf %lf\n",t,betaij,eps,eps2);
		for(i=0;i<3*N;i++){
			C[i] = sqrt(1.0/(1.0+rowsum[i]*betaij*betaij));
		}
		brun += betaij;
		for(i=0;i<3*N;i++){
			Crun[i] += C[i];
		}
		Bfile = fopen(betavt,"a");
		fprintf(Bfile,"%lu 1 %lf\n",t,betaij);
		fclose(Bfile);
		// Cfile = fopen(cvt,"a");
		// fprintf(Cfile,"%lu 1\n",t);
		// for(i=0;i<N;i++){
		// 	fprintf(Cfile,"%d %lf %lf %lf\n",i,C[3*i],C[3*i+1],C[3*i+2]);
		// }
		// fclose(Cfile);
		gwcount++;
	}
	else{
		// printf("opt2\n");
		Bfile = fopen(betavt,"a");
		fprintf(Bfile,"%lu 0 \n",t);
		fclose(Bfile);
		// Cfile = fopen(cvt,"a");
		// fprintf(Cfile,"%lu 0\n",t);
		// fclose(Cfile);
	}
	eps = 0.0; eps2 = 0.0;
	for(i=0;i<3*N;i++){
		rowsum[i] = 0.0;
	}
}
void getNoise(){
	int i;
	for(i = 0; i<N*3; ++i){
		R[i] = gasdev(idum);
	}
}
void updateChains(){
	int i,j,k,l,m,n,a,b,start,startavg,startij,startji,bin,ii,jj;
	double dx,dy,dz,start_time,betaii;
	int offset,offsum,offiter;

	if(itercount==0){
		for(i=0;i<N;i++){
			dx = 0.0; dy = 0.0; dz = 0.0;
			dx = dt*fx[i]+pc*R[3*i];
			dy = dt*fy[i]+pc*R[3*i+1];
			dz = dt*fz[i]+pc*R[3*i+2];
			rx[i] += dx; ry[i] += dy; rz[i] += dz;
			rx[i] -= box_length*round(rx[i]/box_length);
			ry[i] -= box_length*round(ry[i]/box_length);
			rz[i] -= box_length*round(rz[i]/box_length);
		}
	}
	else{
		// -------------------- 1) Find bin locations but do not construct square MM ------------------
		double *dxu = calloc(N, sizeof(double));
		double *dyu = calloc(N, sizeof(double));
		double *dzu = calloc(N, sizeof(double));
		if(t%10==0){
			#pragma omp parallel for schedule(static,1) private(j,bin)
				for(i=0;i<N;i++){
					// for(j=0;j<((int)(i/Nb)*Nb);j++){
					for(j=(((int)(i/Nb)+1)*Nb);j<N;j++){
						// bin = binmod(i,j);
						bin = binmod(j,i);
						start_ij[i*N+j] = 9*bin;
					}
				}
		}
		// start_time = omp_get_wtime();
		#pragma omp parallel
		{
			int i,j,m,n,start,k;
			double betaii;
			double *dxn = calloc(N, sizeof(double));
			double *dyn = calloc(N, sizeof(double));
			double *dzn = calloc(N, sizeof(double));
			// int id = omp_get_thread_num();
			#pragma omp for schedule(static,1)
				for(i=0;i<N;i++){
					m = i*3;
					// for(j=0;j<((int)(i/Nb)*Nb);j++){
					for(j=(((int)(i/Nb)+1)*Nb);j<N;j++){
						n = j*3;
						start = start_ij[i*N+j];
						dxn[i] += dt*(Dij[start]*(fx[j]+p*Cavg[m]*bavg*R[n])+Dij[start+1]*(fy[j]+p*Cavg[m]*bavg*R[n+1])+Dij[start+2]*(fz[j]+p*Cavg[m]*bavg*R[n+2]));
						dyn[i] += dt*(Dij[start+3]*(fx[j]+p*Cavg[m+1]*bavg*R[n])+Dij[start+4]*(fy[j]+p*Cavg[m+1]*bavg*R[n+1])+Dij[start+5]*(fz[j]+p*Cavg[m+1]*bavg*R[n+2]));
						dzn[i] += dt*(Dij[start+6]*(fx[j]+p*Cavg[m+2]*bavg*R[n])+Dij[start+7]*(fy[j]+p*Cavg[m+2]*bavg*R[n+1])+Dij[start+8]*(fz[j]+p*Cavg[m+2]*bavg*R[n+2]));
						dxn[j] += dt*(Dij[start]*(fx[i]+p*Cavg[n]*bavg*R[m])+Dij[start+1]*(fy[i]+p*Cavg[n]*bavg*R[m+1])+Dij[start+2]*(fz[i]+p*Cavg[n]*bavg*R[m+2]));
						dyn[j] += dt*(Dij[start+3]*(fx[i]+p*Cavg[n+1]*bavg*R[m])+Dij[start+4]*(fy[i]+p*Cavg[n+1]*bavg*R[m+1])+Dij[start+5]*(fz[i]+p*Cavg[n+1]*bavg*R[m+2]));
						dzn[j] += dt*(Dij[start+6]*(fx[i]+p*Cavg[n+2]*bavg*R[m])+Dij[start+7]*(fy[i]+p*Cavg[n+2]*bavg*R[m+1])+Dij[start+8]*(fz[i]+p*Cavg[n+2]*bavg*R[m+2]));
					}
					// for(j=((int)(i/Nb)*Nb);j<((int)(i/Nb)+1)*Nb;j++){
					for(j=i+1;j<((int)(i/Nb)+1)*Nb;j++){
						n = j*3;
						// if(i!=j){betaii = bavg;} // i!=j -> Bii = Bij
						// else{betaii = 1.0;} // i=j -> Bii = 1.0
						start = 9*Nb*(i%Nb)+9*(j%Nb);
						// dxn[i] += dt*(Dijself[start+0]*(fx[j]+p*Cavg[m]*betaii*R[n])+Dijself[start+1]*(fy[j]+p*Cavg[m]*bavg*R[n+1])+Dijself[start+2]*(fz[j]+p*Cavg[m]*bavg*R[n+2]));
						// dyn[i] += dt*(Dijself[start+3]*(fx[j]+p*Cavg[m+1]*bavg*R[n])+Dijself[start+4]*(fy[j]+p*Cavg[m+1]*betaii*R[n+1])+Dijself[start+5]*(fz[j]+p*Cavg[m+1]*bavg*R[n+2]));
						// dzn[i] += dt*(Dijself[start+6]*(fx[j]+p*Cavg[m+2]*bavg*R[n])+Dijself[start+7]*(fy[j]+p*Cavg[m+2]*bavg*R[n+1])+Dijself[start+8]*(fz[j]+p*Cavg[m+2]*betaii*R[n+2]));
						// dxu[i] += dt*(Dijself[start+0]*(fx[j]+p*Cavg[m]*betaii*R[n])+Dijself[start+1]*(fy[j]+p*Cavg[m]*bavg*R[n+1])+Dijself[start+2]*(fz[j]+p*Cavg[m]*bavg*R[n+2]));
						// dyu[i] += dt*(Dijself[start+3]*(fx[j]+p*Cavg[m+1]*bavg*R[n])+Dijself[start+4]*(fy[j]+p*Cavg[m+1]*betaii*R[n+1])+Dijself[start+5]*(fz[j]+p*Cavg[m+1]*bavg*R[n+2]));
						// dzu[i] += dt*(Dijself[start+6]*(fx[j]+p*Cavg[m+2]*bavg*R[n])+Dijself[start+7]*(fy[j]+p*Cavg[m+2]*bavg*R[n+1])+Dijself[start+8]*(fz[j]+p*Cavg[m+2]*betaii*R[n+2]));
						dxn[i] += dt*(Dijself[start+0]*(fx[j]+p*Cavg[m]*bavg*R[n])+Dijself[start+1]*(fy[j]+p*Cavg[m]*bavg*R[n+1])+Dijself[start+2]*(fz[j]+p*Cavg[m]*bavg*R[n+2]));
						dyn[i] += dt*(Dijself[start+3]*(fx[j]+p*Cavg[m+1]*bavg*R[n])+Dijself[start+4]*(fy[j]+p*Cavg[m+1]*bavg*R[n+1])+Dijself[start+5]*(fz[j]+p*Cavg[m+1]*bavg*R[n+2]));
						dzn[i] += dt*(Dijself[start+6]*(fx[j]+p*Cavg[m+2]*bavg*R[n])+Dijself[start+7]*(fy[j]+p*Cavg[m+2]*bavg*R[n+1])+Dijself[start+8]*(fz[j]+p*Cavg[m+2]*bavg*R[n+2]));
						dxn[j] += dt*(Dijself[start+0]*(fx[i]+p*Cavg[n]*bavg*R[m])+Dijself[start+1]*(fy[i]+p*Cavg[n]*bavg*R[m+1])+Dijself[start+2]*(fz[i]+p*Cavg[n]*bavg*R[m+2]));
						dyn[j] += dt*(Dijself[start+3]*(fx[i]+p*Cavg[n+1]*bavg*R[m])+Dijself[start+4]*(fy[i]+p*Cavg[n+1]*bavg*R[m+1])+Dijself[start+5]*(fz[i]+p*Cavg[n+1]*bavg*R[m+2]));
						dzn[j] += dt*(Dijself[start+6]*(fx[i]+p*Cavg[n+2]*bavg*R[m])+Dijself[start+7]*(fy[i]+p*Cavg[n+2]*bavg*R[m+1])+Dijself[start+8]*(fz[i]+p*Cavg[n+2]*bavg*R[m+2]));
					}
					dxn[i] += dt*Dtr*(fx[i]+p*Cavg[m]*R[m]);
					dyn[i] += dt*Dtr*(fy[i]+p*Cavg[m+1]*R[m+1]);
					dzn[i] += dt*Dtr*(fz[i]+p*Cavg[m+2]*R[m+2]);
					// start = 9*Nb*(i%Nb)+9*(i%Nb);
					// dxn[i] += dt*Dijself[start]*(fx[i]+p*Cavg[m]*R[m]);
					// dyn[i] += dt*Dijself[start+4]*(fy[i]+p*Cavg[m+1]*R[m+1]);
					// dzn[i] += dt*Dijself[start+8]*(fz[i]+p*Cavg[m+2]*R[m+2]);
				}
			#pragma omp critical
			{
				for(k=0;k<N;k++){
					dxu[k] += dxn[k];
					dyu[k] += dyn[k];
					dzu[k] += dzn[k];
				}
				free(dxn); free(dyn); free(dzn);
			}
		}
		for(i=0;i<N;i++){
			// if(fabs(dxu[i]) > 1.0 || fabs(dyu[i]) > 1.0 || fabs(dzu[i]) > 1.0){
			// 	printf("%lu %d %lf %lf %lf %lf %lf %lf\n",t,i,dxu[i],dyu[i],dzu[i],fx[i],fy[i],fz[i]);
			// 	// exit(1);
			// }
			rx[i] += dxu[i]; ry[i] += dyu[i]; rz[i] += dzu[i];
			// Put beads back in the box
			rx[i] -= box_length*round(rx[i]/box_length);
			ry[i] -= box_length*round(ry[i]/box_length);
			rz[i] -= box_length*round(rz[i]/box_length);
		}


		free(dxu); free(dyu); free(dzu);
		// exit(1);
		// ljtimetest += omp_get_wtime()-start_time;
	}
}
void updateFD(){
	int i,j,k;


	// double dx,dy,dz;
	for(i=0;i<N;i++){
		// dx = 0.0; dy = 0.0; dz = 0.0;
		rx[i] += dt*fx[i]+pc*R[3*i];
		ry[i] += dt*fy[i]+pc*R[3*i+1];
		rz[i] += dt*fz[i]+pc*R[3*i+2];

		// rx[i] += dx; ry[i] += dy; rz[i] += dz;
		rx[i] -= box_length*round(rx[i]/box_length);
		ry[i] -= box_length*round(ry[i]/box_length);
		rz[i] -= box_length*round(rz[i]/box_length);
	}



}
void updateBins(){
	int i,j,k,l,m,n,a,b,start,startavg,startij,startji,bin,ii,jj;
	double dx,dy,dz,start_time,betaii;
	int offset,offsum,offiter;
	// -------------------- 1) Find bin locations but do not construct square MM ------------------
	double *dxu = calloc(N, sizeof(double));
	double *dyu = calloc(N, sizeof(double));
	double *dzu = calloc(N, sizeof(double));
	if(t%10==0){
		#pragma omp parallel for schedule(static,1) private(j,bin)
			for(i=0;i<N;i++){
				for(j=i;j<N;j++){
				// for(j=0;j<i;j++){
				// for(j=(((int)(i/Nb)+1)*Nb);j<N;j++){
					// bin = binmod(i,j);
					bin = binmod(j,i);
					start_ij[i*N+j] = 9*bin;
				}
			}
		// printf("%lu %ld\n",t,start_ij[100*N]);
		// for(i=0;i<N;i++){
		// 	printf("999 %d %d\n",i,start_ij[(N-1)*N+i]);
		// }
		// exit(1);
	}
	start_time = omp_get_wtime();
	#pragma omp parallel
	{
		int i,j,m,n,start,k;
		double betaii;
		double *dxn = calloc(N, sizeof(double));
		double *dyn = calloc(N, sizeof(double));
		double *dzn = calloc(N, sizeof(double));
		// int id = omp_get_thread_num();
		// printf("%lf %lf %lf\n",bavg,Cavg[0],Dij[0]);
		#pragma omp for schedule(static,1)
			for(i=0;i<N;i++){
				m = i*3;
				// printf("%d %d %d\n",id,i,(int)(i/Nb)*Nb);
				// for(j=0;j<i;j++){
				for(j=i;j<N;j++){
					n = j*3;
					start = start_ij[i*N+j];
					dxn[i] += dt*(Dij[start]*(fx[j]+p*Cavg[m]*bavg*R[n])+Dij[start+1]*(fy[j]+p*Cavg[m]*bavg*R[n+1])+Dij[start+2]*(fz[j]+p*Cavg[m]*bavg*R[n+2]));
					dyn[i] += dt*(Dij[start+3]*(fx[j]+p*Cavg[m+1]*bavg*R[n])+Dij[start+4]*(fy[j]+p*Cavg[m+1]*bavg*R[n+1])+Dij[start+5]*(fz[j]+p*Cavg[m+1]*bavg*R[n+2]));
					dzn[i] += dt*(Dij[start+6]*(fx[j]+p*Cavg[m+2]*bavg*R[n])+Dij[start+7]*(fy[j]+p*Cavg[m+2]*bavg*R[n+1])+Dij[start+8]*(fz[j]+p*Cavg[m+2]*bavg*R[n+2]));
					dxn[j] += dt*(Dij[start]*(fx[i]+p*Cavg[n]*bavg*R[m])+Dij[start+1]*(fy[i]+p*Cavg[n]*bavg*R[m+1])+Dij[start+2]*(fz[i]+p*Cavg[n]*bavg*R[m+2]));
					dyn[j] += dt*(Dij[start+3]*(fx[i]+p*Cavg[n+1]*bavg*R[m])+Dij[start+4]*(fy[i]+p*Cavg[n+1]*bavg*R[m+1])+Dij[start+5]*(fz[i]+p*Cavg[n+1]*bavg*R[m+2]));
					dzn[j] += dt*(Dij[start+6]*(fx[i]+p*Cavg[n+2]*bavg*R[m])+Dij[start+7]*(fy[i]+p*Cavg[n+2]*bavg*R[m+1])+Dij[start+8]*(fz[i]+p*Cavg[n+2]*bavg*R[m+2]));
					// dxu[i] += dt*(Dij[start]*(fx[j]+p*Cavg[m]*bavg*R[n])+Dij[start+1]*(fy[j]+p*Cavg[m]*bavg*R[n+1])+Dij[start+2]*(fz[j]+p*Cavg[m]*bavg*R[n+2]));
					// dyu[i] += dt*(Dij[start+3]*(fx[j]+p*Cavg[m+1]*bavg*R[n])+Dij[start+4]*(fy[j]+p*Cavg[m+1]*bavg*R[n+1])+Dij[start+5]*(fz[j]+p*Cavg[m+1]*bavg*R[n+2]));
					// dzu[i] += dt*(Dij[start+6]*(fx[j]+p*Cavg[m+2]*bavg*R[n])+Dij[start+7]*(fy[j]+p*Cavg[m+2]*bavg*R[n+1])+Dij[start+8]*(fz[j]+p*Cavg[m+2]*bavg*R[n+2]));
					// dxu[j] += dt*(Dij[start]*(fx[i]+p*Cavg[n]*bavg*R[m])+Dij[start+1]*(fy[i]+p*Cavg[n]*bavg*R[m+1])+Dij[start+2]*(fz[i]+p*Cavg[n]*bavg*R[m+2]));
					// dyu[j] += dt*(Dij[start+3]*(fx[i]+p*Cavg[n+1]*bavg*R[m])+Dij[start+4]*(fy[i]+p*Cavg[n+1]*bavg*R[m+1])+Dij[start+5]*(fz[i]+p*Cavg[n+1]*bavg*R[m+2]));
					// dzu[j] += dt*(Dij[start+6]*(fx[i]+p*Cavg[n+2]*bavg*R[m])+Dij[start+7]*(fy[i]+p*Cavg[n+2]*bavg*R[m+1])+Dij[start+8]*(fz[i]+p*Cavg[n+2]*bavg*R[m+2]));
				}
				dxn[i] += dt*Dtr*(fx[i]+p*Cavg[m]*R[m]);
				dyn[i] += dt*Dtr*(fy[i]+p*Cavg[m+1]*R[m+1]);
				dzn[i] += dt*Dtr*(fz[i]+p*Cavg[m+2]*R[m+2]);
				// start = 9*Nb*(i%Nb)+9*(i%Nb);
				// dxn[i] += dt*(Dijself[start+0]*(fx[i]+p*Cavg[m]*R[m])+Dijself[start+1]*(fy[i]+p*Cavg[m]*bavg*R[m+1])+Dijself[start+2]*(fz[i]+p*Cavg[m]*bavg*R[m+2]));
				// dyn[i] += dt*(Dijself[start+3]*(fx[i]+p*Cavg[m+1]*bavg*R[m])+Dijself[start+4]*(fy[i]+p*Cavg[m+1]*R[m+1])+Dijself[start+5]*(fz[i]+p*Cavg[m+1]*bavg*R[m+2]));
				// dzn[i] += dt*(Dijself[start+6]*(fx[i]+p*Cavg[m+2]*bavg*R[m])+Dijself[start+7]*(fy[i]+p*Cavg[m+2]*bavg*R[m+1])+Dijself[start+8]*(fz[i]+p*Cavg[m+2]*R[m+2]));
			}
		#pragma omp critical
		{
			for(k=0;k<N;k++){
				dxu[k] += dxn[k];
				dyu[k] += dyn[k];
				dzu[k] += dzn[k];
			}
			free(dxn); free(dyn); free(dzn);
		}
	}
	// printf("0 %lf %lf %lf\n",i,dxu[0],dyu[0],dzu[0]);
	// exit(1);
	for(i=0;i<N;i++){
		// printf("%d %lf %lf %lf\n",i,dxu[i],dyu[i],dzu[i]);
		// if(fabs(dxu[i]) > 1.0 || fabs(dyu[i]) > 1.0 || fabs(dzu[i]) > 1.0){
			// printf("%lu %d %lf %lf %lf %lf %lf %lf\n",t,i,dxu[i],dyu[i],dzu[i],fx[i],fy[i],fz[i]);
			// // exit(1);
		// }
		rx[i] += dxu[i]; ry[i] += dyu[i]; rz[i] += dzu[i];
		// Put beads back in the box
		rx[i] -= box_length*round(rx[i]/box_length);
		ry[i] -= box_length*round(ry[i]/box_length);
		rz[i] -= box_length*round(rz[i]/box_length);
		// Update displacements from bead positions when neighbor lists last updated
		// dr[i] += sqrt(dxu[i]*dxu[i]+dyu[i]*dyu[i]+dzu[i]*dzu[i]);
	}
	free(dxu); free(dyu); free(dzu);
	// exit(1);
	ljtimetest += omp_get_wtime()-start_time;
}
void checkVerlet(){
	// Update when dr > 1/2(rv-rc)
	int i,j;
	double drcheck,dx,dy,dz;
	for(i=0;i<N;i++){
		dx = rx[i] - xo[i]; dy = ry[i] - yo[i]; dz = rz[i] - zo[i];
		dx -= box_length*round(dx/box_length); dy -= box_length*round(dy/box_length); dz -= box_length*round(dz/box_length);
		drcheck = sqrt(dx*dx + dy*dy + dz*dz);
		if(drcheck > rnew){ // Condition met, need to update lists
			verletlist();
			break;
		}
	}
}
void printTrajectory(){
	int i,j,ind;
	// printf("%lu %lu\n",t,t%1000);
	// if(t%1000==0){
		// if(restart==0 || restart==2 || (restart==1 && t!=tstart)){
			xyzfile = fopen(xyz, "a");
			//fprintf(xyzfile,"%d\n",N);
			// fprintf(outputfile,"%lf %lf %lf %lf %lf %lf %lf\n",L1[0],L1[1],L2[0],L2[1],point[1][0],point[1][1],-lz/2.0);
			//fprintf(xyzfile,"Lattice=\"%lf %lf %lf %lf %lf %lf %lf %lf %lf\" Properties=\"\" Time=%lu\n",Lbar[0],Lbar[1],Lbar[2],Lbar[3],Lbar[4],Lbar[5],Lbar[6],Lbar[7],Lbar[8],t);
			fprintf(xyzfile, "%d\n%lu\n", N, t);
			for(i = 0; i<Nc; ++i){
				for(j=0;j<Nb;j++){
						ind = Nb*i +j;
						if(Charge[ind] == 0){
							fprintf(xyzfile, "A %lf %lf %lf\n", rx[Nb*i+j], ry[Nb*i+j], rz[Nb*i+j]);
						}
						else if(Charge[ind] == 1){
							fprintf(xyzfile, "B %lf %lf %lf\n", rx[Nb*i+j], ry[Nb*i+j], rz[Nb*i+j]);
						}
				}
			}
			fclose(xyzfile);
		// }
	// }
}
void printTEA(){
	int i;
	DCfile = fopen(decomposition,"w");
	fprintf(DCfile,"%lu %.12e\n",t,brun/gwcount);
	for(i=0;i<N;i++){
		// fprintf(DCfile,"%d %lf %lf %lf\n",i,Cavg[3*i],Cavg[3*i+1],Cavg[3*i+2]);
		fprintf(DCfile,"%d %.12e %.12e %.12e\n",i,Crun[3*i]/gwcount,Crun[3*i+1]/gwcount,Crun[3*i+2]/gwcount);
	}
	fclose(DCfile);
}
void printMM(){
	int i,j,k,l,m,n,a,b,start;
	// sprintf(matrix, "mm/M%d_%d_%.3lf_%.3lf_%d_%d.txt", Nb, Nc, c_norm, epsilon, trace, itercount);
	MMfile = fopen(matrix, "w");
	fprintf(MMfile,"%lu %ld\n",t,selfcount);
	for(i=0;i<Nb;i++){
		for(j=0;j<Nb;j++){
			start = 9*Nb*i + 9*j;
			fprintf(MMfile,"%d %d\n",i,j);
			for(a=0;a<3;a++){
				for(b=0;b<3;b++){
					// fprintf(MMfile,"%lf ",Dijself[start+3*a+b]);
					fprintf(MMfile,"%.12e ",Dijselfavg[start+3*a+b]/selfcount);
				}
				fprintf(MMfile,"\n");
			}
		}
	}
	for(i=0;i<gd1;i++){
		for(j=0;j<gd1;j++){
			for(k=0;k<gd1;k++){
				start = i*gd1*gd1+j*gd1+k;
				fprintf(MMfile,"%d %d %d %ld\n",i,j,k,count[start]);
				if(count[start]>0){
					for(l=0;l<3;l++){
						for(m=0;m<3;m++){
							// fprintf(MMfile,"%lf ",Dij[9*start+3*l+m]);
							fprintf(MMfile,"%.12e ",Dijavg[9*start+3*l+m]/count[start]);
						}
						fprintf(MMfile,"\n");
					}
				}
				else{
					for(l=0;l<3;l++){
						for(m=0;m<3;m++){
							// fprintf(MMfile,"%lf ",Dij[9*start+3*l+m]);
							fprintf(MMfile,"%.12e ",0.0);
						}
						fprintf(MMfile,"\n");
					}
				}
			}
		}
	}
	for(i=0;i<gd2;i++){
		for(j=0;j<gd2;j++){
			for(k=0;k<gd2;k++){
				start = gn1 + i*gd2*gd2 + j*gd2 + k;
				fprintf(MMfile,"%d %d %d %ld\n",i,j,k,count[start]);
				if(count[start]>0){
					for(l=0;l<3;l++){
						for(m=0;m<3;m++){
							// fprintf(MMfile,"%lf ",Dij[9*start+3*l+m]);
							fprintf(MMfile,"%.12e ",Dijavg[9*start+3*l+m]/count[start]);
						}
						fprintf(MMfile,"\n");
					}
				}
				else{
					for(l=0;l<3;l++){
						for(m=0;m<3;m++){
							// fprintf(MMfile,"%lf ",Dij[9*start+3*l+m]);
							fprintf(MMfile,"%.12e ",0.0);
						}
						fprintf(MMfile,"\n");
					}
				}
			}
		}
	}
	fclose(MMfile);
	// sprintf(matrix, "mm/M2%d_%d_%.3lf_%.3lf_%d_%d.txt", Nb, Nc, c_norm, epsilon, trace, itercount);
	// MMfile = fopen(matrix,"w");
	// // fprintf(MMfile,"%lu %ld\n",t,selfcount);
	// for(i=0;i<9*Nb*Nb;i++){
	// 	fprintf(MMfile,"%.12e\n",Dijselfavg[i]/selfcount);
	// }
	// for(i=0;i<gnt;i++){
	// 	if(count[i]>0){
	// 		for(j=0;j<9;j++){
	// 			fprintf(MMfile,"%.12e\n",Dijavg[9*i+j]/count[i]);
	// 		}
	// 	}
	// 	else{
	// 		for(j=0;j<9;j++){
	// 			fprintf(MMfile,"%.12e\n",0.0);
	// 		}
	// 	}
	// }
}
void resetAverage(){
	int a,i,j,k,start;
	bavg = brun/gwcount;
	brun = 0.0;
	for(i=0;i<3*N;i++){
		Cavg[i] = Crun[i]/gwcount;
		Crun[i] = 0.0;
	}
	gwcount = 0;
	// printf("Self Avg");
	for(i=0;i<9*Nb*Nb;i++){
		Dijself[i] = Dijselfavg[i]/selfcount;
		// printf("%.12e\n",Dijself[i]);
	}
	for(i=0;i<gnt;i++){
		for(j=0;j<9;j++){
			if(count[i]>0){
				Dij[9*i+j] = Dijavg[9*i+j]/count[i];
			}
		}
	}
	// printf("Grid Avg");
	// for(i=0;i<9*gnt;i++){
	// 	printf("%.12e\n",Dij[i]);
	// }
	selfcount = 0.0;
	for(i=0;i<9*Nb*Nb;i++){
		Dijselfavg[i] = 0.0;
	}
	for(i=0;i<9*gnt;i++){
		Dijavg[i] = 0.0;
	}
	for(i=0;i<gnt;i++){
		count[i] = 0;
	}
}
void calcRgCM(){
	int i,j,k,ind2,sample_number;
	double dx,dy,dz,rgxt,rgyt,rgzt,COMx,COMy,COMz,dxo,dyo,dzo,reext,reeyt,reezt;
	for(i=0;i<Nc;i++){
		px[Nb*i] = rx[Nb*i]; py[Nb*i] = ry[Nb*i]; pz[Nb*i] = rz[Nb*i];
		// px[Nb*i] = 0.0; py[Nb*i] = 0.0; pz[Nb*i] = 0.0;
		for(j=0;j<Nb-1;j++){
			dx = rx[Nb*i+j+1] - rx[Nb*i+j]; dy = ry[Nb*i+j+1] - ry[Nb*i+j]; dz = rz[Nb*i+j+1] - rz[Nb*i+j];
			dx -= box_length*round(dx/box_length); dy -= box_length*round(dy/box_length); dz -= box_length*round(dz/box_length);
			px[Nb*i+j+1] = px[Nb*i+j] + dx; py[Nb*i+j+1] = py[Nb*i+j] + dy; pz[Nb*i+j+1] = pz[Nb*i+j] + dz;
		}
	}
	// rgtfile = fopen(rgt,"a");
	// fprintf(rgtfile,"%lu\n",t);
	// sample_number = (t%printperiod)*dt;
	// tprint[sample_number] = t;
	// printf("%d %lu\n",sample_number,tprint[sample_number]);
	for(i=0;i<Nc;i++){
		rgxt = 0.0; rgyt = 0.0; rgzt = 0.0;
		for(j=0;j<Nb;j++){
			for(k=j;k<Nb;k++){
				dx = px[Nb*i+j] - px[Nb*i+k]; dy = py[Nb*i+j] - py[Nb*i+k]; dz = pz[Nb*i+j] - pz[Nb*i+k];
				rgxt += dx*dx; rgyt += dy*dy; rgzt += dz*dz;
			}
		}
		// rgx[i] = sqrt(rgxt/(Nb*Nb)); rgy[i] = sqrt(rgyt/(Nb*Nb)); rgz[i] = sqrt(rgzt/(Nb*Nb));
		rg[i] = sqrt((rgxt+rgyt+rgzt)/(Nb*Nb));
		// fprintf(rgtfile,"%d %lf %lf %lf %lf\n",i,sqrt(rgx/(Nb*Nb)),sqrt(rgy/(Nb*Nb)),sqrt(rgz/(Nb*Nb)),sqrt((rgx+rgy+rgz)/(Nb*Nb)));
	}
	// fclose(rgtfile);
	// cmtfile = fopen(cmt,"a");
	// fprintf(cmtfile,"%lu\n",t);
	for(i=0;i<Nc;i++){
		dxo = 0.0; dyo = 0.0; dzo = 0.0;
		COMx = 0.0; COMy = 0.0; COMz = 0.0;
		for(j=0;j<Nb-1;j++){
			dx = rx[Nb*i+j+1] - rx[Nb*i+j]; dy = ry[Nb*i+j+1] - ry[Nb*i+j]; dz = rz[Nb*i+j+1] - rz[Nb*i+j];
			dx -= box_length*round(dx/box_length); dy -= box_length*round(dy/box_length); dz -= box_length*round(dz/box_length);
			dxo += dx; dyo += dy; dzo += dz;
			// printf("%d %d %lf %lf %lf\n",i,j,dxo,dyo,dzo);
			COMx += dxo; COMy += dyo; COMz += dzo;
		}
		COMx /= Nb; COMy /= Nb; COMz /= Nb;
		COMx += rx[Nb*i]; COMy += ry[Nb*i]; COMz += rz[Nb*i];
		COMx -= box_length*round(COMx/box_length); COMy -= box_length*round(COMy/box_length); COMz -= box_length*round(COMz/box_length);
		comx[i] = COMx; comy[i] = COMy; comz[i] = COMz;
		// fprintf(cmtfile,"%d %lf %lf %lf\n",i,COMx,COMy,COMz);
	}
	for(i=0;i<Nc;i++){
		// reext = 0.0; reeyt = 0.0; reezt = 0.0;
		// for(j=0;j<Nb-1;j++){
		// 	dx = rx[Nb*i+j+1] - rx[Nb*i+j]; dy = ry[Nb*i+j+1] - ry[Nb*i+j]; dz = rz[Nb*i+j+1] - rz[Nb*i+j];
		// 	dx -= box_length*round(dx/box_length); dy -= box_length*round(dy/box_length); dz -= box_length*round(dz/box_length);
		// 	reext += dx; reeyt += dy; reezt += dz;
		// }
		reex[i] = px[Nb*(i+1)-1] - px[Nb*i]; reey[i] = py[Nb*(i+1)-1] - py[Nb*i]; reez[i] = pz[Nb*(i+1)-1] - pz[Nb*i];
	}
	// fclose(cmtfile);
	rgtfile = fopen(rgt,"a");
	cmtfile = fopen(cmt,"a");
	reetfile = fopen(reet,"a");
	//fprintf(rgtfile,"%lu\n",t);
	//fprintf(cmtfile,"%lu\n",t);
	//fprintf(reetfile,"%lu\n",t);
	for(i=0;i<Nc;i++){
		fprintf(rgtfile,"%lu %d %lf\n",t,i,rg[i]);
		fprintf(cmtfile,"%lu %d %lf %lf %lf\n",t,i,comx[i],comy[i],comz[i]);
		fprintf(reetfile,"%lu %d %lf %lf %lf\n",t,i,reex[i],reey[i],reez[i]);
	}
	fclose(rgtfile);
	fclose(cmtfile);
	fclose(reetfile);
}
void calcMSD()
{
  int nsamp,nq,ind1,ind2,count,m,l, wmax,Tmax,q,k,ctemp,n_it , n_nc,nc,ttemp1, i, j, jcount, nb1, n_nb1, nb2, n_nb2,nedot,n_nedot,n_nit,nit, dint,dint2;
  double dcomx,dcomy,dcomz,dx,dy,dz,rxtemp,rytemp,rztemp,MSDtemp,avgmsd;

  Tmax = tmax/printprops;
	wmax = Tmax+1;
	nsamp = Ncharges*wmax;


  double *msd = calloc(nsamp, sizeof(double));
  double *allmsd = calloc(nsamp, sizeof(double));
  //double *avgmsd = calloc(nsamp, sizeof(double));
  double *msdcount = calloc(nsamp, sizeof(double));
	double *wx = calloc(nsamp, sizeof(double));
	double *wy = calloc(nsamp, sizeof(double));
	double *wz = calloc(nsamp, sizeof(double));
	double *vx = calloc(nsamp, sizeof(double));
	double *vy = calloc(nsamp, sizeof(double));
	double *vz = calloc(nsamp, sizeof(double));
	double *avg = calloc(wmax, sizeof(double));
	FILE *datafile;

	//Check if CH File exists//
	datafile = fopen(str2,"r");
	if(!datafile)
	{
			printf("Error - datafile not found, check data string %s\n",str2);
			exit(1);
	}

	//Read in displacements from CH File///
	for(i=0;i<nsamp;i++)
	{
			fscanf(datafile, "%d %d %d %d %lf %lf %lf\n", &ttemp, &ctemp,&dint2, &dint, &rxtemp, &rytemp, &rztemp);

			wx[i] = rxtemp;

			wy[i] = rytemp;

			wz[i] = rztemp;

	}
	fclose(datafile);

	//Sort displacements into correct order [t0,...,tmax,t0,...,tmax]

	for(j=0;j<Ncharges;j++){
		for(i=1;i<Tmax+1;i++){
			ind1 = j*Tmax+i;
			if(j>0){
				ind1 = j*(Tmax+1)+i;
			}
			ind2 = i*Ncharges + j ;
			//printf("ind1: %d ind2: %d\n",ind1,ind2);
			vx[ind1] = wx[ind2] ;
			vy[ind1] = wy[ind2] ;
			vz[ind1] = wz[ind2] ;
		}
	}

		// for(i=0;i<Ncharges*tmax;i++){
		// 	if(i%1000==0){
		// 		printf("i: %d dxt: %lf\n",i,dxt[i]);
		// 	}
		//
		// }

		// for(i=0;i<Ncharges*Tmax;i++){
		// 	printf("vx: %lf\n",vx[i]);
		// }


		//Calculate MSD//

		if(Ncharges>=25){
			nq=25 ;
		}
		else{
			nq = Ncharges;
		}
		//outputfile=fopen(outp1,"w");
		#pragma omp parallel
		{
			double *Ree2 = calloc(nq, sizeof(double));
			double *Rcount = calloc(nq, sizeof(double));

			#pragma omp for

				for(i = 1; i<wmax-MSDstart; ++i) //i is shift
				{
					for(l=0;l<nq;l++){
						Ree2[l] = 0 ;
						Rcount[l] = 0;
					}
					avg[i] = 0 ;
					for(q=0;q<nq;q++){
						for(j = MSDstart; j<wmax-i; ++j)
						{
							dx = 0.0; dy = 0.0; dz = 0.0;
							for(k = 0; k<i; ++k)
							{
								dx += vx[j+k+q*Tmax];
								dy += vy[j+k+q*Tmax];
								dz += vz[j+k+q*Tmax];
							}
							Ree2[q] += dx*dx+dy*dy+dz*dz;
							Rcount[q]++;
						}
						avg[i]+=Ree2[q]/(double)Rcount[q];
					}
				}
		}

		//Write to File
		outputfile=fopen(outp1,"w");
		for(i=1;i<wmax-1;i++){
			fprintf(outputfile,"%f,%lf\n",i*dt*printprops,avg[i]/nq);
		}
		fclose(outputfile);

}
void celllist(){
	int i,j,k,xcel,ycel,zcel,cell;
	// ljtime = omp_get_wtime();
	memset(hoc,-1,rn*rn*rn*sizeof(int));
	// for(i=0;i<rn;i++){
	// 	for(j=0;j<rn;j++){
	// 		for(k=0;k<rn;k++){
	// 			hoc[i][j][k] = -1;
	// 		}
	// 	}
	// }
	// for(i=0;i<rn*rn*rn;i++){
	// 	hoc[i] = -1;
	// }
	// ljtimetest += omp_get_wtime() - ljtime;
	// for(i=0;i<rn*rn*rn;i++){
	// 	printf("%d %d\n",i,hoc[i]);
	// }
	// exit(1);
	for(i=0;i<N;i++){
		xcel = floor((rx[i]+box_side)/rl);
		ycel = floor((ry[i]+box_side)/rl);
		zcel = floor((rz[i]+box_side)/rl);
		// printf("%d %d %d %d %lf %lf %lf\n",i,xcel,ycel,zcel,rx[i],ry[i],rz[i]);
		cell = rn*rn*xcel + rn*ycel + zcel;
		// ll[i] = hoc[xcel][ycel][zcel];
		// hoc[xcel][ycel][zcel] = i;
		ll[i] = hoc[cell];
		hoc[cell] = i;
		// printf("%d %d %d\n",i,cell,hoc[cell]);
	}
	// exit(1);
}
void LJcell(){
	int i,j,jj,k,l,xceli,yceli,zceli,xcelj,ycelj,zcelj,cell;
	double dx, dy, dz, rr, coeff, r6, ratio;
	// double ljtime = omp_get_wtime();
	#pragma omp parallel
	{
		int i,j,jj,k,l,xceli,yceli,zceli,xcelj,ycelj,zcelj,cell;
		double dx,dy,dz,rr,coeff,r6,ratio;
		int id = omp_get_thread_num();
		double *fxn = calloc(N,sizeof(double));
		double *fyn = calloc(N,sizeof(double));
		double *fzn = calloc(N,sizeof(double));
		#pragma omp for schedule(static,1)
			for(i=0;i<N;i++){
				xceli = floor((rx[i]+box_side)/rl);
				yceli = floor((ry[i]+box_side)/rl);
				zceli = floor((rz[i]+box_side)/rl);
				for(j=-1;j<2;j++){
					xcelj = xceli + j;
					xcelj -= rn*floor((double)xcelj/rn);
					for(k=-1;k<2;k++){
						ycelj = yceli + k;
						ycelj -= rn*floor((double)ycelj/rn);
						for(l=-1;l<2;l++){
							zcelj = zceli + l;
							zcelj -= rn*floor((double)zcelj/rn);
							cell = rn*rn*xcelj + rn*ycelj + zcelj;
							jj = hoc[cell];
							while(jj!=-1){
								if(i!=jj){
									dx = rx[i] - rx[jj]; dy = ry[i] - ry[jj]; dz = rz[i] - rz[jj];
									dx -= box_length*round(dx/box_length); dy -= box_length*round(dy/box_length); dz -= box_length*round(dz/box_length);
									rr = dx*dx+dy*dy+dz*dz;
									if(rr<r2cut){
										ratio = 4.00/rr;
										r6 = ratio*ratio*ratio;
										if(r6>10) r6 = 10;
										// coeff = (48*epsilon/rr)*(r6*r6-0.5*r6);
										coeff = (12*epsilon/rr)*(r6*r6-r6);
										// #pragma omp atomic
											// fx[i] += coeff*dx;
											// fy[i] += coeff*dy;
											// fz[i] += coeff*dz;
										fxn[i] += coeff*dx;
										fyn[i] += coeff*dy;
										fzn[i] += coeff*dz;
									}
								}
								jj = ll[jj];
							}
						}
					}
				}
			}
			#pragma omp critical
			{
				for(i=0;i<N;i++){
					fx[i] += fxn[i];
					fy[i] += fyn[i];
					fz[i] += fzn[i];
				}
				free(fxn); free(fyn); free(fzn);
			}
	}
}
