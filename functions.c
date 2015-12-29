#include<stdio.h>
#include<stdlib.h>

/*
 * The following function is the initialization function.
 * It takes paramters from the main file and stores them in global variables declared in the "functions.h" file.
 * The "_d" suffix in variable names below indicates that the values of the variables come with dimensions. i.e, they are real, physical values
 * The function below converts the values into dimensionless quantities by dividing with appropriate constants.
*/


void dostuff( double alpha_d, double gamma_d,double current_d, double ms_d, double eta_d, double beta_d, double vol_d, double *hext_d, double *han_d, double *hde_d,  double *Mf_d, double theta_d, double phi_d)
{
	//Magnitude of the effective magnetic field. Its usage is to normalize the fixed magnetic field vector, Mf
	double mag_mf=0;
	//Dummy variable used for loops	
	int i;
	
	//Operation index. Is used to take in input the operation the user wishes to do.	
	int op;
	//Character which upon input prompt if filled with 'y' or 'n' depending upon whether or not to shutdown pc after simulation
	//Note the user needs to be root for the system to shutdown after completion of simulation
	char shutdown;
	//Function pointer to point to function decided by the "op" variable. Both "op" and "wtd" are assigned values further in the program
	void (*wtd)();
	
	//Alpha is just the damping constant.
	alpha = alpha_d;
	//Constantof proportionality in dm/dt = gama * (m x Heff)
	gama = gamma_d;
	//The current supplied to give rise to Spin Transfer Torque and Field Like Torque
	current = current_d;
	//The magnitude of the soft magnetic field
	ms = ms_d;
	//Used for calculation of aj- the term corresponding to STT
	eta = eta_d;
	//Ratio of aj and bj. ie bj = beta*aj
	beta = beta_d;
	//Volume of the Magnetic Tunneling Junction
	vol = vol_d;
	
	//Initia values of theta and phi
	theta = theta_d;
	phi = phi_d;
	
	//Initial value of time and the size of the time step.


	/*  RATIONALE OF THE TIME STEP SIZE
	    -------------------------------
		
	 * The time step taken in real time is 0.1 Pico Second (10^-13 s)
	 * If t is a time-step then (t*4*pi*ms*gama) is the dimensionless time step
	 * 4*pi*ms*gama is equal to 0.285 * 10 ^ 12 s.
	 * Therefore 1 picoS is equal to 0.0285 in dimensionless time.
	 */

	ti=0; dt=4*M_PI*ms*gama*TIME_STEP;
	printf("Step size is %lf\n", dt);
	
	//Constants contributing to the Effective Field
        // 1. hext ------> The external applied magnetic field
	// 2. han  ------> Field arising from anisotropy. Dependent upon geometry. In this case,in z direction only.
	// 3. hde  ------> Field arising from the magnetization of the fixed magnet, Mf. Dependent on geomtery. In this case, in z direction only.
	for(i=0; i<3; i++)
	{
		hext[i] = hext_d[i];
		han[i] = han_d[i];
		hde[i] = hde_d[i];
	}

	//Making the quantities dimensionless by dividing with (4*pi*ms)
	for(i=0; i<3; i++)
	{
		hext[i] = hext[i]/(4*M_PI*ms);
		han[i] = han[i]/(4*M_PI*ms);
		hde[i] = hde[i]/(4*M_PI*ms);
	}

	//Diving Mf with its magnitude to make it dimensionless and a unit vector
	
	for(i=0; i<3; i++)
	mag_mf = mag_mf + Mf_d[i]*Mf_d[i];

	mag_mf = sqrt(mag_mf);

	for(i=0; i<3; i++)
	mf[i] = Mf_d[i]/mag_mf;
	
	//Display the options of the operations
	displayPrompt();	
	
	scanf("%d", &op);	
	printf("value of op is %d\n", op);
	//The following switch block examines the value of "op" and assigns the address of the appropriate function to "wtd".
	switch(op)
	{
		case 1:
		{printf("1 now\n");
		wtd = &heffOnly;
		break;}
		
		case 2:
		wtd = &heffDamp;
		break;
		
		case 3:
		wtd = &heffDampStt;
		break;

		case 4:
		wtd = &heffDampSttFlt;
		break;

		case 5:
		wtd = &heffDampFluc;
		break;

		case 6:
		wtd = &heffDampSttFluc;
		break;

		case 7:
		wtd = &heffDampSttFltFluc;
		break;
		
		case 8:
		wtd = &flucOnly;
		break;
		
	
		default:
		wtd = &heffDampSttFltFluc;

	}
	printf("\n\n\nValues of theta and phi: %lf, %lf\n\n\n ", theta, phi);
	printf("Do you wish to shut down after completion?[y/n]\n");
	CLEARBUF();
	scanf("%c", &shutdown);
	

	clock_t tic = clock();
	(*wtd)();
	clock_t toc = clock();
	printf("Done computing\n");
	 
   	FILE *tk = fopen("time_taken.d", "w");
   	fprintf(tk, "Elapsed time is %lf seconds\n", (double)(toc-tic)/CLOCKS_PER_SEC);
   	fclose(tk);
	//printf("Value of the stoc constant is %le", -sqrt(2*alpha*KB_TIMES_TEMP)/sqrt(gama*ms*vol));
	
	
	if(shutdown=='y' || shutdown=='Y')
	system("shutdown -h +1");
	
	else;
	return;
	
}



void displayPrompt()
{
	char *s1 = "1. Heff only";
	char *s2 = "2. Heff and damping";
	char *s3 = "3. Heff, damping and STT";
	char *s4 = "4. Heff, damping, STT and FLT";
	char *s5 = "5. Heff, damping, fluctuation";
	char *s6 = "6. Heff, damping, STT and fluctuation";
	char *s7 = "7. Heff, damping, STT, FLT and fluctuation";
	char *s8 = "8. Only Thermal Fluctuations";
	printf("Enter your choice\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n", s1, s2, s3, s4, s5, s6, s7, s8);
}

void heffOnly()
{
	printf("Into 1\n");	
	alpha = 0;
	//Temporary variables for the rk-2	
	double mtheta, mphi;
	//T is equal heff + (bj - alpha*aj)mf
	double T[3];
	//T is equal to alpha*heff + (aj +alpha *bj)mf
	double Tp[3] = {0,0,0};
	//k-term for theta in rk-2
	double kt[2];
	//k -term for phi in rk-2
	double kp[2];
	//Temporary variable
	int i;
	
	double mx, my, mz;
	mx = sin(theta)*cos(phi);
	my = sin(theta)*sin(phi);
	mz = cos(theta);

	//Creates filename
	char **filename = generateFileNames(1);
	

	//The file creating block ends here
	FILE *fp = fopen(*filename, "w");
	if(fp==NULL)
	{
		printf("Could not generate file pointer\n");
		return;

	}
	int j=1;

	for(i=0; i<3; i++)
	printf("%lf\t",heff[i]);
	printf("\n");
	while(j<NO_TIME_STEPS)
	{
		
		
		heffCalc(mx, my, mz);	
		for(i=0; i<3; i++)
		T[i] = heff[i];

 		kp[0] = phiDot(theta, phi, T , Tp) ;
         	kt[0] = thetaDot(theta, phi, T, Tp) ;

         	
         	mphi = phi + dt*(kp[0]) ;
         	mtheta = theta + dt*(kt[0])  ;

         	kp[1] = phiDot(mtheta, mphi, T, Tp);
         	kt[1] = thetaDot(mtheta, mphi, T, Tp);

         	phi = phi + dt/2.*(kp[0]  + kp[1]) ;
         	theta = theta + dt/2.*(kt[0]  + kt[1]);
		if(theta<0)
		{
			theta =-theta;	
			phi = M_PI + phi;
		}
		if(theta>M_PI)
		{
			theta = TWO_PI - theta;
			phi = phi + M_PI;
		}
		
		mx = sin(theta)*cos(phi);
		my = sin(theta)*sin(phi);
		mz = cos(theta);
		
		ti = ti+dt;
		fprintf(fp, "%.4f\t%.4f\t%.4f\t%.4f\n", ti, mx, my, mz);
		//fprintf(fp, "%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n", ti, mx, my, mz, theta);
		if(j%(NO_TIME_STEPS/10)==0)
         	printf("%d percent complete\n", (100*j/NO_TIME_STEPS));
		
		
		j++;

	}
	fclose(fp);
}

void heffDamp()
{
	//Temporary variables for the rk-2
	double mtheta, mphi;
	//T is equal heff + (bj - alpha*aj)mf
	double T[3];
	//T is equal to alpha*heff + (aj +alpha *bj)mf
	double Tp[3];
	//k-term for theta in rk-2
	double kt[2];
	//k -term for phi in rk-2
	double kp[2];
	//Temporary variable
	int i;
	double mx, my, mz;
	mx = sin(theta)*cos(phi);
	my = sin(theta)*sin(phi);
	mz = cos(theta);


	//Creates filename
	char **filename  = generateFileNames(2);
	//The file creating block ends here
	FILE *fp = fopen(*filename, "w");
	if(fp==NULL)
	{
		printf("Could not generate file pointer\n");
		return;

	}
	int j=0;
	while(j<NO_TIME_STEPS)
	{
		heffCalc(mx, my, mz);		
		for(i=0; i<3; i++)
		{
			T[i] = heff[i];
			Tp[i] = alpha*heff[i];
		}

 		kp[0] = phiDot(theta, phi, T , Tp) ;
         	kt[0] = thetaDot(theta, phi, T, Tp) ;

         	
         	mphi = phi + dt*(kp[0]) ;
         	mtheta = theta + dt*(kt[0])  ;

         	kp[1] = phiDot(mtheta, mphi, T, Tp);
         	kt[1] = thetaDot(mtheta, mphi, T, Tp);

         	phi = phi + dt/2.*(kp[0]  + kp[1]) ;
         	theta = theta + dt/2.*(kt[0]  + kt[1]);
		if(theta<0)
		{
			theta =-theta;	
			phi = M_PI + phi;
		}
		if(theta>M_PI)
		{
			theta = TWO_PI - theta;
			phi = phi + M_PI;
		}
		
		mx = sin(theta)*cos(phi);
		my = sin(theta)*sin(phi);
		mz = cos(theta);
	
		ti = ti+dt;
		
		fprintf(fp, "%.4f\t%.4f\t%.4f\t%.4f\n", ti, mx, my, mz);
		//fprintf(fp, "%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n", ti, mx, my, mz, theta);
		if(j%(NO_TIME_STEPS/10)==0)
         	printf("%d percent complete\n", (100*j/NO_TIME_STEPS));
         	j++;

	}
	fclose(fp);
}

void heffDampStt()
{
	
	//Temporary variables for the rk-2
	double mtheta, mphi;
	//T is equal heff + (bj - alpha*aj)mf
	double T[3];
	//T is equal to alpha*heff + (aj +alpha *bj)mf
	double Tp[3];
	//k-term for theta in rk-2
	double kt[2];
	//k -term for phi in rk-2
	double kp[2];
	//Temporary variable
	int i;
	double mx, my, mz;
	mx = sin(theta)*cos(phi);
	my = sin(theta)*sin(phi);
	mz = cos(theta);

	//Creates filename
	char **filename = generateFileNames(3);
	//The file creating block ends here
	FILE *fp = fopen(*filename, "w");
	if(fp ==NULL)
	{
		printf("Could not generate file pointer\n");
		return;

	}
	int j=1;inversion_time[k]=ti;
	while(j<NO_TIME_STEPS)
	{
		heffCalc(mx, my, mz);		
		aj_var = aj(theta, phi);
		//if(j<20)
		//printf("Aj is %e\n", aj_var);		
		for(i=0; i<3; i++)
		{
			T[i] = heff[i] - alpha*mf[i]*aj_var;
			Tp[i] = alpha*heff[i] + aj_var*mf[i];
		}


 		kp[0] = phiDot(theta, phi, T , Tp) ;
         	kt[0] = thetaDot(theta, phi, T, Tp) ;

         	
         	mphi = phi + dt*(kp[0]) ;
         	mtheta = theta + dt*(kt[0])  ;

         	kp[1] = phiDot(mtheta, mphi, T, Tp);
         	kt[1] = thetaDot(mtheta, mphi, T, Tp);

         	phi = phi + dt/2.*(kp[0]  + kp[1]) ;
         	theta = theta + dt/2.*(kt[0]  + kt[1]);
		if(theta<0)
		{
			theta =-theta;	
			phi = M_PI + phi;
		}
		if(theta>M_PI)
		{
			theta = TWO_PI - theta;
			phi = phi + M_PI;
		}
		
		mx = sin(theta)*cos(phi);
		my = sin(theta)*sin(phi);
		mz = cos(theta);
		
		ti = ti+dt;
		
		fprintf(fp, "%.4f\t%.4f\t%.4f\t%.4f\n", ti, mx, my, mz);
		
		if(j%(NO_TIME_STEPS/10)==0)
         	printf("%d percent complete\n", (100*j/NO_TIME_STEPS));
         	j++;
	}
	fclose(fp);
}

void heffDampSttFlt()
{
	
	//Temporary variables for the rk-2
	double mtheta, mphi;
	//T is equal heff + (bj - alpha*aj)mf
	double T[3];
	//T is equal to alpha*heff + (aj +alpha *bj)mf
	double Tp[3];
	//k-term for theta in rk-2
	double kt[2];
	//k -term for phi in rk-2
	double kp[2];
	//Temporary variable
	int i;
	double mx, my, mz;

	mx = sin(theta)*cos(phi);
	my = sin(theta)*sin(phi);
	mz = cos(theta);
	//Creates filename
	char **filename = generateFileNames(4);
	//The file creating block ends here
	FILE *fp = fopen(*filename, "w");
	if(fp==NULL)
	{
		printf("Could not generate file pointer\n");
		return;

	}	
	int j=0;
	aj_var = aj(theta, phi);
	printf("Value of aj in first iteration is %.4f\n", aj_var);
	while(j<NO_TIME_STEPS)
	{
		heffCalc(mx, my, mz);		
		aj_var = aj(theta, phi);
		bj_var = beta*aj_var;		
		for(i=0; i<3; i++)
		{
			T[i] = heff[i] + (bj_var - alpha*aj_var)*mf[i];
			Tp[i] = alpha*heff[i] + (aj_var+ alpha*bj_var)*mf[i];
		}

 		kp[0] = phiDot(theta, phi, T , Tp) ;
         	kt[0] = thetaDot(theta, phi, T, Tp) ;

         	
         	mphi = phi + dt*(kp[0]) ;
         	mtheta = theta + dt*(kt[0])  ;

         	kp[1] = phiDot(mtheta, mphi, T, Tp);
         	kt[1] = thetaDot(mtheta, mphi, T, Tp);

         	phi = phi + dt/2.*(kp[0]  + kp[1]) ;
         	theta = theta + dt/2.*(kt[0]  + kt[1]);	
		if(theta<0)
		{
			theta =-theta;	
			phi = M_PI + phi;
		}
		if(theta>M_PI)
		{
			theta = TWO_PI - theta;
			phi = phi + M_PI;
		}
		
		mx = sin(theta)*cos(phi);
		my = sin(theta)*sin(phi);
		mz = cos(theta);
		

		ti = ti+dt;
		fprintf(fp, "%.4f\t%.4f\t%.4f\t%.4f\n", ti, mx, my, mz);
		
		if(j%(NO_TIME_STEPS/10)==0)
         	printf("%d percent complete\n", (100*j/NO_TIME_STEPS));
         	j++;

	}
	fclose(fp);
}

void heffDampFluc()
{
	//Temporary variables for the rk-2
	double mtheta, mphi;
	//T is equal heff 
	double T[3];
	//T is equal to alpha*heff 
	double Tp[3];
	//k-term for theta in rk-2
	double kt[2];
	//k -term for phi in rk-2
	double kp[2];
        //k -term for stochastic theta
	double st[2];
	//k- term for stochastic phi
	double sp[2];
	//Temporary variable
	int i, j, k;

	//Following are the variables for the ensemble calculation. variable names are self explainatory	
	double *mx = malloc(NO_ENSEMBLES*sizeof(double)), *my=malloc(NO_ENSEMBLES*sizeof(double)), *mz = malloc(NO_ENSEMBLES*sizeof(double));

	double *theta_en = malloc(NO_ENSEMBLES*sizeof(double)), *phi_en = malloc(NO_ENSEMBLES*sizeof(double));

	int *inversion_flag=malloc(NO_ENSEMBLES*sizeof(int));

	double *inversion_time = malloc(NO_ENSEMBLES*sizeof(double));

	double *mxavg = malloc(NO_TIME_STEPS*sizeof(double)), *myavg=malloc(NO_TIME_STEPS*sizeof(double)), *mzavg=malloc(NO_TIME_STEPS*sizeof(double));
	//Creates filename
	char **filename = generateFileNames(5);

	//The file creating block ends here
	FILE *fp = fopen(*filename, "w");
	if(fp==NULL)
	{
		printf("Could not generate file pointer\n");
		return;

	}

	
	for(k=0; k<NO_ENSEMBLES; k++)
	{
		theta_en[k] = theta;
		phi_en[k] = phi;
		mx[k] = sin(theta)*cos(phi);
		my[k] = sin(theta)*sin(phi);
		mz[k] = cos(theta);

	}

	for(j=1; j<NO_TIME_STEPS; j++)
	{

		
		for(k=0; k<NO_ENSEMBLES; k++)
		{
		
			
			heffCalc(mx[k], my[k], mz[k]);
			for(i=0; i<3; i++)
			{
				T[i] = heff[i] ;
				Tp[i] = alpha*heff[i];
			}
 			
			kp[0] = phiDot(theta_en[k], phi_en[k], T , Tp) ;
         		kt[0] = thetaDot(theta_en[k], phi_en[k], T, Tp) ;

			sp[0] = stocPhiDot(theta_en[k], phi_en[k], &seed);
         		st[0] = stocThetaDot(theta_en[k], phi_en[k], &seed);

         	
         		mphi = phi_en[k] + dt*(kp[0]) + sqrt(dt)*sp[0] ;
         		mtheta = theta_en[k] + dt*(kt[0]) + sqrt(dt)*st[0] ;

         		kp[1] = phiDot(mtheta, mphi, T, Tp);
         		kt[1] = thetaDot(mtheta, mphi, T, Tp);
			

			sp[1] = stocPhiDot(mtheta, mphi, &seed);
         		st[1] = stocThetaDot(mtheta, mphi, &seed);

         		phi_en[k] = phi_en[k] + dt/2.*(kp[0]  + kp[1]) + 1/2.*(sp[0] + sp[1])*sqrt(dt);
         		theta_en[k] = theta_en[k] + dt/2.*(kt[0]  + kt[1]) + 1/2.*(st[0] + st[1])*sqrt(dt);

			if(theta_en[k]<0)
			{
				theta_en[k] =-theta_en[k];	
				phi_en[k] = M_PI + phi_en[k];
			}
			if(theta_en[k]>M_PI)
			{
				theta_en[k] = TWO_PI - theta_en[k];
				phi_en[k] = phi_en[k] + M_PI;
			}
		
			mx[k] = sin(theta_en[k])*cos(phi_en[k]);
			my[k] = sin(theta_en[k])*sin(phi_en[k]);
			mz[k] = cos(theta_en[k]);

			if(fabs(mz[k])>0.95 && inversion_flag[k]==0)
			{
				//printf("mz[%d] is %lf and time is %.4f\n", k, mz[k], ti);				
				inversion_time[k]=ti;
				inversion_flag[k]=1;
			}


		}

 		
		mxavg[j] = 0;
 		myavg[j] = 0;
 		mzavg[j] = 0;
		

		for(k=0; k<NO_ENSEMBLES; k++)
 		{
 			mxavg[j]+= mx[k];
 			myavg[j]+= my[k];
 			mzavg[j]+= mz[k];
		
 		}
 		
 		mxavg[j]/=(double)NO_ENSEMBLES;
 		myavg[j]/=(double)NO_ENSEMBLES;
 		mzavg[j]/=(double)NO_ENSEMBLES;
		

		ti = ti+dt;

		if(j%(NO_TIME_STEPS/10)==0)
         	printf("%d percent complete\n", (100*j/NO_TIME_STEPS));

	}

	ti=0.0;

	for(j=0; j<NO_TIME_STEPS; j++)
	{
		fprintf(fp, "%.4f\t%.4f\t%.4f\t%.4f\n", ti, mxavg[j], myavg[j], mzavg[j]);
		ti = ti+dt;
	}

	storeInversionTime(inversion_time, *(filename+1));
	free(mx);
	free(my);
	free(mz);
	free(inversion_flag);
	free(inversion_time);
	free(theta_en);
	free(phi_en);
	fclose(fp);

}
void heffDampSttFluc()
{
	//Temporary variables for the rk-2
	double mtheta, mphi;
	//T is equal heff - alpha*aj*mf
	double T[3];
	//T is equal to alpha*heff + aj*mf
	double Tp[3];
	//k-term for theta in rk-2
	double kt[2];
	//k -term for phi in rk-2
	double kp[2];
        //k -term for stochastic theta
	double st[2];
	//k- term for stochastic phi
	double sp[2];
	//Temporary variable
	int i, j, k;

	//Following are the variables for the ensemble calculation. variable names are self explainatory	
	double *mx = malloc(NO_ENSEMBLES*sizeof(double)), *my=malloc(NO_ENSEMBLES*sizeof(double)), *mz = malloc(NO_ENSEMBLES*sizeof(double));

	double *theta_en = malloc(NO_ENSEMBLES*sizeof(double)), *phi_en = malloc(NO_ENSEMBLES*sizeof(double));

	int *inversion_flag=malloc(NO_ENSEMBLES*sizeof(int));

	double *inversion_time = malloc(NO_ENSEMBLES*sizeof(double));

	double *mxavg = malloc(NO_TIME_STEPS*sizeof(double)), *myavg=malloc(NO_TIME_STEPS*sizeof(double)), *mzavg=malloc(NO_TIME_STEPS*sizeof(double));
	//Creates filename
	char **filename = generateFileNames(6);

	//The file creating block ends here
	FILE *fp = fopen(*filename, "w");
	if(fp==NULL)
	{
		printf("Could not generate file pointer\n");
		return;

	}

	
	for(k=0; k<NO_ENSEMBLES; k++)
	{
		theta_en[k] = theta;
		phi_en[k] = phi;
		mx[k] = sin(theta)*cos(phi);
		my[k] = sin(theta)*sin(phi);
		mz[k] = cos(theta);

	}
	mxavg[0] = mx[0]; myavg[0]=my[0]; mzavg[0]=mz[0];
	for(j=1; j<NO_TIME_STEPS; j++)
	{

		
		for(k=0; k<NO_ENSEMBLES; k++)
		{
		
			heffCalc(mx[k], my[k], mz[k]);			
			aj_var = aj(theta_en[k], phi_en[k]);	
			for(i=0; i<3; i++)
			{
				T[i] = heff[i] - alpha*mf[i]*aj_var;
				Tp[i] = alpha*heff[i] + aj_var*mf[i];
			}
	
 			kp[0] = phiDot(theta_en[k], phi_en[k], T , Tp) ;
         		kt[0] = thetaDot(theta_en[k], phi_en[k], T, Tp) ;

			sp[0] = stocPhiDot(theta_en[k], phi_en[k], &seed);
         		st[0] = stocThetaDot(theta_en[k], phi_en[k], &seed);

         	
         		mphi = phi_en[k] + dt*(kp[0]) + sqrt(dt)*sp[0] ;
         		mtheta = theta_en[k] + dt*(kt[0]) + sqrt(dt)*st[0] ;

         		kp[1] = phiDot(mtheta, mphi, T, Tp);
         		kt[1] = thetaDot(mtheta, mphi, T, Tp);
			

			sp[1] = stocPhiDot(mtheta, mphi, &seed);
         		st[1] = stocThetaDot(mtheta, mphi, &seed);

         		phi_en[k] = phi_en[k] + dt/2.*(kp[0]  + kp[1]) + 1/2.*(sp[0] + sp[1])*sqrt(dt);
         		theta_en[k] = theta_en[k] + dt/2.*(kt[0]  + kt[1]) + 1/2.*(st[0] + st[1])*sqrt(dt);
			
			if(theta_en[k]<0)
			{
				theta_en[k] =-theta_en[k];	
				phi_en[k] = M_PI + phi_en[k];
			}
			if(theta_en[k]>M_PI)
			{
				theta_en[k] = TWO_PI - theta_en[k];
				phi_en[k] = phi_en[k] + M_PI;
			}	
		
			mx[k] = sin(theta_en[k])*cos(phi_en[k]);
			my[k] = sin(theta_en[k])*sin(phi_en[k]);
			mz[k] = cos(theta_en[k]);

			if(fabs(mz[k])>0.95 && inversion_flag[k]==0)
			{
				//printf("mz[%d] is %lf and time is %.4f\n", k, mz[k], ti);				
				inversion_time[k]=ti;
				inversion_flag[k]=1;
			}
		}

		mxavg[j] = 0;
 		myavg[j] = 0;
 		mzavg[j] = 0;
		

		for(k=0; k<NO_ENSEMBLES; k++)
 		{
 			mxavg[j]+= mx[k];
 			myavg[j]+= my[k];
 			mzavg[j]+= mz[k];
		
 		}
 		
 		mxavg[j]/=(double)NO_ENSEMBLES;
 		myavg[j]/=(double)NO_ENSEMBLES;
 		mzavg[j]/=(double)NO_ENSEMBLES;
		

		ti = ti+dt;

		if(j%(NO_TIME_STEPS/10)==0)
         	printf("%d percent complete\n", (100*j/NO_TIME_STEPS));
//
	}

	ti=0.0;

	for(j=0; j<NO_TIME_STEPS; j++)
	{
		fprintf(fp, "%.4f\t%.4f\t%.4f\t%.4f\n", ti, mxavg[j], myavg[j], mzavg[j]);
		ti = ti+dt;
	}

	storeInversionTime(inversion_time, *(filename+1));
	free(mx);
	free(my);
	free(mz);
	free(inversion_flag);
	free(inversion_time);
	free(theta_en);
	free(phi_en);
	fclose(fp);



}
void heffDampSttFltFluc()
{
	//printf("Entering the main func\n");	
	//Temporary variables for the rk-2
	double mtheta, mphi;
	//T is equal heff + (bj - alpha*aj)mf
	double T[3];
	//T is equal to alpha*heff + (aj +alpha *bj)mf
	double Tp[3];
	//k-term for theta in rk-2
	double kt[2];
	//k -term for phi in rk-2
	double kp[2];
        //k -term for stochastic theta
	double st[2];
	//k- term for stochastic phi
	double sp[2];
	//Temporary variable
	int i, j, k;

	
	//Following are the variables for the ensemble calculation. variable names are self explainatory
	double *mx = malloc(NO_ENSEMBLES*sizeof(double)), *my=malloc(NO_ENSEMBLES*sizeof(double)), *mz = malloc(NO_ENSEMBLES*sizeof(double));

	double *theta_en = malloc(NO_ENSEMBLES*sizeof(double)), *phi_en = malloc(NO_ENSEMBLES*sizeof(double));

	int *inversion_flag=malloc(NO_ENSEMBLES*sizeof(int));

	double *inversion_time = malloc(NO_ENSEMBLES*sizeof(double));


	//printf("switching theta : %e and switching phi: %e\n", switching_theta, switching_phi);
	double reverse_switching_polarity;
	if(hext[2]<0)
	reverse_switching_polarity=-1;
	else
	reverse_switching_polarity=1;

	double *mxavg = malloc(NO_TIME_STEPS*sizeof(double)), *myavg=malloc(NO_TIME_STEPS*sizeof(double)), *mzavg=malloc(NO_TIME_STEPS*sizeof(double));
	//Creates filename
	char **filename = generateFileNames(7);

	//The file creating block ends here
	FILE *fp = fopen(*filename, "w");
	if(fp==NULL)
	{
		printf("Could not generate file pointer\n");
		return;

	}	

	double switching_tolerance = 8e-5;
	//double tol_mz = 0.2;
	//FILE *ens = fopen("ens.d", "w");
	
	for(k=0; k<NO_ENSEMBLES; k++)
	{
		theta_en[k] = theta;
		phi_en[k] = phi;
		mx[k] = sin(theta_en[k])*cos(phi_en[k]);
		my[k] = sin(theta_en[k])*sin(phi_en[k]);
		mz[k] = cos(theta_en[k]);
		inversion_time[k] = 0;

	}
	mxavg[0] = mx[0]; myavg[0]=my[0]; mzavg[0]=mz[0];
	for(j=1; j<NO_TIME_STEPS; j++)
	{

		
		for(k=0; k<NO_ENSEMBLES; k++)
		{
		
			heffCalc(mx[k], my[k], mz[k]);			
			aj_var = aj(theta_en[k], phi_en[k]);
			bj_var = beta*aj_var;		
			for(i=0; i<3; i++)
			{
				T[i] = heff[i] + (bj_var - alpha*aj_var)*mf[i];
				Tp[i] = alpha*heff[i] + (aj_var+ alpha*bj_var)*mf[i];
			}

	
 			kp[0] = phiDot(theta_en[k], phi_en[k], T , Tp) ;
         		kt[0] = thetaDot(theta_en[k], phi_en[k], T, Tp) ;

			sp[0] = stocPhiDot(theta_en[k], phi_en[k], &seed);
         		st[0] = stocThetaDot(theta_en[k], phi_en[k], &seed);

         	
         		mphi = phi_en[k] + dt*(kp[0]) + sqrt(dt)*sp[0] ;
         		mtheta = theta_en[k] + dt*(kt[0])+ sqrt(dt)*st[0] ;

         		kp[1] = phiDot(mtheta, mphi, T, Tp);
         		kt[1] = thetaDot(mtheta, mphi, T, Tp);
			

 			sp[1] = stocPhiDot(mtheta, mphi, &seed);
         		st[1] = stocThetaDot(mtheta, mphi, &seed);
			//printf("Deterministic term: %e\n Stochastic term: %e\n", dt/2.*(kp[0]  + kp[1]), 1/2.*(sp[0] + sp[1])*sqrt(4*M_PI*ms*gama*dt));

         		phi_en[k] = phi_en[k] + dt/2.*(kp[0]  + kp[1]) + 1/2.*(sp[0] + sp[1])*sqrt(dt);
         		theta_en[k] = theta_en[k] + dt/2.*(kt[0]  + kt[1]) + 1/2.*(st[0] + st[1])*sqrt(dt);

			if(theta_en[k]<0)
			{
				theta_en[k] =-theta_en[k];	
				phi_en[k] = M_PI + phi_en[k];
			}
			if(theta_en[k]>M_PI)
			{
				theta_en[k] = TWO_PI - theta_en[k];
				phi_en[k] = phi_en[k] + M_PI;
			}
			if(phi_en[k]>TWO_PI)
			{
				phi_en[k] = phi_en[k] - TWO_PI;
			}
			
			if(phi_en[k]<0)
			{
				phi_en[k] = phi_en[k] + TWO_PI;
			}	
		
			mx[k] = sin(theta_en[k])*cos(phi_en[k]);
			my[k] = sin(theta_en[k])*sin(phi_en[k]);
			mz[k] = cos(theta_en[k]);

			/*if(fabs(my[k]-reverse_switching_polarity)<tol_mz && inversion_flag[k]==0)
			{
				inversion_time[k]=ti;
				inversion_flag[k]=1;
			}*/
			
			if(fabs(mz[k])>0.95 && inversion_flag[k]==0)
			{
				//printf("mz[%d] is %lf and time is %.4f\n", k, mz[k], ti);				
				inversion_time[k]=ti;
				inversion_flag[k]=1;
			}

		
		}

		mxavg[j] = 0;
 		myavg[j] = 0;
 		mzavg[j] = 0;
		

		for(k=0; k<NO_ENSEMBLES; k++)
 		{
 			mxavg[j]+= mx[k];
 			myavg[j]+= my[k];
 			mzavg[j]+= mz[k];
		
 		}
 		
 		mxavg[j]/=(double)NO_ENSEMBLES;
 		myavg[j]/=(double)NO_ENSEMBLES;
 		mzavg[j]/=(double)NO_ENSEMBLES;
		
		ti = ti+dt;


		if(j%(NO_TIME_STEPS/10)==0)
         	printf("%d percent complete\n", (100*j/NO_TIME_STEPS));

	}
	ti=0.0;

	for(j=0; j<NO_TIME_STEPS; j++)
	{
		fprintf(fp, "%.4f\t%.4f\t%.4f\t%.4f\n", ti, mxavg[j], myavg[j], mzavg[j]);
		ti = ti+dt;
	}

	storeInversionTime(inversion_time, *(filename+1));
	free(mx);
	free(my);
	free(mz);
	free(inversion_flag);
	free(inversion_time);
	free(theta_en);
	free(phi_en);
	fclose(fp);

}

void flucOnly()
{
		
	//Temporary variables for the rk-2
	double mtheta, mphi;
	//T is equal heff + (bj - alpha*aj)mf
	double T[3];
	//T is equal to alpha*heff + (aj +alpha *bj)mf
	double Tp[3];
	//k-term for theta in rk-2
	double kt[2];
	//k -term for phi in rk-2
	double kp[2];
        //k -term for stochastic theta
	double st[2];
	//k- term for stochastic phi
	double sp[2];
	//Temporary variable
	int i;
	
	double mx[NO_ENSEMBLES], my[NO_ENSEMBLES], mz[NO_ENSEMBLES];
	double mxavg, myavg, mzavg;
	//Creates filename
	char **filename = generateFileNames(8);
	//The file creating block ends here
	FILE *fp = fopen(*filename, "w");
	int j, k;

	double theta_en[NO_ENSEMBLES], phi_en[NO_ENSEMBLES];
	for(k=0; k<NO_ENSEMBLES; k++)
	{
		theta_en[k] = theta;
		phi_en[k] = phi;

	}

	for(j=1; j<NO_TIME_STEPS; j++)
	{
		
		for(k=0; k<NO_ENSEMBLES; k++)
		{

			


			sp[0] = stocPhiDot(theta_en[k], phi_en[k], &seed);
         		st[0] = stocThetaDot(theta_en[k], phi_en[k], &seed);

         	
         		mphi = phi_en[k] + sqrt(dt)*sp[0] ;
         		mtheta = theta_en[k] + sqrt(dt)*st[0] ;
			
			sp[1] = stocPhiDot(mtheta, mphi, &seed);
         		st[1] = stocThetaDot(mtheta, mphi, &seed);

         		phi_en[k] = phi_en[k] +  (1/2.)*(sp[0] + sp[1])*sqrt(dt);
         		theta_en[k] = theta_en[k] + (1/2.)*(st[0] + st[1])*sqrt(dt);

			if(theta_en[k]<0)
			{
				theta_en[k] =-theta_en[k];	
				phi_en[k] = M_PI + phi_en[k];
			}
			if(theta_en[k]>M_PI)
			{
				theta_en[k] = TWO_PI - theta_en[k];
				phi_en[k] = phi_en[k] + M_PI;
			}	
		
			mx[k] = sin(theta_en[k])*cos(phi_en[k]);
			my[k] = sin(theta_en[k])*sin(phi_en[k]);
			mz[k] = cos(theta_en[k]);
		}

		mxavg = 0;
 		myavg = 0;
 		mzavg = 0;

		for(k=0; k<NO_ENSEMBLES; k++)
 		{
 			mxavg+=mx[k];
 			myavg+=my[k];
 			mzavg+=mz[k];
 		}
 		
 		mxavg = mxavg/NO_ENSEMBLES;
 		myavg = myavg/NO_ENSEMBLES;
 		mzavg = mzavg/NO_ENSEMBLES;
		

		ti = ti+dt;
		fprintf(fp, "%.4f\t%.4f\t%.4f\t%.4f\n", ti, mxavg, myavg, mzavg);
		
		if(j%(NO_TIME_STEPS/10)==0)
         	printf("%d percent complete\n", (100*j/NO_TIME_STEPS));

	}
	fclose(fp);


}


double thetaDot(double t, double p, double *T, double *Tp)
{

	double K = 1/(1+alpha*alpha);
	return ( -K*(( sin(p)*T[0] ) - (cos(p)*T[1]) - (cos(p)*cos(t)*Tp[0]) - (sin(p)*cos(t)*Tp[1]) + (sin(t)*Tp[2]) ));

}


double phiDot(double t, double p, double *T, double *Tp)
{
	double K = 1/(1+alpha*alpha);
	return ( -K*( (cos(p)*cos(t)*T[0]/sin(t)) + (T[1]*cos(t)*sin(p)/sin(t))  - T[2] + (Tp[0]*sin(p)/sin(t)) - (Tp[1]*cos(p)/sin(t)) ));	
}

double stocThetaDot(double theta, double phi, long *seed)
{
	double val = gauss(seed)*(sin(phi) - alpha*cos(theta)*cos(phi)) - gauss(seed)*(cos(phi) - alpha*cos(theta)*sin(phi)) + alpha*gauss(seed)*sin(theta);
	double K = 1/(1+alpha*alpha);
	val =  -K*sqrt(2*alpha*KB_TIMES_TEMP)/sqrt(gama*ms*vol)*val;
	//printf("Hlf is %.4f\n", val/4*M_PI*ms);
	return val/(4*M_PI*ms)*sqrt(4*M_PI*ms*gama);
	

}

double stocPhiDot(double theta, double phi, long *seed)
{
	double val = gauss(seed)*(cos(theta)*cos(phi) + alpha*sin(phi)) + gauss(seed)*(cos(theta)*sin(phi) - alpha*cos(phi)) - gauss(seed)  ;
	double K = 1/(1+alpha*alpha);
	val = -K*sqrt(2*alpha*KB_TIMES_TEMP)/sqrt(gama*ms*vol)*val/sin(theta);
	val =  val/(4*M_PI*ms)*sqrt(4*M_PI*ms*gama);
	//printf("Value of stocPhiDot is %e\n", val);
	return val;
	
}
double aj(double t, double p)
{
	double dot_product = sin(t)*cos(p)*mf[0] + sin(t)*sin(p)*mf[1] + cos(t)*mf[2];	
	double g = eta/(1 + 0.1*dot_product);
	double value = (HCROSS*g*current)/(2*ECHARGE*ms*vol);
	value = value/(4*M_PI*ms);
	return value;

}

void heffCalc(double mx, double my, double mz)
{
	heff[0] = hext[0] - (hde[0] - han[0])*mx;
	heff[1] = hext[1] - (hde[1] - han[1])*my;
	heff[2] = hext[2] - (hde[2] - han[2])*mz;

}


void storeInversionTime(double *inversion_time, char *filename)
{
	
	FILE *fp = fopen(filename, "w");
	int i;
	for(i=0; i<NO_ENSEMBLES; i++)
	{
		fprintf(fp, "%d\t%.4f\n", i+1, inversion_time[i]);

	}
	fclose(fp);

}

void storeDwellingTime(double *dwelling_time, char *filename)
{
	FILE *fp = fopen(filename, "w");
	int i;
	for(i=0; i<NO_ENSEMBLES; i++)
	{
		fprintf(fp, "%d\t%.4f\n", i+1, dwelling_time[i]);

	}
	fclose(fp);
	

}



