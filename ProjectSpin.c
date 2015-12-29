#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"ProjectSpin.h"


void ProjectSpin(double alphaDim, double gamaDim, double currentDim, double MsDim,
	 			 double etaDim ,double betaDim, long double volumeDim, double  *hexDim, 
	 			 double *h_anisoDim, double *h_demagDim, double *MfDim, double thetaDim, double phiDim ,
	 			 long No_of_Ensemble, double time_of_sim,int choice, double toleranceDim, double tempDim,
	 			 double time_step_sizeDim)
	{

	
	//variables required
	int i;
	
	//Variable to store dimension less form of given argument to this function 
	 alpha = alphaDim;
	 gama = gamaDim;
	 current = currentDim;
	 Ms = MsDim;
	 eta = etaDim;
	 beta = betaDim;
	 volume = volumeDim;
	 NO_OF_ENSEMBLE = No_of_Ensemble;
	 NO_OF_TIME_STEPS = (long)(time_of_sim/time_step_sizeDim);
	 //newly added
	 temp = tempDim;
	 tolerance = toleranceDim;
	 time_step_size = 4 * M_PI * Ms* gama * time_step_sizeDim;
	 //NO_OF_TIME_STEPS /= time_step_size;
	 KB_TIMES_TEMPERATURE = KB * temp;
	 //done
	 printf("Number offtime steps : %ld\n",NO_OF_TIME_STEPS);
	 theta = thetaDim;
	 phi = phiDim;
	 double M[3];
	 M[0] = sin(theta)*cos(phi);
	 M[1] = sin(theta)*sin(phi);
	 M[2] = cos(theta);
	 
	 		//heff dimensionless
	 for(i = 0 ; i<3 ; i++){
	 	h_demag[i] = h_demagDim[i];
	 	h_aniso[i] = h_anisoDim[i];;
	 	hex[i] = hexDim[i];
	 	//heff[i] = ( hexDim[i]  - (   h_demagDim[i] - h_anisoDim[i] ) * M[i] )/(4 * M_PI * Ms);
	 }
	 
		//calculating magnitude of Mf
	double magMf = 0;
	for(i = 0 ; i < 3 ; i++){
		magMf += (MfDim[i] * MfDim[i]);
	}
	magMf = sqrt(magMf);
	//printf("Magnitude of Mf: %lf\n",magMf);
		//Mf unit vector
	for(i = 0 ; i < 3 ; i++){
		Mf[i] = MfDim[i]/magMf;
		
	}

	
	//T initial and time steps
	  ti = 0;
	  dt = time_step_size;
	  printf("%lf\n",dt);
	// dt = dimension time step * 4 * PI * gama * ms
	void(*wtd)();
	switch(choice){
		case 1: wtd = &Heffonly;
				break;
				
		case 2: wtd = &ThermalFluc;
				break;
				
		case 3: wtd = &Heff_Damp;
				break;
				
		case 4: wtd = &Heff_Damp_STT;
				break;
				
		case 5: wtd = &Heff_Damp_Fluc;
				break;
				
		case 6: wtd = &Heff_Damp_STT_Fluc;
				break;
				
		case 7: wtd = &Heff_Damp_STT_FLT;
				break;
				
		case 8: wtd = &Heff_Damp_STT_FLT_Fluc;
				break;
				
		default : wtd = &Heff_Damp_STT;
				break;
	};
	
	(*wtd)();
}


//1. 	Only H-effective term;
void Heffonly(){
	printf("Heffective only:\n");
	return;
}

//2.	Thermal Fluctuation only
void ThermalFluc(){
	return;
}

//3.	H-effective and  Damping 
void Heff_Damp(){
	return;
}

//4.	H-effective , Damping and Spin Transfer Torque
void Heff_Damp_STT(){
	return;
}

//5.	Heffective , Damping and Fluctuation
void Heff_Damp_Fluc(){
	return;
}

//6.	Heffective , Damping , spin Transfer Torque and Fluctuation
void Heff_Damp_STT_Fluc(){
	return;
}

//7.	H-effective , Damping, Spin transfer Torque and Field like Torque
void Heff_Damp_STT_FLT(){
	return;
}

//8.	H-effective , Damping, Spin transfer Torque , Field like Torque and Fluctuation
void  Heff_Damp_STT_FLT_Fluc(){

	//printf("H-effective , Damping, Spin transfer Torque , Field like Torque and Fluctuation:\n");
	printf("Values of theta and phi are %.4f, %.4f\n", theta, phi);
	
	//Temporary variables reqired:
	int i,j, k;
	// T is equal heff + (bj - alpha*aj)mf
	double T[3];
	//Tp is equal to alpha*heff + (aj +alpha *bj)mf
	double Tp[3];
	//k-term for theta in rk-2
	double k_theta[2];
	//k -term for phi in rk-2
	double k_phi[2];
	//k -term for stochastic theta
	double stoch_theta[2];
	//k -term for stochastic phi
	double stoch_phi[2];
	//Mz_previous
	double Mz_previous = cos(theta) ;
	//double theta_en[NO_OF_ENSEMBLE], phi_en[NO_OF_ENSEMBLE];
	//temporary variable for rk2 method
	double Mphi,Mtheta;
	
	/* edited
	double Mx[NO_OF_ENSEMBLE],My[NO_OF_ENSEMBLE],Mz[NO_OF_ENSEMBLE];
	*/
	//edited
	double *Mx;
	Mx = (double *)malloc(NO_OF_ENSEMBLE * sizeof(double));
	double *My;
	My = (double *)malloc(NO_OF_ENSEMBLE * sizeof(double));
	double *Mz ;
	Mz = (double *)malloc(NO_OF_ENSEMBLE * sizeof(double));
	double *theta_en;
	theta_en = (double *) malloc(NO_OF_ENSEMBLE * sizeof(double));
	double *phi_en;
	phi_en = (double *) malloc(NO_OF_ENSEMBLE * sizeof(double));
	//till here
	
	//int inversion_flag[NO_OF_ENSEMBLE];
	/*edited
	double inversion_time[NO_OF_ENSEMBLE];
	double dwelling_time[NO_OF_ENSEMBLE];
	double ta[NO_OF_ENSEMBLE], tb[NO_OF_ENSEMBLE];
	double dwelling_flag[NO_OF_ENSEMBLE];
	*/
	double *inversion_time;
	inversion_time = (double *)malloc(NO_OF_ENSEMBLE * sizeof(double));
	
	int *inversion_flag;
	inversion_flag = (int *)malloc(NO_OF_ENSEMBLE * sizeof(int));
	
	/*
	double *dwelling_time;
	dwelling_time = (double *)malloc(NO_OF_ENSEMBLE * sizeof(double));
	*/
	/*
	double *ta, *tb;
	ta = (double *)malloc(NO_OF_ENSEMBLE * sizeof(double));
	tb = (double *)malloc(NO_OF_ENSEMBLE * sizeof(double));
	*/
	
	
	//int temp,temp2;
	
	double *Mx_avg, *My_avg, *Mz_avg;
	Mx_avg = (double *)malloc(NO_OF_TIME_STEPS * sizeof(double));
	My_avg = (double *)malloc(NO_OF_TIME_STEPS * sizeof(double));
	Mz_avg = (double *)malloc(NO_OF_TIME_STEPS * sizeof(double));
	
	// *****************************************Edit from here*************
	//double tol_Mz = 8e-2;
	printf("\nHere**********************\n");
	char **filename = generateFileNames(7);
	//printf("\nHere**********************\n");
	FILE *fp = fopen(*filename,"w");
	//printf("\nHere**********************\n");
	//FILE *fp2 = fopen(*(filename+1),"w");
	if(fp == NULL){
		printf("Could not create file pointer \n");
		return;
	} 
	//FILE *ts = fopen("tswitch.d", "w");
	//FILE *dw = fopen("dwel_time.d", "w");
		
	
	// edit this by taking value from the file
	//double tol_mz;
	/*
	double reverse_switching_polarity;
	if(theta<=1.57)
	reverse_switching_polarity=-1;
	else
	reverse_switching_polarity=1;
	*/
	
	for(k=0; k < NO_OF_ENSEMBLE; k++)
	{
		theta_en[k] = theta;
		phi_en[k] = phi;
		Mx[k] = sin(theta_en[k])*cos(phi_en[k]);
		My[k] = sin(theta_en[k])*sin(phi_en[k]);
		Mz[k] = cos(theta_en[k]);
		inversion_time[k] = 0;
		inversion_flag[k] = 0;
	}
	
	Mx_avg[0] = Mx[0];
	My_avg[0] = My[0];
	Mz_avg[0] = Mz[0];
	
	
	for( j = 1 ; j < NO_OF_TIME_STEPS ; j++){
	
		for( k = 0 ; k < NO_OF_ENSEMBLE; k++){			
			
			aj = aj_value(theta_en[k], phi_en[k]);
			bj = beta * aj;
			
		
				heff[0] = ( hex[0]  - (h_demag[0] - h_aniso[0] ) * Mx[k] )/(4 * M_PI * Ms);
				heff[1] = ( hex[1]  - (h_demag[1] - h_aniso[1] ) * My[k] )/(4 * M_PI * Ms);
				heff[2] = ( hex[2]  - (h_demag[2] - h_aniso[2] ) * Mz[k] )/(4 * M_PI * Ms);
			
			for(i = 0; i<3 ; i++){
				T[i] = heff[i] + ( bj - alpha * aj ) * Mf[i];
				Tp[i] = alpha * heff[i] + ( aj + alpha*bj ) * Mf[i];
	
				//printf("T[%d] is %lf \n", i, T[i]);
				//printf("%lf, %lf\n", aj, bj);
			}//end of for loop(i)
			
			
			k_theta[0] = thetaDot( theta_en[k], phi_en[k], T, Tp );
			//printf("%lf\n",k_theta[0]);
				k_phi[0] = phiDot( theta_en[k], phi_en[k], T, Tp );
				//printf("%lf\n",k_phi[0]);
			stoch_theta[0] = stochasticThetaDot( theta_en[k], phi_en[k], &seed );
			//printf("%lf\n",stoch_theta[0]);
				stoch_phi[0] = stochasticPhiDot( theta_en[k], phi_en[k], &seed );
				//printf("%lf\n",stoch_phi[0]);
			
			Mtheta = theta + dt*k_theta[0] + sqrt(dt) * stoch_theta[0];			
				Mphi = phi + dt*k_phi[0] + sqrt(dt) * stoch_phi[0];
				
			k_theta[1] = thetaDot( Mtheta, Mphi, T, Tp );
				k_phi[1] = phiDot( Mtheta, Mphi, T, Tp );
			stoch_theta[1] = stochasticThetaDot( Mtheta, Mphi, &seed );
				stoch_phi[1] = stochasticPhiDot( Mtheta, Mphi, &seed );
				
			
			theta_en[k] = theta_en[k] + dt/2.*(k_theta[0]  + k_theta[1]) + 1/2.*(stoch_theta[0] + stoch_theta[1])*sqrt(dt);
         		phi_en[k] = phi_en[k] + dt/2.*(k_phi[0]  + k_phi[1]) + 1/2.*(stoch_phi[0] + stoch_phi[1])*sqrt(dt);
		
			if(theta_en[k]<0)
			{
				theta_en[k] =-theta_en[k];	
				phi_en[k] = M_PI + phi_en[k];
			}
			if(theta_en[k]>M_PI)
			{
				theta_en[k] = 2*M_PI - theta_en[k];
				phi_en[k] = phi_en[k] + M_PI;
			}	
			
			//added
			if(phi_en[k]<0)
			{
				phi_en[k] = 2*M_PI + phi_en[k];
			}
			if(phi_en[k]> 2*M_PI)
			{
				phi_en[k] = phi_en[k] - 2*M_PI;
			}	
         		
         	
         	Mx[k] = sin(theta_en[k])*cos(phi_en[k]);
         	My[k] = sin(theta_en[k])*sin(phi_en[k]);
         	Mz[k] = cos(theta_en[k]);
         	
         	//if(abs(Mz[k]-reverse_switching_polarity)<tol_mz && inversion_flag[k]==0)
         	if(inversion_flag[k]==0 && fabs(Mz[k] - Mz_previous)>1e-7 && fabs(Mz[k])>tolerance)
			{
					
				inversion_time[k]=ti;				
				inversion_flag[k]=1;
				//temp = 1;

			}
			/*
			//if(Mz[k]<0 && dwelling_flag[k]==0)
			if(Mz[k]<0 && temp2==0)
			{
				ta[k] = ti;
				//dwelling_flag[k]=1;				
				temp2 = 1;
			}
//			if(sin(theta_en[k])<0 && abs(Mz[k])<tol_mz && dwelling_flag[k]==1)
			if(sin(theta_en[k])<0 && abs(Mz[k])<tol_mz && temp2==1)
			{
				tb[k] = ti;
				//dwelling_flag[k] = -1;
				temp2 = -1;

			}         	
         	*/
         	
         	

		}//end of forloop(k)
		
		Mx_avg[j] = 0;
        My_avg[j] = 0;
        Mz_avg[j] = 0;
        
		for(k=0; k<NO_OF_ENSEMBLE; k++)
 		{
 			Mx_avg[j]+= Mx[k];
 			My_avg[j]+= My[k];
 			Mz_avg[j]+= Mz[k];
		
 		}
 		
	   	Mx_avg[j] /= (double)NO_OF_ENSEMBLE;
      	My_avg[j] /= (double)NO_OF_ENSEMBLE;
      	Mz_avg[j] /= (double)NO_OF_ENSEMBLE;
        
		ti = ti + dt;
		
		//fprintf(fp, "%e\t%lf\t%lf\t%lf\n", ti/(4*M_PI*gama*Ms), Mx_avg, My_avg, Mz_avg);
		
		
		if(j%(NO_OF_TIME_STEPS/100)==0)
         	printf("%ld percent complete\n", (100*j/NO_OF_TIME_STEPS));
	}//end of forloop(j)	
	
	
	ti = 0;
	
	for(j=0; j<NO_OF_TIME_STEPS ; j++)
	{
		fprintf(fp, "%.4f\t%.4f\t%.4f\t%.4f\n", ti, Mx_avg[j], My_avg[j], Mz_avg[j]);
		ti = ti+dt;
	}
	storeInversionTime(inversion_time,*(filename + 1));
	free(Mx);
	free(My);
	free(Mz);
	free(inversion_flag);
	free(inversion_time);
	free(theta_en);
	free(phi_en);
	fclose(fp);
	//fclose(ts);
	//fclose(dw);
	return;
}

// Theta_dot 
double thetaDot(double theta, double phi, double *T, double *Tp){
	double K = 1/(1+alpha*alpha);
	double val =  ( -K *
			 ( 
			 	  (T[0]*sin(phi)) 
			 	- (T[1] * cos(phi)) 
			 	- (Tp[0]*cos(theta)*cos(phi)) 
			 	- (Tp[1]*cos(theta)*sin(phi)) 
			 	+ (Tp[2]*sin(theta)) 
			  ) 
			 
			);
	
	//printf("Theta_dot is %lf\n", val);
	return val;
}

// Phi_dot
double phiDot(double theta, double phi, double *T, double *Tp){
	double K = 1/(1+alpha*alpha);
	double sin_theta = sin(theta);
	double cos_theta = cos(theta);
	double sin_phi = sin(phi);
	double cos_phi = cos(phi);
	//printf("Value of theta and phi is is %.4f, %.4f\n", theta, phi);
	return ( -K*
				( 	( T[0]*cos_theta*cos_phi/sin_theta )
				   -( T[1] * cos_theta*sin_phi/sin_theta )
				   -( T[2] )	
				   +( Tp[0] * sin_phi / sin_theta )
				   -( Tp[1] * cos_phi / sin_theta )
				)
			);
}


//Stochastic Theta_dot
double stochasticThetaDot(double theta, double phi, long *seed){
	double sin_theta = sin(theta);
	double cos_theta = cos(theta);
	double sin_phi = sin(phi);
	double cos_phi = cos(phi);
	 double K = 1/(1+alpha*alpha);
	 double val =   gauss(seed)*sin_phi - gauss(seed)*cos_phi - alpha*gauss(seed)*cos_theta*cos_phi
	 			  - alpha *gauss(seed)*cos_theta*sin_phi + alpha * gauss(seed)*sin_theta;
	 val =  -K*sqrt(2*alpha*KB_TIMES_TEMPERATURE)/sqrt(gama*Ms*volume)*val*10e+6;
	 return val/(4*M_PI*Ms);
}


//stochastic Phi_dot
double stochasticPhiDot(double theta, double phi, long *seed){
	double sin_theta = sin(theta);
	double cos_theta = cos(theta);
	double sin_phi = sin(phi);
	double cos_phi = cos(phi);
	 double K = 1/(1+alpha*alpha);
	 double val = 	  gauss(seed)*cos_theta*cos_phi/sin_theta 
	 				+ gauss(seed)*cos_theta*sin_phi/sin_theta 
	 				- gauss(seed)
	 				+ alpha*gauss(seed)*sin_phi / sin_theta
	 				- alpha *cos_phi / sin_theta ;
	val = -K*sqrt(2*alpha*KB_TIMES_TEMPERATURE)/sqrt(gama*Ms*volume)*val*10e+6/sin_theta;
	val =  val/(4*M_PI*Ms);
	return val;
}

//Calculate Aj
double aj_value(double theta, double phi){
	double sin_theta = sin(theta);

	//dot product of M.Mf
	double unknownConst = 0.1;
	double dotProduct = sin_theta*cos(phi)*Mf[0] + sin_theta*sin(phi)*Mf[1] + cos(theta)*Mf[2];
	//double dotProduct =M[0]*Mf[0] + M[1]*Mf[1] + M[2]*Mf[2];
	//printf("dotProduct = %lf\n",dotProduct);
	double g = eta/(1 + unknownConst*dotProduct);
	//printf("g = %lf\n",g);
	double val = (HCROSS * g * current)/(2 * ECHARGE * Ms * volume);
	val = val/(4*M_PI*Ms);
	return val;
}

void storeInversionTime(double *inversion_time, char *filename)
{
	
	FILE *fp = fopen(filename, "w");
	int i;
	for(i=0; i<NO_OF_ENSEMBLE; i++)
	{
		fprintf(fp, "%d\t%.4f\n", i+1, inversion_time[i]);

	}
	fclose(fp);

}











