/********************************************************************************
 *	Project Name: Spintronics and Spin Transfer Torque.							*
 *	Written By : Ajinkya S. Gavane , Sudarshan M. Ghonge						*
 *	----------------------------------------------------------------------------*	
 *	Some related Theory and constant and variables used in the code.			*
 *																				*
 * NOTE: 'x' means vector product, '*'Multiplication and '.' means Scalar		*
 * product of vectors. Thus terms associated with 'x' and '.' are vectors 		*
 * similarly terms associated with '*' are scalar. 								*					
 *																				*
 *  Landou-Lifshift-Gilbert(LLG) Equation: 										*
 *  dM/dt = ( -gamma/(1 + alpha^2) )*(M x h_eff) 								*
 *			- ( gamma/ (1 + alpha^2) )*aj*( M_cap x (F x Mf) ) 					*
 *			- ( gamma * alpha/(a + alpha ^2) )( M_cap x (M x h_eff) )			*
 *			- ( gamma/ (1 + alpha^2) )*bj*( M x Mf ) 							*
 *			+ ( gamma * alpha/(a + alpha ^2) ) * aj * ( M x Mf) 				*
 *			- ( gamma * alpha/(a + alpha ^2) ) * bj * ( M_cap x ( M x Mf ) )	*
 *         ---------------------------------------------                        *
 *																				*
 *	Landou-Lifshift Equation: 													*
 *	dM/dt = - gamma * ( M x h_eff ) 			//precision term				*
 *			+ gamma * ( M x h_thermal)			//stochastic term				*
 *			+ alpha * ( M_cap x (dM/dt) ) 		// damping or dissipition		*
 * 			- gamma * aj * ( M x ( m * Mf) )	//Spin Transfer torque(STT)		*
 *			- gamma * bj * ( M x Mf )			//Field like torque(FLT)		*
 *     		--------------------------------------------						*
 *																				*
 *	M_cap = (M / |Ms| )		//Ms is saturation magnetization, M is a vector		*
 *	h_eff = - dovE(M) / dovM	//partial differentiation of E(M) w.r.t M		*
 *	^ Effective field which depends on free energy density						*
 * 	E(M) = - (1/2)*H_an*Ms*(M_cap.e_cap)^2  -  h_ext . M 						*							
 *			+ 2 * PI * Ms * Ms * (M_cap . z_cap)^2		  //Z_cap is direction	* 
 ********************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"ProjectSpin.c"


int main(int argc, char *argv[]){
	//Required Variable constants
	//newly added
	double toleranceDim,tempDim,time_step_sizeDim;
	//done
	double alphaDim, gamaDim, currentDim, MsDim, etaDim ,betaDim;
	long double volumeDim;
	double hexDim[3], MfDim[3],h_anisoDim[3],h_demagDim[3];
	double thetaDim, phiDim;
	int choice;
	long No_of_Ensemble;
	double time_of_sim;
	double constants[21];
	//printf("Value of time steps is %d and that of no of ensembles is %d\n", NO_OF_TIME_STEPS, NO_OF_ENSEMBLE);
	//Reading constants value stored in the file called Constants.d
	FILE *fp = fopen("Constants.d","r");
	if(fp == NULL){
		printf("Get the file containing the Constants, named as -Constants.d- !!");	
		return EXIT_FAILURE;	
	} 
	
	
		fscanf(fp,"%*s %lf",&alphaDim);
		fscanf(fp,"%*s %lf",&gamaDim);
		fscanf(fp,"%*s %lf",&currentDim);
		fscanf(fp,"%*s %lf",&MsDim);
		fscanf(fp,"%*s %lf",&etaDim);
		fscanf(fp,"%*s %lf",&betaDim);
		fscanf(fp,"%*s %Lf",&volumeDim);
		fscanf(fp,"%*s %lf",&hexDim[0]);
		fscanf(fp,"%*s %lf",&hexDim[1]);
		fscanf(fp,"%*s %lf",&hexDim[2]);
		fscanf(fp,"%*s %lf",&h_anisoDim[0]);
		fscanf(fp,"%*s %lf",&h_anisoDim[1]);
		fscanf(fp,"%*s %lf",&h_anisoDim[2]);
		fscanf(fp,"%*s %lf",&h_demagDim[0]);
		fscanf(fp,"%*s %lf",&h_demagDim[1]);
		fscanf(fp,"%*s %lf",&h_demagDim[2]);
		fscanf(fp,"%*s %lf",&MfDim[0]);
		fscanf(fp,"%*s %lf",&MfDim[1]);
		fscanf(fp,"%*s %lf",&MfDim[2]);
		fscanf(fp,"%*s %lf",&thetaDim);
		fscanf(fp,"%*s %lf",&phiDim);		
		fscanf(fp,"%*s %ld",&No_of_Ensemble);
		fscanf(fp,"%*s %lf",&time_of_sim);
		fscanf(fp,"%*s %d",&choice);
		fscanf(fp,"%*s %lf",&tempDim);
		fscanf(fp,"%*s %lf",&toleranceDim);
		fscanf(fp,"%*s %lf",&time_step_sizeDim);
	fclose(fp);
	
	//Go to the calculation file, where main code for the project exis..
	
	
	/* To check if the file is being read...*
	printf("alpha = %lf\n",alphaDim);
	printf("gammaDim = %lf\n",gamaDim);
	printf("currentDim = %lf\n",currentDim);
	printf("MsDim = %lf\n",MsDim);
	printf("etaDim = %lf\n",etaDim);
	printf("betaDim = %lf\n",betaDim);
	printf("volumeDim = %Lf\n",volumeDim);
	printf("hexDim_x = %lf\n",hexDim[0]);
	printf("hexDim_y = %lf\n",hexDim[1]);
	printf("hexDim_z = %lf\n",hexDim[2]);
	printf("h_anisoDim_x = %lf\n",h_anisoDim[0]);
	printf("h_anisoDim_y = %lf\n",h_anisoDim[1]);
	printf("h_anisoDim_z = %lf\n",h_anisoDim[2]);
	printf("h_demagDim_x = %lf\n",h_demagDim[0]);
	printf("h_demagDim_y = %lf\n",h_demagDim[1]);
	printf("h_demagDim_z = %lf\n",h_demagDim[2]);
	printf("MfDim_x = %lf\n",MfDim[0]);
	printf("MfDim_y = %lf\n",MfDim[1]);
	printf("MfDim_z = %lf\n",MfDim[2]);
	printf("thetaDim = %lf\n",thetaDim);
	printf("phiDim = %lf\n",phiDim);
	printf("choice = %d\n",choice);
	printf("No_of_Ensemble = %ld\n",No_of_Ensemble);
	printf("time_step = %ld\n",time_step);
	/*
	*/ 
	
	ProjectSpin(alphaDim, gamaDim, currentDim, MsDim, etaDim ,betaDim, volumeDim, hexDim,
				 h_anisoDim, h_demagDim, MfDim, thetaDim, phiDim ,No_of_Ensemble, time_of_sim, choice,
				 toleranceDim,tempDim,time_step_sizeDim);
	
	return EXIT_SUCCESS;
}
