#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include"gauss.c"

#define HCROSS 1.054e-27
#define ECHARGE 1.602e-20		//in cgs units
#define KB 1.3806466e-16  //KB is Boltzman constant multiplied by room temperature of 300K

//Global Variables
	//newly added
double temp, tolerance;
	//done

static double	alpha,	gama, current,	Ms, eta, beta, KB_TIMES_TEMPERATURE,time_step_size ;
static double volume;
static double 	heff[3], Mf[3], h_aniso[3], h_demag[3], hex[3];
static double 	theta, phi;
static long NO_OF_ENSEMBLE, NO_OF_TIME_STEPS;

//T initial and time steps
double ti, dt;

#include"generateName.c"

// global constants..
static double aj, bj;

//For generating random number;
long seed = -10000;

//Required function definitions: 

//1. Returns aj value depending on m-cap;
double aj_value(double , double );

//2.
double Bj();

//3. First derivative of theta w.r.t. 
double thetaDot(double , double , double * , double *);

//4.
double phiDot(double , double , double *, double *);

//5.
double stochasticThetaDot(double , double , long *);

//6.
double stochasticPhiDot(double , double , long *);

/************************
 * Simulation Types: 	*
 ************************/
 

//1. 	Only H-effective term;
void Heffonly();

//2.	Thermal Fluctuation only
void ThermalFluc();

//3.	H-effective and  Damping 
void Heff_Damp();

//4.	H-effective , Damping and Spin Transfer Torque
void Heff_Damp_STT();

//5.	Heffective , Damping and Fluctuation
void Heff_Damp_Fluc();

//6.	Heffective , Damping , spin Transfer Torque and Fluctuation
void Heff_Damp_STT_Fluc();

//7.	H-effective , Damping, Spin transfer Torque and Field like Torque
void Heff_Damp_STT_FLT();

//8.	H-effective , Damping, Spin transfer Torque , Field like Torque and Fluctuation
void  Heff_Damp_STT_FLT_Fluc();

void storeInversionTime(double *, char *);
//#include"ProjectSpin.c"

