#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>


#define HCROSS 1.054e-27
#define ECHARGE 1.602e-19
#define NO_OF_TIME_STEPS 20
#define NO_OF_ENSEMBLE 10
#define KB_TIMES_TEMPERATURE 414.194e-16  //what and how??
#define TWO_PI M_PI*2
//Global Variables
static double	alpha,	gamma, current,	Ms, eta, beta , volume;
static double 	heff[3], Mf[3];
static double  	 h_aniso[3], h_demag[3], hex[3];
static double 	theta, phi;
long NO_OF_ENSEMBLE, NO_OF_TIME_STEPS;

//T initial and time steps
double ti, dt;

// global constants..
static double aj, bj;

//For generating random number;
long seed = -10000;

//Required function definitions: 

//1. Returns aj value depending on m-cap;
double Aj();

//2.
double Bj();

//3. First derivative of theta w.r.t. 
double thetaDot();

//4.
double phiDot();

//5.
double stochasticThetaDot();

//6.
double stochasticPhiDot();

/************************
 * Simulation Types: 	*
 ************************/
 

//1. 	Only H-effective term;
void Heffonly(){

}

//2.	Thermal Fluctuation only
void ThermalFluc(){

}

//3.	H-effective and  Damping 
void Heff_Damp()

//4.	H-effective , Damping and Spin Transfer Torque
void Heff_Damp_STT()

//5.	Heffective , Damping and Fluctuation
void Heff_Damp_Fluc()

//6.	Heffective , Damping , spin Transfer Torque and Fluctuation
void Heff_Damp_STT_Fluc()

//7.	H-effective , Damping, Spin transfer Torque and Field like Torque
void Heff_Damp_STT_FLT()

//8.	H-effective , Damping, Spin transfer Torque , Field like Torque and Fluctuation
void  Heff_Damp_STT_FLT_Fluc()

#include"ProjectSpin.c"

