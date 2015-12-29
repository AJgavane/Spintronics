/*
 * Simulation for Spintronics and Spin Transfer Torque
 * Code writen by Sudarshan Ghonge and Ajinkya Gavane
 * Department of Physics,
 * Birla Institute Institute of Technology and Science, Pilani. Hyderabad Campus
 *
 * 
 * The following is a program to simulate the dynamics of a soft iron layer in Magnetic Tunneling Junction.
 * The Landu-Lifhittz-Gilbert Equation has been used.
 * Refer http://en.wikipedia.org/wiki/Landau%E2%80%93Lifshitz%E2%80%93Gilbert_equation
 
 * The terms used are 
 * 1. Effective Magnetic Field
 * 2. Damping
 * 3. Spin Transfer Torque
 * 4. Field Like Toruqe
 * 5. Stochastic Field simulated using white Gaussian Noise

 * The load of solving for three quantities namely mx, my and mz -the three components of the magnetization of the soft iron magnet
 * has been reduced by assuming that none of the terms lead to dissipation and the magnitude remains constant and hence the dynamics of only the
 * angular variables has been computed.

 * The dynamics have been done on the time variable using Runge-Kutta 2 method.
 

 */


#include<stdio.h>
#include"gauss.c"
#include<time.h>
#include<string.h>
#define HCROSS 1.054e-27   //The reduced planck's constant in cgs units
#define ECHARGE 1.602e-20  //The electric charge in emu
#define NO_TIME_STEPS 30000//Number of time steps
#define NO_ENSEMBLES 100000 //Number of ensembles
#define KB_TIMES_TEMP 414.194e-16 //Value of the product of the "Kb", Boltzmann Constant and temperature , T=300 K
#define TWO_PI 2*M_PI   //Value of 2*pi
#define CLEARBUF() char ch; while(ch = getchar() != '\n' && ch!= EOF); //CLEARBUF is used to clear the input buffer. 
#define TIME_STEP 1e-13

void dump_line(FILE * fp) {
    int ch;
    while ((ch = fgetc(fp)) != EOF && ch != '\n') {
        /* null body */;
    }
}


/*

	Declaration of static variables which retain values of the parameters and the function dostuff() makes them dimensionless

*/
static double alpha, gama, current, ms, eta, beta, vol;
static double heff[3], mf[3], hext[3], han[3], hde[3];
static double theta, phi;

double ti, dt;

static double aj_var, bj_var;

#include"generateName.c"

long seed = -10000;  //Seed for generating random numbers

/*
	Returns value of aj which is instantaneously dependent on m-hat
*/
double aj(double, double);

/*
	Returns value of heff using current value of mx, my and mz
*/
void heffCalc(double, double, double);

/*
	The following two functions return the dimensionless value of the deterministic part of the dynamics of theta and phi.
	The value returned depends on the passed value of theta, phi, T and Tp.
*/

double thetaDot(double , double , double *, double *);
double phiDot(double , double , double *, double *);
/*
	The following two functions return the dimensionless value of the stochastic part of the dynamics of theta and phi.
	The value returned depends on the passed value of theta, phi and on the random number generator.
*/

double stocThetaDot(double , double , long *);
double stocPhiDot(double , double , long *);

/*
	The following stores the inversion time and dwelling time of the soft magnet
*/

void storeInversionTime(double *, char*);
void storeDwellingTime(double *, char *);


/*
	Below are various modes of the simulation 
*/

/*
 * 1. Simulates the presence of only the Effective field
 * 2. Simulates the presence of the Effective field in the presence of Damping
 * 3. Simulates the presence of the Effective field, damping and Spin Transfer Torque
 * 4. Simulates the presence of all deterministic terms: Effective field, damping, Spin Transfer Torque and Field Like Torque
 * 5. Simulates the presence of the Effective field, damping and Fluctuating field
 * 6. Simulates the presence of the Effective field, damping, Spin Transfer Torque and Fluctuating field
 * 7. Simulates the presence of the Effective field, damping, STT, FLT and Fluctuating field
 * 8. Only the Fluctuating field.
 */

//1.
void heffOnly();
//2.
void heffDamp();
//3.
void heffDampStt();
//4.
void heffDampSttFlt();
//5
void heffDampFluc();
//6
void heffDampSttFluc();
//7
void heffDampSttFltFluc();
//8
void flucOnly();

void displayPrompt();

#include"functions.c"
