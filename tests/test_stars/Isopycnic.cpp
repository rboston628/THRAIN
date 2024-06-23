//**************************************************************************************
//							ISOPYCNIC STAR
// Isopycnic.cpp
//		A simple stellar object in Newtonian physics of uniform density
//				rho = constant
//		Thi exploits the known analytic solution to Lane-Emden for n=0
//				y = 1 - x^2/6
//		Used as stellsr background for mode testing, to minimize background errors
//		THIS IS NOT PART OF THRAIN -- ONLY FOR CERTAIN TESTS WITH NONRADIAL MODES
//**************************************************************************************


#ifndef ISOPYCNICCLASS
#define ISOPYCNICCLASS
#include "Isopycnic.h"

//initalize polytrope from length
Isopycnic::Isopycnic(std::size_t L)
	: len(L), GG(1.0), Gamma(5./3.), 
	x2r( sqrt(6.0/(4.*m_pi)) ),
	x2m( 4.*m_pi/3. * pow(x2r, 3) )
{
	name = "isopycnic";
	indexFit = len/2;
	for(std::size_t X=1; X<len; X++){
		//as we scan through x,y,z, set matching point where y[X] = 0.5
		if(y(X-1)>0.5 & y(X+1)<=0.5) indexFit = X;
	}
	indexFit /= 2;
}

Isopycnic::~Isopycnic(){}

// normalized radius r/R
double Isopycnic::x(std::size_t X){
	return double(X)/double(len - 1);
}

// solution in terms of normalied radius
double Isopycnic::y(std::size_t X){
	return 1.0 - x(X)*x(X);
}

// derivative of y with respect to normalized radius
double Isopycnic::z(std::size_t X){
	return -2.0 * x(X);
}

//Here we define functions to access radius, pressure, etc.
double Isopycnic::rad(std::size_t X){
	return x2r * x(X);
}
double Isopycnic::rho(std::size_t X){
	return 1.0;
}
double Isopycnic::drhodr(std::size_t X){
	return 0.0;
}
inline double Isopycnic::P(std::size_t X){
	return y(X);
}
double Isopycnic::dPdr(std::size_t X){
	return z(X)/x2r;
}
double Isopycnic::Phi(std::size_t X){
	return 1.0 -y(X);
}
 double Isopycnic::dPhidr(std::size_t X){
	return -z(X)/x2r;
}
double Isopycnic::mr(std::size_t X){
	return x2m * pow(x(X), 3);
}

double Isopycnic::Schwarzschild_A(std::size_t X, double GamPert){
	if(GamPert == 0.0) return -0.6*z(X)/(y(X)*x2r);
	else               return -z(X)/(y(X)*x2r)/GamPert;
}

double Isopycnic::getAstar(std::size_t X, double GamPert){
	return -getVg(X, GamPert);
}

double Isopycnic::getU(std::size_t X){
	return 3.0;
}
double Isopycnic::getVg(std::size_t X, double GamPert){
	double gam1 = (GamPert == 0.0 ? Gamma1(0) : GamPert);
	double x2 = x(X) * x(X);
	return 2.0 * x2/(1.-x2)/gam1;
}
double Isopycnic::getC(std::size_t X){
	return 1.0;
}
double Isopycnic::Gamma1(std::size_t X){
	return Gamma;
}

double Isopycnic::sound_speed2(std::size_t X, double GamPert){
	if(GamPert == 0.0) return Gamma  *y(X);
	else               return GamPert*y(X);
}

double Isopycnic::Radius(){return x2r;}	//total radius
double Isopycnic::Mass(){ return x2m;}//total mass


// **************************  CENTRAL BOUNDARY **************************************
// the following provide coefficients for central expansions of A*, Vg, U, c1 in trms of x=r/R
//	Note that only even powers are needed
//	Up to order 2N, require terms 0, 2, 4, ... , 2N
//	The expansion coefficients must be in powers of x=r/R, not in powers of xi = r/Rn
//  Note that because of extremely simple structure, coefficiens available to arbitrary precision
void Isopycnic::setupCenter(){}

void Isopycnic::getAstarCenter(double *Ac, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	if(maxPow>=0) Ac[0] = 0.0;
	for(int k=1; k<=maxPow/2; k++){
		Ac[k] = -2./Gam1;
	}
}

void Isopycnic::getVgCenter(double *Vc, int& maxPow, double g){
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	if(maxPow>=0) Vc[0] = 0.0;
	for(int k=1; k<=maxPow/2; k++){
		Vc[k] = 2./Gam1;
	}
}

void Isopycnic::getUCenter(double *Uc, int& maxPow){
	if(maxPow>=0) Uc[0] = 3.0;
	for(int k=1; k<= maxPow/2; k++){
		Uc[k] = 0.0;
	}
}

void Isopycnic::getC1Center(double *cc, int& maxPow){
	if(maxPow>=0) cc[0] = 1.0;
	for(int k=1; k<= maxPow/2; k++){
		cc[k] = 0.0;
	}
}


// **************************  SURFACE BOUNDARY  ***************************************
// the following provide coefficients for surface expansions of A*, Vg, U, c1 in terms of t=1-r/R
//	Note that A*, Vg require a power -1
//	Note that up to order N, we only require:
//		A*, Vg, c1 up to order N-1
//		U          up to order N
//  Note that because of extremely simple structure, coefficiens available to arbitrary precision
void Isopycnic::setupSurface(){}

void Isopycnic::getAstarSurface(double *As, int& maxPow, double g){
	//we make use of the fact that A* and Vg are simply related in polytropes
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	int O=1;
	if(maxPow>= -1) As[O-1] = -1./Gam1;
	if(maxPow>=  0) As[O  ] = 1.5/Gam1;
	for(int k=1; k<= maxPow; k++){
		As[O+k] = -pow(2, -k-1)/Gam1;
	}
}

void Isopycnic::getVgSurface(double *Vs, int& maxPow, double g){
	// coefficients of Vg must extend up to maxPow-1
	double Gam1 = (g==0.0 ? Gamma1(0) : g);
	int O=1;
	if(maxPow>= -1) Vs[O-1] = 1./Gam1;
	if(maxPow>=  0) Vs[O  ] = -1.5/Gam1;
	for(int k=1; k<= maxPow; k++){
		Vs[O+k] = pow(2, -k-1)/Gam1;
	}
}

void Isopycnic::getUSurface(double *Us, int& maxPow){
	// coefficients of U must extend up to order maxPow
	if(maxPow>=0) Us[0]  = 3.0;
	for(int k=1; k<= maxPow; k++){
		Us[k] = 0.0;
	}
}

void Isopycnic::getC1Surface(double *cs, int& maxPow){
	// coefficients of c1 are only needed up to order maxPow-1
	if(maxPow>=0) cs[0]  = 1.0;
	for(int k=1; k<= maxPow; k++){
		cs[k] = 0.0;
	}
}

#endif