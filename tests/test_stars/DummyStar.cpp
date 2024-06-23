//**************************************************************************************
//							DUMMY STAR
// DummyStar.cpp
//		A simple object for use in testing.
//		It is only there to be used when something requries a Star object
//		THIS IS NOT PART OF THRAIN -- ONLY FOR CERTAIN TESTS
//**************************************************************************************


#ifndef DUMMYSTARCLASS
#define DUMMYSTARCLASS
#include "DummyStar.h"

//initalize polytrope from length
DummyStar::DummyStar(std::size_t L) : len(L) {
	indexFit = len/4;
	name = "dummy_star";
}

DummyStar::~DummyStar(){}

//Here we define functions to access radius, pressure, etc.
double DummyStar::rad(std::size_t X)    {return double(X)/double(len-1);}
double DummyStar::mr(std::size_t X)     {return rad(X);}
double DummyStar::rho(std::size_t X)    {return 1.0;}
double DummyStar::drhodr(std::size_t X) {return 0.0;}
double DummyStar::P(std::size_t X)      {return 1.0;}
double DummyStar::dPdr(std::size_t X)   {return 0.0;}
double DummyStar::Phi(std::size_t X)    {return 0.0;}
double DummyStar::dPhidr(std::size_t X) {return 0.0;}

double DummyStar::Schwarzschild_A(std::size_t X, double GamPert){return 0.0;}
double DummyStar::getAstar(std::size_t X, double GamPert){return 0.0;}
double DummyStar::getU(std::size_t X){return 3.0;}
double DummyStar::getVg(std::size_t X, double GamPert){return 0.0;}
double DummyStar::getC(std::size_t X){return 1.0;}
double DummyStar::Gamma1(std::size_t X){return 1.0;}

double DummyStar::sound_speed2(std::size_t X, double GamPert){return 1.0;}

double DummyStar::Radius(){return 1.0;}	//total radius
double DummyStar::Mass(){  return 1.0;} //total mass
double DummyStar::Gee(){   return 1.0;} //constant G


// **************************  CENTRAL BOUNDARY **************************************
// the following provide coefficients for central expansions of A*, Vg, U, c1 in trms of x=r/R
//	Note that only even powers are needed
//	Up to order 2N, require terms 0, 2, 4, ... , 2N
//	The expansion coefficients must be in powers of x=r/R, not in powers of xi = r/Rn
//  Note that because of extremely simple structure, coefficiens available to arbitrary precision

void DummyStar::getAstarCenter(double *Ac, int& maxPow, double g){
	for (int i=0; i<=maxPow/2; i++) Ac[i] = 0.0;
}

void DummyStar::getVgCenter(double *Vc, int& maxPow, double g){
	getAstarCenter(Vc, maxPow, g);
}

void DummyStar::getUCenter(double *Uc, int& maxPow){
	if(maxPow>=0) Uc[0] = 3.0;
	for(int k=1; k<= maxPow/2; k++){
		Uc[k] = 0.0;
	}
}

void DummyStar::getC1Center(double *cc, int& maxPow){
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

void DummyStar::getAstarSurface(double *As, int& maxPow, double g){
	for (int i=0; i<=maxPow; i++) As[i] = 0.0;
}

void DummyStar::getVgSurface(double *Vs, int& maxPow, double g){
	getAstarSurface(Vs, maxPow, g);
}

void DummyStar::getUSurface(double *Us, int& maxPow){
	if(maxPow>=0) Us[0] = 3.0;
	for(int k=1; k<= maxPow; k++){
		Us[k] = 0.0;
	}
}

void DummyStar::getC1Surface(double *cs, int& maxPow){
	if(maxPow>=0) cs[0] = 1.0;
	for(int k=1; k<= maxPow; k++){
		cs[k] = 0.0;
	}
}

#endif