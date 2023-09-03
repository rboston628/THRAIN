//**************************************************************************************
//							STAR ABSTRACT CLASS
// Star.h
//		A virtual class star for later models to inherit
//		Serves as the basis in code for other stellar models
//		Can be used in code as a generic star, without regard to particular physics
//**************************************************************************************

#ifndef STARH
#define STARH

#include "../constants.h"

class Star {
public:	
	//this is an identifier for the model of star
	//must be chosen to make UNIX accessible filename -- no spaces!
	std::string name;	
	virtual std::string graph_title() =0;
	//this is the index where shooting methods for modes should join
	std::size_t indexFit;

	// boundary conditions
	//must return stellar structure expanded in powers of x=r/R near center
	virtual void getAstarCenter(double *, int&, double g=0) =0;
	virtual void getUCenter( double*, int&) =0;
	virtual void getVgCenter(double*, int&, double g=0) =0;
	virtual void getC1Center(double*, int&) =0;
	//must return stellar structure expanded in powers of t=1-r/R near surface
	virtual void getAstarSurface(double *, int&, double g=0) =0;
	virtual void getUSurface( double*, int&) =0;
	virtual void getVgSurface(double*, int&, double g=0) =0;
	virtual void getC1Surface(double*, int&) =0;

	//distance from center of star
	virtual double rad(std::size_t) =0;
	//density, pressure, potential, and their derivatives
	virtual double rho(std::size_t) =0, drhodr(std::size_t) =0;
	virtual double   P(std::size_t) =0,   dPdr(std::size_t) =0;
	virtual double Phi(std::size_t) =0, dPhidr(std::size_t) =0;
	//interior mass to r
	virtual double  mr(std::size_t) =0;
	//Schwarzschild discriminant
	virtual double Schwarzschild_A(std::size_t, double g=0.0) =0;
	virtual double getAstar(std::size_t, double g=0.0) = 0;
	virtual double getU( std::size_t) = 0;
	virtual double getVg(std::size_t, double g=0.0) = 0;
	virtual double getC( std::size_t) = 0;
	virtual double Gamma1(std::size_t) =0;
	virtual double sound_speed2(std::size_t, double g=0.0) =0;
	
	//destructor
	virtual ~Star(){};
	//return length of star
	virtual std::size_t length() =0;
	//dimensionfull quantities of interest
	virtual double Radius() =0;
	virtual double Mass() =0;
	virtual double Gee() =0;
	virtual double light_speed2() =0;
	//print relevant values of the star in .txt and gnuplot
	virtual void writeStar(const char *const c = NULL);
	virtual void printStar(const char *const c = NULL);
	virtual void printBV(  const char *const c = NULL, double const gam1=0.0);
	virtual void printCoefficients(const char *const c = NULL, double const gam1=0.0);
	virtual double SSR();
	
	//allow Modes to access private members of Star
	friend class ModeDriver;
};

#endif