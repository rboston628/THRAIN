//**************************************************************************************
//							The MODE Object
//	Mode.h
//  	Inhereits form ModeBase, contained in ModeDriver.g
//		Equation agnostic 
//			-- all information specific to physics supplied by the driver
// 		Capable of nonradial, Cowling, 1PN modes by supplying different drivers
//  Reece Boston Mar 24, 2022
//**************************************************************************************

#ifndef MODEH
#define MODEH

#include "../constants.h"
#include "../STARS/Star.h"
#include "../MODES/ModeDriver.h"

template <size_t numvar> 
class Mode : public ModeBase {
public:
	static const unsigned int  num_var = numvar;
	//retrieve quantum numbers of pulsation
	int  modeOrder() const override {return k;};
	void modeNumbers(int& K, int& L, int& M) const override {K=k; L=l; M=m;};
	//constructors: differ only in initial choice of frequency
	//initialize from quantum numbers and a background star
	Mode<numvar>(int k, int l, int m, ModeDriver*);
	//initialize from frequency guess and a background star
	Mode<numvar>(double omega2, int l, int m, ModeDriver*);
	//initialize from brackets on frequency and a background star
	Mode<numvar>(double omegMin, double omegMax, int l, int m, ModeDriver*);
	//owns raw arrays; copying would double-free
	Mode<numvar>(Mode<numvar> const&) = delete;
	Mode<numvar>& operator=(Mode<numvar> const&) = delete;
	//deconstructor
	virtual ~Mode();
	
	//ways to access the frequency
	double getOmega2() const override;
	double getFreq() const override;
	double getPeriod() const override;
	double SSR() const override;
	double tidal_overlap() const override;

	//methods to access mode functions
	double getRad(std::size_t x) const override { return rad[x];};
	double getY(int i, std::size_t x) const override { return pow(rad[x],l-2)*y[x][i];};
	double getYtilde(int i, std::size_t x) const override { return y[x][i];};
	
protected:
	//the equilibrium star
	Star *const star;
	//the mode driver specifying equations and boundaries
	ModeDriver *const driver;
	void basic_setup();	//handles the normal setup common to all constructors
	int k = 0,l,m;		//the mode numbers
	double cee2 = 0.0, Gee = 0.0;	//for specifying units
	double Gamma1 = 0.0;//adiabatic exponent for perturbations
	std::size_t len = 0, len_star = 0;	//length of mode arrays, and length of background model arrays
	double omega2 = 0.0;//frequency squared
	double *rad = nullptr;	//the radial coordinate, x = r/R
	//perturbation quantities
	double **y = nullptr;	//the eigenfunctions y_1, y_2, ..., y_N
	double yCenter[numvar], ySurface[numvar];//central and surface values of y_1,...y_N

	bool converged = false;	//a flag for indicating if calculation converged -- no longer needed
	void convergeBisect(double);		//find frequency using a bisection search
	void convergeNewton(double, int);	//find frequency using a Newton method

	double calculateWronskian(double w2);
	void linearMatch(double w2, double y0[numvar], double ys[numvar]);
	std::function<double(double)> wronskian {};
	int verifyMode();
	void converge();
	//
	void RK4step(std::size_t, double w2, double dx, double yin[numvar], double yout[numvar]);
	void RK4out(std::size_t, double, double y0[numvar]);
	void RK4in( std::size_t, double, double ys[numvar]);
	
	double boundaryMatrix[numvar][numvar]; //the BCs specified in matrix form
	int    indexOrder[numvar];             //the order in which to multipy factors to match
	
	std::size_t xfit = 0;    //the coordinate where inner/outer solutions are fit
	double omeg2freq = 0.0;  //converstion from dimensionless frequency to angular frequency
	double R = 0.0;          //the radius of the star

public:
	//file output methods to write and view plots of mode
	void writeMode(std::string const& c = "") const override;
};

//template classes in C++ cannot be split into multiple files, so we must include Mode.pp
#include "Mode.cpp"

#endif