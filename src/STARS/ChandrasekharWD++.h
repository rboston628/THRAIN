//**************************************************************************************
//							CHANDRASEKHAR WHITE DWARF
// ChandrasekharWD++.cpp
// 		This is a model of a WD based on the equation (13)  of Chandrasekhar 1935
//			based on a cold degenerate electron gas equation of state
//			Chandrasekhar defines a variable y in equation (12)
//		This model assumes T=0, and ignores Coulombic and other effects
//		The surface is not treated in any special way
//		Updated to include a chemical profile, indicated by mu_electron
// Reece Boston, Mar 24, 2022
// **************************************************************************************/

#ifndef ChandrasekharWDH
#define ChandrasekharWDH

#include "Star.h"
#include "../../lib/chandra.h"


class ChandrasekharWD : public Star {
public:

	std::string graph_title() override {
		return strmakef("Chandrasekhar WD with y_0=%1.2f", Y0);
	}
	
	//the constructors
	ChandrasekharWD(double, std::size_t,               double MU0, double K, double AC, double AS);
	ChandrasekharWD(double, std::size_t, const double, double MU0, double K, double AC, double AS);
	ChandrasekharWD(double, std::size_t, const double A0, const double B0);
	virtual ~ChandrasekharWD();   //destructor
	std::size_t length() override {return len;}
	//these three functions specify units
	double Radius() override;	// total radius
	double Mass() override;     // total mass
	double Gee() override;      // G in units
	
	double rad(std::size_t) override;
	double rho(std::size_t) override, drhodr(std::size_t) override;
	double   P(std::size_t) override,   dPdr(std::size_t) override;
	double Phi(std::size_t) override, dPhidr(std::size_t) override;
	double mr(std::size_t) override;
	
	double Schwarzschild_A(std::size_t, double GamPert=0.0) override;
	double getAstar(std::size_t, double GamPert=0.0) override;
	double getVg(std::size_t, double GamPert=0.0) override;
	double getU(std::size_t) override;
	double getC(std::size_t) override;
	double Gamma1(std::size_t) override;
	double sound_speed2(std::size_t, double GamPert=0.0) override;
	double Ledoux(std::size_t, double GamPert=0.0);
	
private:
	void basic_setup();
 	void init_arrays();
	double Y0;		// central value of y, y^2=1+x^2
	double X0, X02, Y02;
	double Rn;		// radius scale
	double A0, B0;  // pressure and density scales
	std::size_t len;
	double dx;
	//lane-emden solution functions
	double *xi;	//normalized radius
	double *x;  //the relativity factor x = pF/mc
	double *y;	//Chandrasekhar's y, y^2=1+x^2
	double *z;	//derivative (dy/dxi)  note: dx/dxi = (dy/dxi) * y/x
	double *mass;
	double *f;
		
	//parameters of the chemical profile
	double mu0, k, acore, aswap;
	void chemical_gradient(const double, const double, double&, double&);
	
	
	//to handle chemical profile
	double* mue;   //mean atomic mass per electron
	double* dmue;  //derivative of above
	//integrate using basic RK4
	double RK4integrate(const std::size_t, double);
	std::size_t RK4integrate(const std::size_t, double, int);
		
	//methods for handling the BCs
	double yc[4], xc[3], fc[2];	//series coefficients of y,x,f near center
	void setupCenter();		//prepare values near center
	void setupSurface();	//prepare values near surface
public:
	//methods to find central, surface power series expansions of key variables in pulsation
	void getAstarCenter(double *, int&, double g=0) override;
	void getUCenter(double*, int&) override;
	void getVgCenter(double*, int&, double g=0) override;
	void getC1Center(double*, int&) override;
	void getAstarSurface(double *, int&, double g=0) override;
	void getUSurface(double*, int&) override;
	void getVgSurface(double*, int&, double g=0) override;
	void getC1Surface(double*, int&) override;

	//a particular output generation for this model of white dwarf
	void writeStar(char const *const c=NULL) override;
	void printDeg(char const *const c);
 	void printChem(char const *const c);
};

#endif