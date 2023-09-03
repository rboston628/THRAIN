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
	ChandrasekharWD(double, int,               double MU0, double K, double AC, double AS);
	ChandrasekharWD(double, int, const double, double MU0, double K, double AC, double AS);
	virtual ~ChandrasekharWD();   //destructor
	int length() override {return len;}
	//these three functions specify units
	double Radius() override;	//total radius
	double Mass() override;//total mass
	double Gee() override; //{return G_CGS;};
	//in Newtonian, light speed is infinity...
	virtual double light_speed2() override {
		return std::numeric_limits<double>::infinity();
	}
	
	double rad(int) override;
	double rho(int) override, drhodr(int) override;
	double   P(int) override,   dPdr(int) override;
	double Phi(int) override, dPhidr(int) override;
	double mr(int) override;
	
	double Schwarzschild_A(int, double GamPert=0.0) override;
	double getAstar(int, double GamPert=0.0) override;
	double getVg(int, double GamPert=0.0) override;
	double getU(int) override;
	double getC(int) override;
	double Gamma1(int) override;
	double sound_speed2(int, double GamPert=0.0) override;
	double Ledoux(int, double GamPert=0.0);
	
private:
	double Y0;		// central value of y, y^2=1+x^2
	double X0, X02, Y02;
	//double A0;		//pressure scale
	//double B0;		//density scale
	double Rn;		//radius scale
	int len;
	double dx;
	//lane-emden solution functions
	double *xi;	//normalized radius
	double *x;  //the relativity factor x = pF/mc
	double *y;	//Chandrasekhar's y, y^2=1+x^2
	double *z;	//derivative (dy/dxi)  note: dx/dxi = (dy/dxi)/x
	double *mass;
	double *f;
		
	//parameters of the chemical profile
	double mu0, k, acore, aswap;
	void chemical_gradient(const double, const double, double&, double&);
	
	
	//to handle chemical profile
	double* mue;   //mean atomic mass per electron
	double* dmue;  //derivative of above
	//integrate using basic RK4
	double RK4integrate(const int, double);
	int RK4integrate(const int, double, int);
	
	//the T=0 Fermi function
	//double factor_f(double x);
	
	//methods for handling the BCs
	double yc[4], xc[3], fc[2];	//series coefficients of y,x,f near center
	void setupCenter();		//prepare values near center
	void setupSurface();	//prepare values near surface
public:
	//methods to find central, surfae power series expansions of key variables in pulsation
	void getAstarCenter(double *, int&, double g=0) override;
	void getUCenter(double*, int&) override;
	void getVgCenter(double*, int&, double g=0) override;
	void getC1Center(double*, int&) override;
	void getAstarSurface(double *, int&, double g=0) override;
	void getUSurface(double*, int&) override;
	void getVgSurface(double*, int&, double g=0) override;
	void getC1Surface(double*, int&) override;

	//a particular output generation for this model of white dwarf
	void writeStar(const char *const c=NULL) override;
};

#endif