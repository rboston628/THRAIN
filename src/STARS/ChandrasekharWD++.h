//**************************************************************************************
//							CHANDRASEKHAR WHITE DWARF
// ChandrasekharWD++.cpp
// 		This is a model of a WD based on the equation (13)  of Chandrasekhar 1935
//			based on a cold degenerate electron gas equation of state
//			Chandrasekhar defines a variable y in equation (12)
//		This model assumes T=0, and ignores Coulombic and other effects
//		The surface is not treated in any special way
//		Updated to include a chemical profile, indicated by mu_electron
// Reece Boston, Sep 03, 2023
// **************************************************************************************/

#ifndef ChandrasekharWDH
#define ChandrasekharWDH

#include "Star.h"
#include "../../lib/chandra.h"
#include <functional>

namespace Calculation {
	struct InputData;
}

class ChandrasekharWD : public Star {
public:

	static int read_star_input(FILE* input_file, Calculation::InputData&);

	std::string graph_title() override final {
		return strmakef("Chandrasekhar WD with y_0=%1.2f", Y0);
	}
	
	//the constructors
	typedef std::function<void(double const, double const, double&, double&)> ChemicalGrad;
	ChandrasekharWD(double, std::size_t,               ChemicalGrad);
	ChandrasekharWD(double, std::size_t, const double, ChemicalGrad);
	ChandrasekharWD(double, std::size_t, double const A0, double const B0);

	virtual ~ChandrasekharWD();   //destructor
	std::size_t length() override final {return len;}
	//these three functions specify units
	double Radius() override final;	//total radius
	double Mass() override final;//total mass
	double Gee() override final; //{return G_CGS;};
	
	double rad(std::size_t) override final;
	double rho(std::size_t) override final;
	double   P(std::size_t) override final;
	double Phi(std::size_t) override final;
	double  mr(std::size_t) override final;
	// derivs
	double drhodr(std::size_t) override final;
	double   dPdr(std::size_t) override final;
	double dPhidr(std::size_t) override final;
	
	double Schwarzschild_A(std::size_t, double GamPert=0.0) override final;
	double getAstar(std::size_t, double GamPert=0.0) override final;
	double getVg(std::size_t, double GamPert=0.0) override final;
	double getU(std::size_t) override final;
	double getC(std::size_t) override final;
	double Gamma1(std::size_t) override final;
	double sound_speed2(std::size_t, double GamPert=0.0) override final;
	double Ledoux(std::size_t, double GamPert=0.0);
	
private:
	void basic_setup();
	void init_arrays();
	std::size_t len;
	double Y0;		// central value of y, y^2=1+x^2
	double X0, X02, Y02;
	double Rn;		// radius scale
	double A0, B0;  // pressure and density scales
	//lane-emden solution functions
	// @xi	normalized radius
	// @x   the relativity factor x = pF/mc
	// @y   Chandrasekhar's y, y^2=1+x^2
	// @z   derivative (dy/dxi)  note: dx/dxi = (dy/dxi)/x
	// @f   the factor f(X)
	enum VarName {xi=0, y, z, x, f, numvar};
	double **Y;
	double *mass;
	double *x3; // holds x^3, to avoid repeated calls to pow
	double set_mass(double const y[numvar], double const x);
		
	// the chemical profile
	ChemicalGrad chemical_gradient;

	//to handle chemical profile
	double* mue;   //mean atomic mass per electron
	double* dmue;  //derivative of above

	//integrate using basic RK4
	void centerInit(double ycenter[numvar]);
	void RK4step(double dx, double yin[numvar], double yout[numvar]);
	enum class SurfaceBehavior : bool {CONTINUE_FULL_LENGTH=false, STOP_AT_ZERO=true};
	double RK4integrate(double, SurfaceBehavior);
		
	//methods for handling the BCs
	double yc[4], xc[3], fc[2];	//series coefficients of y,x,f near center
	void setupCenter();		//prepare values near center
	void setupSurface();	//prepare values near surface
public:
	//methods to find central, surface power series expansions of key variables in pulsation
	void getAstarCenter(double *, int&, double g=0) override final;
	void getUCenter(double*, int&) override final;
	void getVgCenter(double*, int&, double g=0) override final;
	void getC1Center(double*, int&) override final;
	void getAstarSurface(double *, int&, double g=0) override final;
	void getUSurface(double*, int&) override final;
	void getVgSurface(double*, int&, double g=0) override final;
	void getC1Surface(double*, int&) override final;

	//a particular output generation for this model of white dwarf
	void writeStar(std::string const& calcname = "") override final;
	void printDeg(std::string const& calcname);
	void printChem(std::string const& calcname);
};

namespace Chandrasekhar {

	struct constant_mu {
		double mu0;
		void operator()(double const, double const, double& mu, double& dmu){
			mu  = mu0;
			dmu = 0.0;
		}
	};

	struct sigmoidal_in_logf {
		double k, F0, mu0, mu1;
		void operator()(double const x, double const dydxi, double& mu, double& dmu){
			double F  = Chandrasekhar::factor_f(x);
			double z = -log(F/F0);
			double zeec = -log(1e-5);     // TODO why 1e-5?
			double EXP = exp(k*(z-zeec)); // the expoential
			mu =  mu1 + (mu0-mu1)/(1.+EXP);
			//   dmudx  = k f'(x)/f(x) * EXP/(1+EXP)^2           = 8k x^4/y/f       * EXP/(1+EXP)^2
			//   dmudxi = dmudx * dx/dxi = dmudx * (y/x * dydxi) = 8k x^3/f * dydxi * EXP/(1+EXP)^2
			dmu = 8.*k*pow(x,3)/F * (mu0-mu1)*EXP*pow(1.+EXP,-2) * dydxi;
		}
	};

}

#endif