//**************************************************************************************
//							POLYTROPIC STAR
// Polytrope.h
//		A simple stellar object in Newtonian physics obeying everywhere a polytropic
//			equation of state of the form 
//				P~rho^Gamma
//		Solves Lane-Emden equation, y=theta, x=xi, z=dy/dxi in more usual notation
//		See Hansen & Kawaler Chapter 7 for further information
//**************************************************************************************

#ifndef POLYTROPEH
#define POLYTROPEH

#include "Star.h"

namespace Calculation {
	struct InputData;
}

class Polytrope : public Star {
public:

	static int read_star_input(FILE* input_file, Calculation::InputData&);

	std::string graph_title() override {
		return strmakef("polytrope n=%1.2f", n);
	}
	
	Polytrope(double, double, double, std::size_t);	//constructor, M, R, n and length
	Polytrope(double, std::size_t);					//constructor, n and length
	Polytrope(double, std::size_t, const double);	//constructor, n and length, dx
	virtual ~Polytrope();  //destructor
	std::size_t length() override {return len;}
	//these three functions specify units
	double Radius() override;	//total radius
	double Mass() override;		//total mass
	double Gee() override {return GG;};
	
	//these return value of indicated variable -- used in testing
	double getX(std::size_t X){     return Y[X][x];}
	double getY(std::size_t X){     return Y[X][y];}
	double getYderiv(std::size_t X){return Y[X][z];}
		
	double rad(std::size_t) override;
	double rho(std::size_t) override, drhodr(std::size_t) override;
	double   P(std::size_t) override,   dPdr(std::size_t) override;
	double Phi(std::size_t) override, dPhidr(std::size_t) override;
	double  mr(std::size_t) override;
	
	double Schwarzschild_A(std::size_t, double GamPert=0.0) override;
	double getAstar(std::size_t, double GamPert=0.0) override;
	double getVg(std::size_t, double GamPert=0.0) override;
	double getU(std::size_t) override;
	double getC(std::size_t) override;
	double Gamma1(std::size_t) override;
	double sound_speed2(std::size_t, double GamPert=0.0) override;
	
		
	void writeStar(const char *const c=NULL) override;
	
private:
	void basic_setup();
	void init_arrays();
	std::size_t len;
	double n;		//polytropic index
	double Gamma;	//polytropic exponent
	//the values of K, rho0, P0 relate to the choice of units and scale factors
	//going to set all to 1
	//to rescale, basically "finish" the part of Rn we left off, sqrt[Pc/(Grho^2)]
	// radius scales like finished Rn
	// mass scales like rho0*Rn3
	double rho0;	//central density
	double P0;		//central pressure
	double Rn;		//radius scale factor
	//lane-emden solution functions
	enum VarName {x=0, y, z, numvar};
	double **Y;   //at each point X, variables x, y, z
	double *mass;
	double *base; //= pow(y,n-1), avoids repeated calls to pow(y,n)
	double GG;
	double set_mass(double const y[numvar]);
		
	//integrate Lane-Emden using basic RK4
	void centerInit(double ycenter[numvar]);
	void RK4step(double dx, double yin[numvar], double yout[numvar]);
	enum class SurfaceBehavior : bool {CONTINUE_FULL_LENGTH=false, STOP_AT_ZERO=true};
	double RK4integrate(const std::size_t, double, SurfaceBehavior);
	
	//methods for handling the BCs
	double ac[4], as[6];	//expansion coefficients of y near center, surface
	void setupCenter();		//prepare values of ac[]
	void setupSurface();	//prepare values of as[]
public:
	//methods to find central, surfae power series expansions of key variables in pulsation
	void getAstarCenter(double *, int&, double g=0) override;
	void getUCenter(double*, int&) override;
	void getVgCenter(double*, int&, double g=0) override;
	void getC1Center(double*, int&) override;
	void getAstarSurface(double*, int&, double g=0) override;
	void getUSurface(double*, int&) override;
	void getVgSurface(double*, int&, double g=0) override;
	void getC1Surface(double*, int&) override;
};

#endif