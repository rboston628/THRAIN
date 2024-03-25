//**************************************************************************************
//							ISOPYCNIC STAR
// Isopycnic.h
//		A simple stellar object in Newtonian physics of uniform density
//				rho = constant
//		Thi exploits the known analytic solution to Lane-Emden for n=0
//				y = 1 - x^2/6
//		Used as stellsr background for mode testing, to minimize background errors
//		THIS IS NOT PART OF THRAIN -- ONLY FOR CERTAIN TESTS WITH NONRADIAL MODES
//**************************************************************************************

#include "Star.h"

#ifndef ISOPYCNICH
#define ISOPYCNICH

class Isopycnic : public Star {
public:
	std::string graph_title() override {
		return std::string("uniform density");
	}
	
	Isopycnic(std::size_t);			//constructor, n and length
	virtual ~Isopycnic();	//destructor
	std::size_t length() override {return len;}
	//these three functions specify units
	double Radius() override;		//total radius
	double Mass() override;			//total mass
	double Gee() override {return GG;};
		
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
	
private:
	std::size_t len;
	double Gamma;	//polytropic exponent
	double rho0;	//central density
	double P0;		//central pressure
	double Rn;		//radius scale factor
	//lane-emden solution functions
	double *x;	//normalized radius (xi)
	double *y;	//lane-emden solution (theta)
	double *z;	//derivative (dtheta/dxi)
	double *mass;
	double GG;
	
	//integrate Lane-Emden using basic RK4
	std::size_t populateValues(const std::size_t, double);
	
	//methods for handling the BCs
	void setupCenter();		//for conformity
	void setupSurface();	//for conformity
	char isoname[40];		//a name of the star
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