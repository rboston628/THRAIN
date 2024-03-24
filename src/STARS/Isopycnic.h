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

#ifndef ISOPYCNICH
#define ISOPYCNICH

#include "Star.h"

class Isopycnic : public Star {
public:
	std::string graph_title() override {
		return std::string("uniform density");
	}
	
	Isopycnic(std::size_t);			//constructor, n and length
	virtual ~Isopycnic();	//destructor
	int length() override {return len;}
	//these three functions specify units
	double Radius() override;		//total radius
	double Mass() override;			//total mass
	double Gee() override {return GG;};

	double rad(int) override;
	double rho(int) override, drhodr(int) override;
	double   P(int) override,   dPdr(int) override;
	double Phi(int) override, dPhidr(int) override;
	double  mr(int) override;
	
	double Schwarzschild_A(int, double GamPert=0.0) override;
	double getAstar(int, double GamPert=0.0) override;
	double getVg(int, double GamPert=0.0) override;
	double getU(int) override;
	double getC(int) override;
	double Gamma1(int) override;
	double sound_speed2(int, double GamPert=0.0) override;
	
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
	int populateValues(const int, double);
	
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