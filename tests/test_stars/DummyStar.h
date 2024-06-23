//**************************************************************************************
//							DUMMY STAR
// DummyStar.h
//		A simple object for use in testing.
//		It is only here to be used when something requries a Star object
//      for initializatio, but not any of itrs properties.
//		THIS IS NOT PART OF THRAIN -- ONLY FOR CERTAIN TESTS
//**************************************************************************************

#include "../../src/STARS/Star.h"

#ifndef DUMMYSTARH
#define DUMMYSTARH

class DummyStar : public Star {
public:
	std::string graph_title() override {
		return std::string("dummy");
	}
	
	DummyStar(std::size_t);	//constructor
	virtual ~DummyStar();	//destructor
	std::size_t length() override {return len;}
	//these three functions specify units
	double Radius() override;		//total radius
	double Mass() override;			//total mass
	double Gee() override;
		
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