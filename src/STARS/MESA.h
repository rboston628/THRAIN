//**************************************************************************************
//							MESA Wrapper
// MESA.h
//		This is a simple wrapper, intended to work with data generated using MESA
//	
// Reece Boston, Sep 02, 2023
//**************************************************************************************

#ifndef MESACLASSH
#define MESACLASSH

#include "Star.h"

class MESA : public Star {
public:	
	MESA(const char*, std::size_t);	//constructor
	virtual ~MESA(); //destructor - clears all space in memory
	std::size_t length() override {return len;}
	double Mass() override;
	double Radius() override;
	double Gee() override;
	std::string graph_title () override {
		return strmakef("MESA model w/ Mass=%lg", Mtot);
	}
	
	double rad(std::size_t) override;
	double rho(std::size_t) override, drhodr(std::size_t) override;
	double   P(std::size_t) override,   dPdr(std::size_t) override;
	double Phi(std::size_t) override, dPhidr(std::size_t) override;
	double mr(std::size_t)  override;
	
	double Gamma1(std::size_t) override;
	double Schwarzschild_A(std::size_t, double GamPert=0.0) override;
	double getAstar(std::size_t, double GamPert=0.0) override;
	double getVg(std::size_t, double GamPert=0.0) override;
	double getU(std::size_t) override;
	double getC(std::size_t) override;
	double sound_speed2(std::size_t, double GamPert=0) override;

private:
	std::size_t Ntot, len, subgrid;  //Ntot = grid size from MESA
	double G, c2;
	//in units of g, cm, erg/s, respectively
	double Mtot, Rtot, Ltot;
	double Mstar, Rstar, Lstar;
	double Dscale, Pscale, Gscale;
	//arrays for density, pressure, mass, and gravitational field
	double *radi;	
	Splinor *dens, *pres, *mass, *grav, *Gam1, *BVfq;
	Splinor *aSpline, *vSpline, *uSpline, *cSpline;
	
	//for BCs of modes
	double nc;// effective polytropic index at center
	double ac[8]; //coefficients of theta near center
	double dc[3],pc[3];
	double A0[3],V0[3],c0[3],U0[3];
	//surface
	double ds[5], ps[5], ts[5], dels;
	double A1[5], V1[5], c1[5], U1[5];
	
	void setupCenter();
	void setupSurface();
	
public:
	void getAstarCenter(double *, int&, double g=0) override;
	void getUCenter(double*, int&) override;
	void getVgCenter(double*, int&, double g=0) override;
	void getC1Center(double*, int&) override;
	void getAstarSurface(double *, int&, double g=0) override;
	void getUSurface(double*, int&) override;
	void getVgSurface(double*, int&, double g=0) override;
	void getC1Surface(double*, int&) override;
	void writeStar(const char *const c=NULL) override;
	double SSR() override;
private:
	void printBV(const char *const, double const g=0) override;
	void printCoefficients(const char *const, double const g=0) override;
};

#endif