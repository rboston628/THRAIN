//**************************************************************************************
//							ZERO AGE MAIN SEQUENCE
// ZAMS.h
//	A Newtonian star at the beginning of the main sequence lifetime
//	Chemically homogeneous (X,Y,Z)
//	Integrates with Schwarzschild's dimensionless logarithmic variables
//  This code draws on inspiration from:
//		ZAMS.f by C.J. Hansen & S.D. Kawaler -- http://astro.if.ufrgs.br/evol/evolve/hansen/
//  	StatStar by Bradley W. Carroll & Dale A. Ostlie,  2007.
//  Reece Boston, Mar 24, 2022
// **************************************************************************************/

#ifndef ZAMSH
#define ZAMSH

#include "Star.h"
#include "../../lib/stellar.h"
#include "../../lib/chandra.cpp"
#include "../../lib/stellar.cpp"

extern PartialPressure ideal, rad_gas, coul, deg_zero, deg_finite, deg_trap, deg_partial;

class ZAMS : public Star {
public:
	void graph_title(char* buff){
		sprintf(buff, "Simple WD, M=%1.2f, R=%1.2f, T_{eff}=%1.2le", Msolar, Rsolar*exp(logY[Ntot-1][radi]), Teff);
	}
	//constructor
	ZAMS(
		double M, 	//the mass, in solar units
		double R,   //the radius, in solar units
		double T,	//effective temperature in K
		int N,
		double X, double Y //the metallicity Z is implied
	);
	virtual ~ZAMS();	//destructor
	int length(){return Ntot;}
	//these three functions specify units
	double Radius();	//total radius
	double Mass();	//total mass
	double Gee();
	//in Newtonian, light speed is infinity... just use the max value to represent this
	virtual double light_speed2();
	
	double rad(int);
	double rho(int), drhodr(int);
	double   P(int),   dPdr(int);
	double Phi(int), dPhidr(int);
	double mr(int);
	
	double Schwarzschild_A(int, double GamPert=0.0);
	double getAstar(int, double GamPert=0.0);
	double getVg(int, double GamPert=0.0);
	double getU(int);
	double getC(int);
	double Gamma1(int);
	double sound_speed2(int, double GamPert=0.0);
	double Ledoux(int, double GamPert=0.0);
	double BruntVaisala(int, double GamPert=0.0);
		
private:
	void setup();
	void initFromPolytrope();

	int Ntot, Ncore, Natm;//number of grid points
	
	//surface values
	double Msolar, Mstar; //in solar and CGS units
	double Lsolar, Lstar; //in solar and CGS units
	double Rsolar, Rstar; //in solar and CGS units
	double Teff;		  //effective temperature (K)
	//scale values
	double Dscale, Pscale, Tscale;

	//solution functions
	double *logQ;	   // the independent variable
	StellarVar*  logY; // log density, radius, pressure, mass, temperature, luminosity
	StellarVar* dlogY; // dlogY/dlogQ, as above
	void setupGrid(double, int);

	StellarVar Yscale, logYscale;
	StellarVar Ysolar;
	StellarVar Ystar;
	
	double opacity(const StellarVar&, const Abundance&);
	//
	StellarVar Ystart0, YstartS;
	
	//abundance
	Abundance  Xtot, Xmass, dXelem;
	//
	EOS core_pressure, atm_pressure;
	EOS* getEOS(const StellarVar& Y, const Abundance& X);
	
	
	//methods to calculate each of the three sections
	static const int numv=3;
	void   joinAtCenter(double x[numv], double f[numv], double& F);
	double calculateCore( const double x[numv], int Nmax);
	int    firstCoreStep( const double x[numv], double& rholast, int Nmax);
	void   calculateAtmosphere( const double x[numv]);
	int    firstAtmosphereStep( const double x[numv], double& rholast);
	StellarVar dYdR(       const StellarVar&, const double& delta, const double epsilon=1.0);
	StellarVar dlogYdlogR( const StellarVar&, const double& delta, const double epsilon=1.0);
	StellarVar dYdM(       const StellarVar&, const double& delta, const double epsilon=1.0);
	StellarVar dlogYdlogM( const StellarVar&, const double& delta, const double epsilon=1.0);
	
	//the physics
	double energyProduction(const StellarVar&, const Abundance&);
	double equationOfState( const StellarVar&, const Abundance&, double& rholast);
	double energyTransport( const StellarVar&, const Abundance&);

	Abundance findAbundance(const double, const double, Abundance&);
	Abundance massFraction();
	
	//thermodynamic variables
	double *adiabatic_1, *nabla, *nabla_ad, *brunt_vaisala, *ledoux, *kappa;
	void populateBruntVaisala();
		
	//methods for handling the BCs
	double nc;
	double ac[8];
	void setupCenter();
	void setupSurface();
public:
	void getAstarCenter(double *, int&, double g=0);
	void getUCenter( double*, int&);
	void getVgCenter(double*, int&, double g=0);
	void getC1Center(double*, int&);
	void getAstarSurface(double *, int&, double g=0);
	void getUSurface( double*, int&);
	void getVgSurface(double*, int&, double g=0);
	void getC1Surface(double*, int&);
	
	//a particular output generation for this model of a WD
	void writeStar(char *c=NULL);
	double SSR();
};

#endif