//**************************************************************************************
//							SIMPLE WHITE DWARF STAR
// SimpleWD.h                                                             Reece Boston
// Uses simple C/O core with He/H envelope and atmosphere                 Mar 24, 2022
//  User must specify total H, He, C, and O fractions within the star
//  Layer stratification is determined using the method described in 
//		- Reece Boston, PhD Thesis UNC, 2022
//  	- Boston, Clemens, Evans (2022)
//  The interior equation of state includes pressure contributions from
//		P = P_deg + P_ions + P_coulomb + P_rad
//		the electron degeneracy is approximated as T=0
//	Integrates with Schwarzschild's dimensionless logarithmic variables
// **************************************************************************************/

#ifndef SIMPLEWDH
#define SIMPLEWDH

#include "Star.h"
#include "../../lib/stellar.h"

extern PartialPressure ideal, rad_gas, coul, deg_zero, deg_finite, deg_trap, deg_partial;

class SimpleWD : public Star {
public:
	std::string graph_title() override {
		return strmakef("Simple WD, M=%1.3lf, R=%1.3lf, T_{eff}=%1.3le", Msolar, Rsolar*exp(logY[Ntot-1][radi]), Teff);
	}
	//constructor
	SimpleWD(
		double M, 	//the mass, in solar units
		double T,	//effective temperature in K
		std::size_t N
	);
	virtual ~SimpleWD();	//destructor
	std::size_t length() override {return Ntot;}
	//these three functions specify units
	double Radius() override;	//total radius
	double Mass() override;	    //total mass
	double Gee() override;
	//in Newtonian, light speed is infinity...
	virtual double light_speed2() override;
	
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
	double Ledoux(std::size_t, double GamPert=0.0);
	double BruntVaisala(std::size_t, double GamPert=0.0);
	Abundance Xmass;
		
private:
	void setup();
	void initFromChandrasekhar();
	StellarVar Ystart0, YstartS;
	double Y0; //the value y0 from the Chandrasekhar model
	void rescaleR();


	std::size_t Ntot, Ncore, Natm;//number of grid points
	
	//surfce values
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
	void setupGrid(double, std::size_t);
	void expandGrid(std::size_t);

	StellarVar Yscale, logYscale;
	StellarVar Ysolar;
	StellarVar Ystar;
	
	//methods to calculate the two regions
	static const int numv=3;
	void        joinAtCenter(double x[numv], double f[numv], double& F);
	double      calculateCore( const double x[numv], std::size_t Nmax);
	std::size_t firstCoreStep( const double x[numv], double& rholast, std::size_t Nmax);
	void        calculateAtmosphere( const double x[numv]);
	std::size_t firstAtmosphereStep( const double x[numv], double& rholast);
	StellarVar dYdR(       const StellarVar&, const double& delta, const double epsilon=1.0);
	StellarVar dlogYdlogR( const StellarVar&, const double& delta, const double epsilon=1.0);
	StellarVar dYdM(       const StellarVar&, const double& delta, const double epsilon=1.0);
	StellarVar dlogYdlogM( const StellarVar&, const double& delta, const double epsilon=1.0);
	
	//the physics
	double opacity(const StellarVar&, const Abundance&);
	double energyProduction(const StellarVar&, const Abundance&);
	double equationOfState( const StellarVar&, const Abundance&, double& rholast, int X=0);
	double energyTransport( const StellarVar&, const Abundance&);
	
	//abundances
	Abundance  Xtot;
	Abundance *Xelem, *dXelem;
	Abundance massFraction();
	//
	EOS core_pressure, atm_pressure;
	EOS* getEOS(const StellarVar& Y, const Abundance& X);
	
	double zy, zc, zo;
	double by, bc, bo;
	double my, mc, mo;
	Abundance findAbundance(const double, const double, Abundance&);
	
	//thermodynamic variables
	double *adiabatic_1, *nabla, *nabla_ad, *brunt_vaisala, *ledoux, *kappa;
	void populateBruntVaisala();
		
	//methods for handling the BCs
	double nc;
	double ac[8];
	double tc[4], dc[4], pc[4];
	double A0[4], V0[4], c0[4], U0[4];
	double ts[5], ds[5], ps[5];
	double A1[5], V1[5], c1[5], U1[5];
	void setupCenter();
	void setupSurface();
public:
	void getAstarCenter(double *, int&, double g=0) override;
	void getUCenter( double*, int&) override;
	void getVgCenter(double*, int&, double g=0) override;
	void getC1Center(double*, int&) override;
	void getAstarSurface(double *, int&, double g=0) override;
	void getUSurface( double*, int&) override;
	void getVgSurface(double*, int&, double g=0) override;
	void getC1Surface(double*, int&) override;
	
	//a particular output generation for this model of a WD
	void writeStar(const char *const c=NULL) override;
	double SSR() override;
private:
	void printChem(const char *const c);
	void printBV(const char *const c, const double g=0)  override;
	void printOpacity(const char *const c);
	void printBigASCII(const char *const c);
	void printCoefficients(const char *const c, const double g=0) override;
};

#endif