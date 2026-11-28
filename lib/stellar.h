#ifndef STELLARHELPSH
#define STELLARHELPSH

#include "../src/constants.h"
#include "chandra.h"

#include <array>

struct Abundance;

namespace chemical {
	enum elem {h=0, he, c, o};
	std::array<elem, 4> constexpr allElements = {h, he, c, o};
	const double Z[] = {1., 2.,  6.,  8.};
	const double A[] = {1., 4., 12., 16.};
	
	double partial_mean_A(elem i, Abundance const& X);
	double partial_mean_Z(elem i, Abundance const& X);
	double partial_mu_e(  elem i, Abundance const& X);
	double partial_mean_coul(elem i, Abundance const& X);
}

struct Abundance {
	double H1, He4, C12, O16;
	int e;
	double mean_A() const {
		using namespace chemical;
		return 1./(H1 + He4/chemical::A[he] + C12/chemical::A[c] + O16/chemical::A[o]);
	}
	double mean_Z() const {
		using namespace chemical;
		return H1 + He4*chemical::Z[he] + C12*chemical::Z[c] + O16*chemical::Z[o];
	}
	double mu_e(double y=1.0) const {
		using namespace chemical;
		return 1./(H1*y + 0.5*(He4+C12+O16)*y);
	}
	double mean_weight(double y=1.0) const {
		using namespace chemical;
		return 1./(1./mu_e(y) + 1./mean_A());
	}
	double mean_coulomb() const{
		using namespace chemical;
		return H1 
			+ He4*pow(chemical::Z[he],2)*pow(chemical::A[he],-1./3.) 
			+ C12*pow(chemical::Z[ c],2)*pow(chemical::A[ c],-1./3.) 
			+ O16*pow(chemical::Z[ o],2)*pow(chemical::A[ o],-1./3.) ;
	}
	void print() const {
		ThrainLogger::info("%0.2f\t%0.2f\t%0.2f\t%0.2f\n", H1, He4, C12, O16);
	}
	Abundance() :H1(0), He4(0), C12(0), O16(0), e(0) {};
	Abundance(double x1, double x2, double x3, double x4) 
		: H1(x1),He4(x2),C12(x3),O16(x4), e(0){};
	
	Abundance operator*(const double& x) {
		return Abundance(H1*x, He4*x, C12*x, O16*x);
	}
	Abundance operator+(const Abundance& x){
		return Abundance(H1+x.H1, He4+x.He4, C12+x.C12, O16+x.O16);
	}
	Abundance operator-(const Abundance& x){
		return Abundance(H1-x.H1, He4-x.He4, C12-x.C12, O16-x.O16);
	}
	Abundance& operator=(const Abundance& x) = default;
	bool operator==(const Abundance& x) const {
		return (H1==x.H1 && He4==x.He4 && C12==x.C12 && O16==x.O16);
	};
	double  operator[](int n) const {
		switch(n){
			case 0: return H1;  break;
			case 1: return He4; break;
			case 2: return C12; break;
			case 3: return O16; break;
			default: return 0.0;
		}
	}
	double& operator[](int n){
		switch(n){
			case 0: return H1;  break;
			case 1: return He4; break;
			case 2: return C12; break;
			case 3: return O16; break;
			default: return nothing; 
		}
	}
	void enforce(){
		if(H1 <0.0) H1  = 0.0;
		if(He4<0.0) He4 = 0.0;
		if(C12<0.0) C12 = 0.0;
		if(O16<0.0) O16 = 0.0;
		double const tot = H1+He4+C12+O16;
		H1 /= tot;
		He4/= tot;
		C12/= tot;
		O16/= tot;
	}
	private:
		static double nothing;
};

enum VarName {dens=0, radi, pres, mass, temp, lumi, num_var};

class StellarVar {
private:
	double var[num_var];
public:
	//access the elements
	double& operator[](const int n){ return (var[n]);};
	double  operator[](const int n) const { return (var[n]);};
	double& operator[](const VarName n){ return (var[n]);};
	double  operator[](const VarName n) const { return var[n];};
	//constructors
	StellarVar(const double x[num_var]){
		for(int j=0; j<num_var; j++) var[j] = x[j];
	};
	StellarVar(double x1, double x2, double x3, double x4, double x5, double x6){
		var[dens]=x1; var[radi]=x2; var[pres]=x3; var[mass]=x4; var[temp]=x5; var[lumi]=x6;
	};
	StellarVar(double x2, double x3, double x4, double x5, double x6){
		var[dens]=0.; var[radi]=x2; var[pres]=x3; var[mass]=x4; var[temp]=x5; var[lumi]=x6;
	}
	StellarVar() {
		for(int j=0; j<num_var; j++) var[j] = 0.0;\
	};
	
	//mathematical operations
	StellarVar operator*(const double &x) const {
		StellarVar val;
		for(int j=0; j<num_var; j++) val[j] = (*this)[j]*x;
		return val;
	}
	StellarVar operator/(const double &x) const {
		StellarVar val;
		for(int j=0; j<num_var; j++) val[j] = (*this)[j]/x;
		return val;
	}
	StellarVar operator+(const StellarVar &x) const {
		StellarVar val;
		for(int j=0; j<num_var; j++) val[j] = (*this)[j] + x[j];
		return val;
	}
	StellarVar operator-(const StellarVar &x) const {
		StellarVar val;
		for(int j=0; j<num_var; j++) val[j] = (*this)[j] - x[j];
		return val;
	}
	StellarVar operator*(const StellarVar &x) const {
		StellarVar val;
		for(int j=0; j<num_var; j++) val[j] = (*this)[j]*x[j];
		return val;
	}
	StellarVar operator/(const StellarVar &x) const {
		StellarVar val;
		for(int j=0; j<num_var; j++) val[j] = (*this)[j]/x[j];
		return val;
	}
		
	//assignment
	StellarVar& operator=(double const x[num_var]){
		for(int j=0; j<num_var; j++) var[j] = x[j];
		return *this;
	}
	StellarVar& operator=(StellarVar const& x){
		for(int j=0; j<num_var; j++) var[j] = x[j];
		return *this;
	}
};

//exponent and logarithm
StellarVar exp(const StellarVar &x);
StellarVar log(const StellarVar &x);

// **************************************************************************
// PARTIAL PRESSURE
// A structure describing some interaction within a stellar plasma
// A bundle of functions for pressure, energy, and partial derivatives
// Each function must take rho, T, and chemical composition X as arguments
// **************************************************************************
struct PartialPressure {
	double P(double const rho, double const T, Abundance const& X) const {return pfunc(rho,T,X);}
	double partialRho(double const rho, double const T, Abundance const& X) const {return prho(rho,T,X);}
	double partialT(double const rho, double const T, Abundance const& X) const {return pT(rho,T,X);}
	double partialX(chemical::elem i, double const rho, double const T, Abundance const& X) const {return pX(i,rho,T,X);}
	double U(double const rho, double const T, Abundance const& X) const {return ufunc(rho,T,X);};
	double UpartialRho(double const rho, double const T, Abundance const& X) const {return urho(rho,T,X);}
	double UpartialT(double const rho, double const T, Abundance const& X) const {return uT(rho,T,X);}
	double UpartialX(chemical::elem i, double const rho, double const T, Abundance const& X) const {return uX(i,rho,T,X);}
	double operator()(double const rho, double const T, Abundance const& chem) const {return P(rho,T,chem);}
	using funcptr = double (*)(double const, double const, Abundance const&);
	using funcptrX = double (*)(chemical::elem, double const, double const, Abundance const&);
	PartialPressure(funcptr p, funcptr pr, funcptr pt, funcptrX px,
		funcptr u, funcptr ur, funcptr ut, funcptrX ux)
		: pfunc(p), prho(pr), pT(pt), pX(px), ufunc(u), urho(ur), uT(ut), uX(ux) {};
private:
	funcptr pfunc, prho, pT;
	funcptrX pX;
	funcptr ufunc, urho, uT;
	funcptrX uX;
};

struct OneValueCache {
	double x, y;
	Abundance chem;
	double val;
	bool operator==(OneValueCache const& other) const {
		return (x==other.x && y==other.y && chem==other.chem);
	}
	bool cache_match(double const x, double const y, Abundance const& chem) const {
		return (this->x==x && this->y==y && this->chem==chem);
	}
	void recache(double const val, double const x, double const y, Abundance const& chem) {
		this->x = x; this->y = y; this->chem = chem; this->val = val;
	}
	OneValueCache() : x(0.0), y(0.0), chem(), val(nan("")) {};
};

// **************************************************************************
// ELECTRON PRESSURE
// The specific partial pressure of the plasma's electron gas
// This is described by an internal state, usually the degeneracy parameter
// The internal state must be described by rho, T, and chemical composition X
// Each other function must take state, T, and chemical composition X as arguments
// **************************************************************************
struct ElectronPressure : public PartialPressure {
	using funcptr = double (*)(double const, double const, Abundance const&);
	using funcptrX = double (*)(chemical::elem, double const, double const, Abundance const&);
	ElectronPressure(funcptr p, funcptr pr, funcptr pt, funcptrX px,
		funcptr u, funcptr ur, funcptr ut, funcptrX ux, funcptr es, funcptr rf) :
		 PartialPressure(p, pr, pt, px, u, ur, ut, ux), internal_state(es), rhofunc(rf), id(++next_id) {};
	ElectronPressure(ElectronPressure const&) = default;
	ElectronPressure& operator=(ElectronPressure const&) = default;
	bool operator==(ElectronPressure const& other) const { return id == other.id; }
	bool operator!=(ElectronPressure const& other) const { return id != other.id; }
	double state(double const rho, double const T, Abundance const& chem) { return internal_state(rho, T, chem); }
	double density(double const state, double const T, Abundance const &chem) { return rhofunc(state, T, chem); }
private:
	// the internal state, eg degeneracy, as function of rho, T, and chemical composition
	funcptr internal_state;
	// the total density, as function of internal state, T, and chemical composition
	funcptr rhofunc;
	// use this id for easier comparison of electron pressure objects
	static unsigned char next_id;
	unsigned char id;
};

extern ElectronPressure no_electrons;

class EOS {
public:
	EOS() = default;
	EOS(std::vector<PartialPressure> const& p) : pressure(p), density_cache(), electron_state_cache() {};
	EOS(EOS const &) = default;
	EOS& operator=(EOS const &) = default;
	double invert(double const rho_last, double const P, double const T, Abundance const& chem);
	double operator()(double const rho, double const T, Abundance const& chem);
	~EOS() = default;
	// other stuff
	double U(double const rho, double const T, Abundance const& chem);
	double Gamma1(double const rho, double const T, Abundance const& chem);
	double Gamma3(double const rho, double const T, Abundance const& chem);
	double nabla_ad(double const rho, double const T, Abundance const& chem);
	double chiRho(double const rho, double const T, Abundance const& chem);
	double chiT(double const rho, double const T, Abundance const& chem);
	double chiY(chemical::elem i, double const rho, double const T, Abundance const& chem);
	// calculate the Ledoux term -- see Brassard etc 1991, eqn (13, 14)
	double Ledoux(double const rho, double const T, Abundance const& chem, Abundance const& dchem);
	
	double electronState(double const rho, double const T, Abundance const& chem);

	// partial derivatives
	double partialRho(double const rho, double const T, Abundance const& chem);
	double partialT(double const rho, double const T, Abundance const& chem);
	
	PartialPressure operator[](int n);
	void push_back(PartialPressure const&);
	void push_back(ElectronPressure const&);
private:
	// the partial pressures are stored as function pointers in a vector object
	std::vector<PartialPressure> pressure {};
	ElectronPressure electron_pressure = no_electrons;
	// density cache to speed up repeated calls to find density
	OneValueCache density_cache;
	// electron state cache to speed up repeated calls to find degeneracy parameter
	OneValueCache electron_state_cache;
};

//************************************************************************************
//	PARTIAL PRESSURE FUNCTIONS for EQUATION OF STATE
//		Each must take arguments double const rho, double const T, Abundance const& X
//		Output must be a double, the pressure 
//		For references, see:
//		*	Chandrasekhar 1939
//		*	Cox & Giulu 1980
//		*	Pichon 1989
//		*	Tassoul, Fontaine, Winget, 1990
//		*	Van Horn 1968
//		*	Koester 1976
//************************************************************************************

//IDEAL GAS
double pressure_ideal(double const rho, double const T, Abundance const& X);
double partialRho_ideal(double const rho, double const T, Abundance const& X);
double partialT_ideal(double const rho, double const T, Abundance const& X);
double partialX_ideal(chemical::elem, double const rho, double const T, Abundance const& X);
double energy_ideal(double const rho, double const T, Abundance const& X);
double UpartialRho_ideal(double const rho, double const T, Abundance const& X);
double UpartialT_ideal(double const rho, double const T, Abundance const& X);
double UpartialX_ideal(chemical::elem, double const rho, double const T, Abundance const& X);


//RADIATION GAS
double pressure_rad(double const rho, double const T, Abundance const& X);
double partialRho_rad(double const rho, double const T, Abundance const& X);
double partialT_rad(double const rho, double const T, Abundance const& X);
double partialX_rad(chemical::elem, double const rho, double const T, Abundance const& X);
double energy_rad(double const rho, double const T, Abundance const& X);
double UpartialRho_rad(double const rho, double const T, Abundance const& X);
double UpartialT_rad(double const rho, double const T, Abundance const& X);
double UpartialX_rad(chemical::elem, double const rho, double const T, Abundance const& X);


//COULOMB PRESSURE
double pressure_coul(double const rho, double const T, Abundance const& X);
double partialRho_coul(double const rho, double const T, Abundance const& X);
double partialT_coul(double const rho, double const T, Abundance const& X);
double partialX_coul(chemical::elem, double const rho, double const T, Abundance const& X);
double energy_coul(double const rho, double const T, Abundance const& X);
double UpartialRho_coul(double const rho, double const T, Abundance const& X);
double UpartialT_coul(double const rho, double const T, Abundance const& X);
double UpartialX_coul(chemical::elem, double const rho, double const T, Abundance const& X);


//************************************************************************************
//	ELECTRON DEGENERACY PRESSURE TERMS
//  All of the following describe electron degeneracy pressure
//  They differ by evaluation technique used, or domain of validity
//************************************************************************************

//ELECTRON DEGENERACY PRESSURE -- ZERO TEMP
double density_deg_zero(double const X0, double const T, Abundance const& X);
double pressure_deg_zero(double const rho, double const T, Abundance const& X);
double partialRho_deg_zero(double const rho, double const T, Abundance const& X);
double partialT_deg_zero(double const rho, double const T, Abundance const& X);
double partialX_deg_zero(chemical::elem, double const rho, double const T, Abundance const& X);
double energy_deg_zero(double const rho, double const T, Abundance const& X);
double UpartialRho_deg_zero(double const rho, double const T, Abundance const& X);
double UpartialT_deg_zero(double const rho, double const T, Abundance const& X);
double UpartialX_deg_zero(chemical::elem, double const rho, double const T, Abundance const& X);

//ELECTRON DEGENERACY PRESSURE -- INTEGRATED GAUSS-LEGENDRE
// these are defined with template methods in stellar.cpp

//ELECTRON DEGENERACY PRESSURE -- QUADRATURE SAGAR 1991
// these are defined with templar methods in stellar.cpp

//ELECTRON DEGENERACY PRESSURE -- SWAP METHOD
//pressure of electron gas under different degrees of degeneracy
//  uses a combination of methods depending on degeneracy
//	following Cox and Giuli 1980 -- see equation 24.51, the definition of w
double pressure_deg_swap(double const rho, double const T, Abundance const& X);
double partialRho_deg_swap(double const rho, double const T, Abundance const& X);
double partialT_deg_swap(double const rho, double const T, Abundance const& X);
double partialX_deg_swap(chemical::elem i, double const rho, double const T, Abundance const& X);
double energy_deg_swap(double const rho, double const T, Abundance const& X);
double UpartialRho_deg_swap(double const rho, double const T, Abundance const& X);
double UpartialT_deg_swap(double const rho, double const T, Abundance const& X);
double UpartialX_deg_swap(chemical::elem i, double const rho, double const T, Abundance const& X);

#endif

