#ifndef STELLARHELPSCPP
#define STELLARHELPSCPP

#include "stellar.h"
#include "rootfind.h"
#include "logger.h"
#include <array>

double Abundance::nothing{0.0};

double chemical::partial_mean_A(elem i, Abundance const& X){
	return -pow(X.mean_A(),2)*
			(1./chemical::A[i] - 1./chemical::A[X.e]);
}

double chemical::partial_mean_Z(elem i, Abundance const& X){
	return chemical::Z[i] - chemical::Z[X.e];
}

double chemical::partial_mu_e(elem i, Abundance const& X){
	return -pow(X.mu_e(),2)*
			(chemical::Z[i]/chemical::A[i] - chemical::Z[X.e]/chemical::A[X.e]);
}

double chemical::partial_mean_coul(elem i , Abundance const& X){
	return pow(chemical::Z[i],2)*pow(chemical::A[i],-1./3.) 
			- pow(chemical::Z[X.e],2)*pow(chemical::A[X.e],-1./3.);
}

//exponent and logarithm
StellarVar exp(StellarVar const &x){
	StellarVar val;
	for(int j=0; j<num_var; j++) val[j] = std::exp(x[j]);
	return val;
}
StellarVar log(const StellarVar &x){
	StellarVar val;
	for(int j=0; j<num_var; j++) val[j] = std::log(x[j]);
	return val;
}

double EOS::electronState(double const rho, double const T, Abundance const& chem){
	if (!electron_state_cache.cache_match(rho, T, chem)) {
		electron_state_cache.recache(electron_pressure.state(rho,T,chem), rho, T, chem);
	}
	return electron_state_cache.val;
}

double EOS::operator()(double const rho, double const T, Abundance const& chem){
	double P = electron_pressure(electronState(rho,T,chem),T,chem);
	for(PartialPressure const &p : pressure) P += p(rho,T,chem);
	return P;
}

double EOS::U(double const rho, double const T, Abundance const& chem){
	double U = electron_pressure.U(electronState(rho,T,chem),T,chem);
	for(PartialPressure p : pressure) U += p.U(rho,T,chem);
	return U;
}

double EOS::Gamma1(double const rho, double const T, Abundance const& chem){
	double state = electronState(rho,T,chem);
	double P=electron_pressure(state,T,chem), U=electron_pressure.U(state,T,chem);
	double dPdT = electron_pressure.partialT(state,T,chem), dUdT = electron_pressure.UpartialT(state,T,chem);
	double dPdr = electron_pressure.partialRho(state,T,chem), dUdr = electron_pressure.UpartialRho(state,T,chem);
	for(PartialPressure const &p : pressure){
		P += p(rho,T,chem);
		U += p.U(rho,T,chem);
		dPdT += p.partialT(rho,T,chem);
		dUdT += p.UpartialT(rho,T,chem);
		dPdr += p.partialRho(rho,T,chem);
		dUdr += p.UpartialRho(rho,T,chem);
	}
	return rho/P*dPdr + U/P*dPdT/dUdT*(1. + P/U - rho/U*dUdr);
}
double EOS::Gamma3(double const rho, double const T, Abundance const& chem){
	double state = electronState(rho,T,chem);
	double dPdT = electron_pressure.partialT(state,T,chem), dUdT = electron_pressure.UpartialT(state,T,chem);
	for(PartialPressure const &p : pressure) dPdT += p.partialT(rho,T,chem);
	for(PartialPressure const &p : pressure) dUdT += p.UpartialT(rho,T,chem);
	return dPdT/dUdT + 1.;
}
double EOS::nabla_ad(double const rho, double const T, Abundance const& chem){
	return (Gamma3(rho,T,chem)-1.)/Gamma1(rho,T,chem);
}
double EOS::chiRho(double const rho, double const T, Abundance const& chem){
	double state = electronState(rho,T,chem);
	double P=electron_pressure(state,T,chem), dPdr = electron_pressure.partialRho(state,T,chem);
	for(PartialPressure const &p : pressure) {
		P += p(rho,T,chem);
		dPdr += p.partialRho(rho,T,chem);
	}
	return rho/P*dPdr;
}
double EOS::chiT(double const rho, double const T, Abundance const& chem){
	double state = electronState(rho,T,chem);
	double P=electron_pressure(state,T,chem), dPdT = electron_pressure.partialT(state,T,chem);
	for(PartialPressure const &p : pressure) {
		P += p(rho,T,chem);
		dPdT += p.partialT(rho,T,chem);
	}
	return T/P*dPdT;
}
double EOS::chiY(chemical::elem i, double const rho, double const T, Abundance const& chem){
	double state = electronState(rho,T,chem);
	double P = electron_pressure(state,T,chem), dPdY = electron_pressure.partialX(i, state,T,chem);
	for(PartialPressure const &p : pressure) {
		P += p(rho,T,chem);
		dPdY += p.partialX(i, rho,T,chem);
	}
	return chem[i]/P*dPdY;
}
//calculate the Ledoux term -- see Brassard etc 1991, eqn (13, 14)
// expanding out the logarithmic derivatives
// B = -1/chiT*Sum[ (delP/delX)_rho,T * dX/dP ]
double EOS::Ledoux(double const rho, double const T, Abundance const& chem, Abundance const& dchem){
	double state = electronState(rho,T,chem);
	double P = electron_pressure(state,T,chem), dPdT = electron_pressure.partialT(state,T,chem);
	double dPdY = 0.0;
	for (chemical::elem const &e : chemical::allElements) {
		dPdY += electron_pressure.partialX(e, state,T,chem)*dchem[e];
	}
	// the partialX terms are not stored in totals -- compute
	for(PartialPressure const &p : pressure){
		P += p(rho,T,chem);
		dPdT += p.partialT(rho,T,chem);
		for (chemical::elem const &e : chemical::allElements) {
			dPdY += p.partialX(e, rho,T,chem)*dchem[e];
		}
	}
	double chiT = T/P*dPdT;
	return -dPdY/chiT;
}

double EOS::partialRho(double const rho, double const T, Abundance const& chem){
	double dPdr = electron_pressure.partialRho(electronState(rho,T,chem),T,chem);
	for(PartialPressure const &p : pressure) dPdr += p.partialRho(rho,T,chem);
	return dPdr;
}
	
double EOS::partialT(double const rho, double const T, Abundance const& chem){
	double dPdT = electron_pressure.partialT(electronState(rho,T,chem),T,chem);
	for(PartialPressure const &p : pressure) dPdT += p.partialT(rho,T,chem);
	return dPdT;
}

//retrieve a particular partial-pressure term
PartialPressure EOS::operator[](int n){
	return pressure[n];
}

// Invert the usual EOS P = P(rho,T,X), to instead solve for rho given P,T,X
// uses a combination of Newton method and bisection search
double EOS::invert(double const rho_last, double const P, double const T, Abundance const& chem){	
	// check the cache first
	if (density_cache.cache_match(P, T, chem)) {
		return density_cache.val;
	}
	// if not, we have to rootfind the density from the EOS
	double startingGuess;
	std::function<double(double)> pressureDifference;
	// if a non-zero electron pressure is present, use that to find the density
	// this will use the electron degeneracy as the rootfinding variable
	// this avoids a double-inversion when the electron pressure searches for eta
	if (electron_pressure != no_electrons) {
		startingGuess = electron_pressure.state(rho_last, T, chem);
		pressureDifference = [this,P,T, &chem](double eta)->double {
			// find the density from the electron pressure's internal state (eg degeneracy parameter)
			double rho = electron_pressure.density(eta, T, chem);
			double Ptot = electron_pressure(eta, T, chem);
			for(PartialPressure const &p : pressure) Ptot += p(rho, T, chem);
			// return the difference between the target pressure and the total pressure
			return P - Ptot;
		};
	}
	// otherwise we use the density as an input to the EOS and rootfind for pressure
	else {
		startingGuess = rho_last;
		pressureDifference = [this,P,T, &chem](double rho)->double {
			return P - (*this)(rho, T, chem);
		};
	}

	double x=startingGuess, xmin, xmax; // NOTE cached value not guaranteed to be initialized
	if (rootfind::bisection_find_brackets_newton(pressureDifference, x, xmin, xmax)) {
		ThrainLogger::error("EOS::invert: failed to find brackets for P=%g, T=%g, rho_last=%g\n", P, T, rho_last);
		ThrainLogger::error("\tbrackets: [%g, %g]\n", xmin, xmax);
		return std::nan("");
	}
	double y = rootfind::bisection_search(pressureDifference, x, xmin, xmax);
	if (!std::isfinite(y)) {
		ThrainLogger::error("EOS::invert: failed to find root for P=%g, T=%g, rho_last=%g\n", P, T, rho_last);
		ThrainLogger::error("\tbrackets: [%g, %g]\n", xmin, xmax);
		return std::nan("");
	} // NOTE y ~= 0 if this worked

	// if we used internal state, convert back to density
	if (electron_pressure != no_electrons) {
		double density = electron_pressure.density(x, T, chem);
		electron_state_cache.recache(x, density, T, chem);
		x = density;
	}
	// store the resulting density (x) in the cache for next time
	density_cache.recache(x, P, T, chem);
	return x;
}

//  Adding new terms to the EOS will add them to the vector of function pointers
void EOS::push_back(PartialPressure const& p){
	pressure.push_back(p);
}

void EOS::push_back(ElectronPressure const& p){
	if (electron_pressure == no_electrons) {
		electron_pressure = p;
	} else {
		ThrainLogger::error("EOS::push_back: electron pressure already set; cannot add another\n");
	}
}


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
PartialPressure ideal = {
	pressure_ideal, partialRho_ideal, partialT_ideal, partialX_ideal,
	energy_ideal, UpartialRho_ideal, UpartialT_ideal, UpartialX_ideal};
//pressure of an ideal gas
double pressure_ideal(double const rho, double const T, Abundance const& X){
	return N_Avogadro*boltzmann_k/X.mean_A()*rho*T;
}
double partialRho_ideal(double const rho, double const T, Abundance const& X){
	return N_Avogadro*boltzmann_k/X.mean_A()*T;
}
double partialT_ideal(double const rho, double const T, Abundance const& X){
	return N_Avogadro*boltzmann_k/X.mean_A()*rho;
}
double partialX_ideal(chemical::elem i, double const rho, double const T, Abundance const& X){
	return -N_Avogadro*boltzmann_k/pow(X.mean_A(),2)*rho*T*chemical::partial_mean_A(i,X);
}
double energy_ideal(double const rho, double const T, Abundance const& X){
	return 1.5*pressure_ideal(rho,T,X);
}
double UpartialRho_ideal(double const rho, double const T, Abundance const& X){
	return 1.5*partialRho_ideal(rho,T,X);
}
double UpartialT_ideal(double const rho, double const T, Abundance const& X){
	return 1.5*partialT_ideal(rho,T,X);
}
double UpartialX_ideal(chemical::elem i, double const rho, double const T, Abundance const& X){
	return 1.5*partialX_ideal(i, rho,T,X);
}

//RADIATION GAS
//standard radiation pressure
PartialPressure rad_gas = {
	pressure_rad, partialRho_rad, partialT_rad, partialX_rad,
	energy_rad, UpartialRho_rad, UpartialT_rad, UpartialX_rad};
double pressure_rad(double const rho, double const T, Abundance const& X){
	return radiation_a/3.*pow(T,4);
}
double partialRho_rad(double const rho, double const T, Abundance const& X){
	return 0.0;
}
double partialT_rad(double const rho, double const T, Abundance const& X){
	return 4./3.*radiation_a*pow(T,3);
}
double partialX_rad(chemical::elem i, double const rho, double const T, Abundance const& X){
	return 0.0;
}
double energy_rad(double const rho, double const T, Abundance const& X){
	return radiation_a*pow(T,4);
}
double UpartialRho_rad(double const rho, double const T, Abundance const& X){
	return 0.0;
}
double UpartialT_rad(double const rho, double const T, Abundance const& X){
	return 4.*radiation_a*pow(T,3);
}
double UpartialX_rad(chemical::elem i, double const rho, double const T, Abundance const& X){
	return 0.0;
}


//COULOMB PRESSURE
//pressure contribution of ideal gas due to coloumb  correction
//  following EOS of Tassoul, Fontaine, Winget 1990
//  see also Van Horn 1968, Koester 1976
//  According to Huebner 2014, eq 9.31, GC should have a power of 4 and not a 5...
//  According to Cassisi et al 2007, eq 7, power should be 5
PartialPressure coul = {
	pressure_coul, partialRho_coul, partialT_coul, partialX_coul,
	energy_coul, UpartialRho_coul, UpartialT_coul, UpartialX_coul};
double pressure_coul(double const rho, double const T, Abundance const& X){
	double GC = 2.275e5*(X.mean_coulomb())*pow(rho,1./3.)/T;
	if(T<=0.0) GC = 0.0;
	double coul = -0.3*pow(GC,1.5)/(pow(GC,0.5)+1.03921);
	return coul*pressure_ideal(rho,T,X);
}
double partialRho_coul(double const rho, double const T, Abundance const& X){
	double GC = 2.275e5*(X.mean_coulomb())*pow(rho,1./3.)/T;
	if(T<=0.0) GC = 0.0;
	double denom = (pow(GC,0.5)+1.03921);
	double coul = -0.3*pow(GC,1.5)/denom;
	double dGCdr = 2.275e5*(X.mean_coulomb())*(1./3.)*pow(rho,-2./3.)/T;
	return coul*partialRho_ideal(rho,T,X)
			 + pressure_ideal(rho,T,X)*( -0.45*pow(GC,0.5)/denom + 0.15*GC/denom/denom )*dGCdr;
}
double partialT_coul(double const rho, double const T, Abundance const& X){
	double GC = 2.275e5*(X.mean_coulomb())*pow(rho,1./3.)/T;
	if(T<=0.0) GC = 0.0;
	double denom = (pow(GC,0.5)+1.03921);
	double coul = -0.3*pow(GC,1.5)/denom;
	double dGCdt = -2.275e5*(X.mean_coulomb())*pow(rho,1./3.)/T/T;
	return coul*partialT_ideal(rho,T,X)
			 + pressure_ideal(rho,T,X)*( -0.45*pow(GC,0.5)/denom + 0.15*GC/denom/denom  )*dGCdt;
}
double partialX_coul(chemical::elem i, double const rho, double const T, Abundance const& X){
	double GC = 2.275e5*(X.mean_coulomb())*pow(rho,1./3.)/T;
	if(T<=0.0) GC = 0.0;
	double denom = (pow(GC,1./2.)+1.03921);
	double coul = -0.3*pow(GC,3./2.)/denom;
	double dGCdX = 2.275e5*pow(rho,1./3.)/T*chemical::partial_mean_coul(i, X);
	return coul*partialX_ideal(i, rho,T,X)
			 + pressure_ideal(rho,T,X)*( -0.45*pow(GC,0.5)/denom + 0.15*GC/denom/denom  )*dGCdX;
}
double energy_coul(double const rho, double const T, Abundance const& X){
	return 3.*pressure_coul(rho,T,X);
}
double UpartialRho_coul(double const rho, double const T, Abundance const& X){
	return 3.*partialRho_coul(rho,T,X);
}
double UpartialT_coul(double const rho, double const T, Abundance const& X){
	return 3.*partialT_coul(rho,T,X);
}
double UpartialX_coul(chemical::elem i, double const rho, double const T, Abundance const& X){
	return 3.*partialX_coul(i, rho,T,X);
}

//************************************************************************************
//	ELECTRON PARTIAL PRESSURES and EQUATIONS OF STATE
//  These functions describe the pressure and energy of a degenerate electron gas
//	Each must take arguments double const rho, double const T, Abundance const& X
//************************************************************************************

namespace {
	double zero(double const, double const, Abundance const&) {return 0.0;}
	double zero_x(chemical::elem, double const, double const, Abundance const&) {return 0.0;}
	double rho_out(double rho, double T, Abundance const& X) { return rho;}
}

unsigned char ElectronPressure::next_id{0};
ElectronPressure no_electrons = {
	zero, zero, zero, zero_x,
	zero, zero, zero, zero_x,
	rho_out, rho_out
};

//ELECTRON DEGENERACY PRESSURE -- ZERO TEMP
//pressure of completely degenerate electron gas (T=0)
//	following Chandrasekhar 1939
namespace {
double calc_x0(double const rho, double const T, Abundance const& X){
	return pow(rho/Chandrasekhar::B0/X.mu_e(), 1./3.);
}
}
ElectronPressure deg_zero = {
	pressure_deg_zero, partialRho_deg_zero, partialT_deg_zero, partialX_deg_zero,
	energy_deg_zero, UpartialRho_deg_zero, UpartialT_deg_zero, UpartialX_deg_zero,
	rho_out, rho_out
};
double density_deg_zero(double const rho, double const T, Abundance const& X){
	using namespace Chandrasekhar;
	double X0 = calc_x0(rho, T, X);
	return B0*X.mu_e()*pow(X0,3);
}
double pressure_deg_zero(double const rho, double const T, Abundance const& X){
	using namespace Chandrasekhar;
	double X0 = calc_x0(rho, T, X);
	return A0*factor_f(X0);
}
double partialRho_deg_zero(double const rho, double const T, Abundance const& X){
	using namespace Chandrasekhar;
	double X0 = calc_x0(rho, T, X);
	return A0*8.*pow(X0,4)/sqrt(1.+X0*X0) * (1./(3.*B0*X.mu_e()*X0*X0));
}
double partialT_deg_zero(double const rho, double const T, Abundance const& X){
	return 0.0;
}
double partialX_deg_zero(chemical::elem i, double const rho, double const T, Abundance const& X){
	using namespace Chandrasekhar;
	double X0 = calc_x0(rho, T, X);
	return A0*8.*pow(X0,4)/sqrt(1.+X0*X0) 
		* pow(rho/B0,1./3.)*(-1./3.*pow(X.mu_e(),-4./3.))*chemical::partial_mu_e(i,X);
}
double energy_deg_zero(double const rho, double const T, Abundance const& X){
	using namespace Chandrasekhar;
	double X0 = calc_x0(rho, T, X);
	return A0*factor_g(X0);
}
double UpartialRho_deg_zero(double const rho, double const T, Abundance const& X){
	using namespace Chandrasekhar;
	double X0 = calc_x0(rho, T, X);
	return A0*24.*X0*X0*(sqrt(1.+X0*X0)-1.) * (1./(3.*B0*X.mu_e()*X0*X0));
}
double UpartialT_deg_zero(double const rho, double const T, Abundance const& X){
	return 0.0;
}
double UpartialX_deg_zero(chemical::elem i, double const rho, double const T, Abundance const& X){
	using namespace Chandrasekhar;
	double X0 = calc_x0(rho, T, X);
	return A0*24.*X0*X0*(sqrt(1.+X0*X0)-1.)
		* pow(rho/B0,1./3.)*(-1./3.*pow(X.mu_e(),-4./3.))*chemical::partial_mu_e(i,X);
}

namespace {
// DEGENERATE ELECTRONS WITH FERMI-DIRAC STATISTICS
// these static methods help abstract out Fermi-Dirac statistics for degenerate electrons
// this allows for different means of calculaing the F_k functions with unified handling
struct FermiDiracSet {
	using FD = double (*)(double, double);
	FD F[3], dFdEta[3], dFdBeta[3];
};

FermiDiracSet const fermi_sagar {
	{FermiDirac::FermiDiracSagar1Half, FermiDirac::FermiDiracSagar3Half, FermiDirac::FermiDiracSagar5Half},
	{FermiDirac::FermiDiracSagar1HalfDelEta, FermiDirac::FermiDiracSagar3HalfDelEta, FermiDirac::FermiDiracSagar5HalfDelEta},
	{FermiDirac::FermiDiracSagar1HalfDelBeta, FermiDirac::FermiDiracSagar3HalfDelBeta, FermiDirac::FermiDiracSagar5HalfDelBeta}
};

FermiDiracSet const fermi_trap {
	{FermiDirac::FermiDirac1Half, FermiDirac::FermiDirac3Half, FermiDirac::FermiDirac5Half},
	{FermiDirac::FermiDirac1HalfDelEta, FermiDirac::FermiDirac3HalfDelEta, FermiDirac::FermiDirac5HalfDelEta},
	{FermiDirac::FermiDirac1HalfDelBeta, FermiDirac::FermiDirac3HalfDelBeta, FermiDirac::FermiDirac5HalfDelBeta}
};

// calculate the coefficient in front of the Fermi-Dirac F_k functions for density
double fermi_dirac_rhocoeff(double const beta, double const mue) {
	using namespace Chandrasekhar;
	return 3.*sqrt(2.)*B0*mue*pow(beta, 1.5);
}
double fermi_dirac_Pcoeff(double const beta) {
	using namespace Chandrasekhar;
	return 16.*sqrt(2.)*A0*pow(beta, 2.5);
}
double fermi_dirac_Ucoeff(double const beta) {
	using namespace Chandrasekhar;
	return 24.*sqrt(2.)*A0*pow(beta, 2.5);
}

double relativistic_beta(double const T) {
	using namespace Chandrasekhar;
	return boltzmann_k*T/(electron.mass_CGS*C_CGS*C_CGS);
}

template <FermiDiracSet const *const fd>
double fd_rho(double const eta, double const T, Abundance const &X) {
	double beta = relativistic_beta(T);
	return fermi_dirac_rhocoeff(beta, X.mu_e()) * ( fd->F[0](eta, beta) + beta * fd->F[1](eta, beta) );
}

template <FermiDiracSet const *const fd>
double solveEta(double const rho, double const T, Abundance const &X) {
	using namespace Chandrasekhar;

	// unhappy path: rho or beta is invalid
	if (std::isnan(rho) || std::isnan(T) || rho < 0.0) {
		ThrainLogger::error("ERROR in solveEta: invalid rho = %le, T = %le\n", rho, T);
		return nan("");
	}

	double beta = relativistic_beta(T);
	double mue = X.mu_e();
	// approximate eta for degenerate cases: as in Cox and Giuli eq 24.29', pg 793 (but using 24.327b)
	double const eta_coeff = pow(8.*B0*B0*mue*mue, -1./3.);
	double const eta_deg = eta_coeff*pow(rho, 2./3.)/beta;
	// approximate eta for non-degenerate cases: F_{1/2}(eta) --> Gamma(3/2) * exp(eta)
	double const rho_coeff = fermi_dirac_rhocoeff(beta, mue);
	double const eta_nondeg = log(2.*rho/(rho_coeff*sqrt(math::PI)));
	// determine which to use
	double const guess = (eta_deg > 1.0 ? std::max(eta_deg, eta_nondeg) : eta_nondeg);

	// unhappy path: the guess is invalid
	if (!std::isfinite(guess)) {
		ThrainLogger::error("ERROR in solveEta: invalid guess = %le\n", guess);
		return nan("");
	}

	// now find the value
	std::function<double(double)> densDiff = [rho, T, &X](double eta)->double {
		return fd_rho<fd>(eta, T, X) - rho;
	};
	double eta_lo, eta_hi;
	if (rootfind::bisection_find_brackets_move(densDiff, guess, eta_lo, eta_hi)) {
		ThrainLogger::error("ERROR in solveEta: failed to find brackets for rho=%g, beta=%g, mue=%g\n", rho, beta, mue);
		return nan("");
	}
	double eta = guess;
	double y = rootfind::bisection_search(densDiff, eta, eta_lo, eta_hi);
	if (std::isnan(y)) {
		ThrainLogger::error("ERROR in solveEta: failed to find root for rho=%g, beta=%g, mue=%g\n", rho, beta, mue);
		return nan("");
	}
	return eta;
}

template <FermiDiracSet const *const fd>
double fd_pressure(double const eta, double const T, Abundance const &X) {	
	double const beta = relativistic_beta(T);
	return fermi_dirac_Pcoeff(beta) * ( fd->F[1](eta, beta) + 0.5 * beta * fd->F[2](eta, beta) );
}

template <FermiDiracSet const *const fd>
double fd_energy(double const eta, double const T, Abundance const &X) {
	double const beta = relativistic_beta(T);
	return fermi_dirac_Ucoeff(beta) * ( fd->F[1](eta, beta) + beta * fd->F[2](eta, beta) );
}

template <FermiDiracSet const *const fd>
double fd_p_del_rho(double const eta, double const T, Abundance const &X) {
	double const beta = relativistic_beta(T);
	double const mue = X.mu_e();
	double dfde12 = fd->dFdEta[0](eta, beta), dfde32 = fd->dFdEta[1](eta, beta), dfde52 = fd->dFdEta[2](eta, beta);
	return fermi_dirac_Pcoeff(beta) / fermi_dirac_rhocoeff(beta, mue) * 
		(dfde32 + 0.5 * beta * dfde52) / (dfde12 + beta * dfde32 );
}

template <FermiDiracSet const *const fd>
double fd_p_del_T(double const eta, double const T, Abundance const &X) {
	double const beta = relativistic_beta(T);
	double f12 = fd->F[0](eta, beta), f32 = fd->F[1](eta, beta), f52 = fd->F[2](eta, beta);
	double dfde12 = fd->dFdEta[0](eta, beta), dfde32 = fd->dFdEta[1](eta, beta), dfde52 = fd->dFdEta[2](eta, beta);
	double dfdb12 = fd->dFdBeta[0](eta, beta), dfdb32 = fd->dFdBeta[1](eta, beta), dfdb52 = fd->dFdBeta[2](eta, beta);
	return 16.*sqrt(2.)*Chandrasekhar::A0*(beta/T)*(
		2.5*pow(beta,1.5)*(	f32 + 0.5*beta*f52 )
		+ pow(beta,2.5)*( dfdb32 + 0.5*f52 + 0.5*beta*dfdb52 )
		- pow(beta,2.5)*( dfde32 + 0.5*beta*dfde52 ) * ( 
				1.5/beta * (f12 + beta*f32) + dfdb12 + f32 + beta * dfdb32
			) / ( dfde12 + beta*dfde32 )
	);
}

template <FermiDiracSet const *const fd>
double fd_p_del_X(chemical::elem i, double const eta, double const T, Abundance const &X) {
	double const beta = relativistic_beta(T);
	double const mue = X.mu_e();
	double dfde12 = fd->dFdEta[0](eta, beta), dfde32 = fd->dFdEta[1](eta, beta), dfde52 = fd->dFdEta[2](eta, beta);
	return -16.*sqrt(2.)*Chandrasekhar::A0*pow(beta,2.5)
		*(fd->F[0](eta, beta)+beta*fd->F[1](eta, beta))
		*(dfde32 + 0.5*beta*dfde52) / (dfde12 + beta * dfde32)
		*chemical::partial_mu_e(i,X)/mue;
}

template <FermiDiracSet const *const fd>
double fd_u_del_rho(double const eta, double const T, Abundance const &X) {
	double const beta = relativistic_beta(T);
	double const mue = X.mu_e();
	double dfde12 = fd->dFdEta[0](eta, beta), dfde32 = fd->dFdEta[1](eta, beta), dfde52 = fd->dFdEta[2](eta, beta);
	double dfdb12 = fd->dFdBeta[0](eta, beta), dfdb32 = fd->dFdBeta[1](eta, beta), dfdb52 = fd->dFdBeta[2](eta, beta);
	return fermi_dirac_Ucoeff(beta) / fermi_dirac_rhocoeff(beta, mue) * 
	  ( fd->dFdEta[1](eta, beta) + beta * fd->dFdEta[2](eta, beta) ) /
		( fd->dFdEta[0](eta, beta) + beta * fd->dFdEta[1](eta, beta) );
}

template <FermiDiracSet const *const fd>
double fd_u_del_T(double const eta, double const T, Abundance const &X) {
	double const beta = relativistic_beta(T);
	return 24.*sqrt(2.)*Chandrasekhar::A0*(beta/T)*(
		2.5*pow(beta,1.5)*(	fd->F[1](eta, beta) + beta*fd->F[2](eta, beta) )
		+ pow(beta,2.5)*(fd->dFdBeta[1](eta, beta) + (1.+beta)*fd->dFdBeta[2](eta, beta))
	);
}

template <FermiDiracSet const *const fd>
double fd_u_del_X(chemical::elem i, double const eta, double const T, Abundance const &X) {
	double const beta = relativistic_beta(T);
	double const mue = X.mu_e();
	return -24.*sqrt(2.)*Chandrasekhar::A0*pow(beta,2.5)
		*(fd->F[0](eta, beta)+beta*fd->F[1](eta, beta))
		*(fd->dFdEta[1](eta, beta)+beta*fd->dFdEta[2](eta, beta))
		/(fd->dFdEta[0](eta, beta)+    beta*fd->dFdEta[1](eta, beta))
		*chemical::partial_mu_e(i,X)/mue;
}

template <FermiDiracSet const *const fd>
ElectronPressure make_fd_partial_pressure() {
	return {
		fd_pressure<fd>, fd_p_del_rho<fd>, fd_p_del_T<fd>, fd_p_del_X<fd>,
		fd_energy<fd>, fd_u_del_rho<fd>, fd_u_del_T<fd>, fd_u_del_X<fd>,
		solveEta<fd>, fd_rho<fd>
	};
}
} // namespace


//ELECTRON DEGENERACY PRESSURE -- FINITE TEMP
// Solve for the Fermi-Dirac integral functions by numerical integration
// the exact method splits the domain into several regions, 
// each integrated with Gauss-Legendre 8-point quadrature
ElectronPressure deg_integrate = make_fd_partial_pressure<&fermi_trap>();


//ELECTRON DEGENERACY PRESSURE -- SAGAR QUADRATURE
// Solve the Fermi-Dirac integrals using a quadrature method due to Sagar,
// See papers:
// - Robin P Sagar, Computer Physics Communications vol 66, 1991
// - Robin P Sagar, ApJ 376, 1991
ElectronPressure deg_sagar = make_fd_partial_pressure<&fermi_sagar>();
	

//ELECTRON DEGENERACY PRESSURE -- SWAP METHOD
// Allow for swapping between different degeneracy pressures based on a parameter w
// For high degeneracy, will use the zero-T approximation
// For intermediate degeneracy, will use the Sagar quadrature
// For low degeneracy, will use the ideal gas approximation
// See Cox and Giuli 1980, quation 24.51, for the definition of w
ElectronPressure deg_swap = {
	pressure_deg_swap, partialRho_deg_swap, partialT_deg_swap, partialX_deg_swap,
	energy_deg_swap, UpartialRho_deg_swap, UpartialT_deg_swap, UpartialX_deg_swap,
	rho_out, rho_out
};
namespace {
std::array<double, 2> findW(double const rho, double const T, Abundance const& X){
	double const w_coeff = 0.5*pow(planck_h_CGS,3)/pow(2.*m_pi*electron.mass_CGS*boltzmann_k,1.5);
	double const ne = rho/proton.mass_CGS/X.mu_e();
	return {w_coeff*ne*pow(T, -1.5), ne};
}
}
double pressure_deg_swap(double const rho, double const T, Abundance const& X){
	std::array<double, 2> const wne = findW(rho,T,X);
	if      (wne[0] >10.0) return pressure_deg_zero(rho,T,X);
	else if (wne[0] >= 0.1) {
		double eta = solveEta<&fermi_sagar>(rho,T,X);
		return fd_pressure<&fermi_sagar>(rho,T,X);
	}
	else if (wne[0] <  0.1) return wne[1]*boltzmann_k*T;
	return 0.0;
}
double partialRho_deg_swap(double const rho, double const T, Abundance const& X){
	std::array<double, 2> const wne = findW(rho,T,X);
	if      (wne[0] >10.0) return partialRho_deg_zero(rho,T,X);
	else if (wne[0] >= 0.1) {
		double eta = solveEta<&fermi_sagar>(rho,T,X);
		return fd_p_del_rho<&fermi_sagar>(rho,T,X);
	}
	else if (wne[0] <  0.1) return boltzmann_k*T/proton.mass_CGS/X.mu_e();
	return 0.0;
}
double partialT_deg_swap(double const rho, double const T, Abundance const& X){
	std::array<double, 2> const wne = findW(rho,T,X);
	if      (wne[0] >10.0) return partialT_deg_zero(rho,T,X);
	else if (wne[0] >= 0.1) {
		double eta = solveEta<&fermi_sagar>(rho,T,X);
		return fd_p_del_T<&fermi_sagar>(rho,T,X);
	}
	else if (wne[0] <  0.1) return wne[1]*boltzmann_k;
	return 0.0;
}
double partialX_deg_swap(chemical::elem i, double const rho, double const T, Abundance const& X){
	std::array<double, 2> const wne = findW(rho,T,X);
	if      (wne[0] >10.0) 	return partialX_deg_zero(i,rho,T,X);
	else if (wne[0] >= 0.1) {
		double eta = solveEta<&fermi_sagar>(rho,T,X);
		return fd_p_del_X<&fermi_sagar>(i,eta,T,X);
	}
	else if (wne[0] <  0.1) return rho*boltzmann_k*T*pow(X.mu_e(),-2)*chemical::partial_mu_e(i,X);
	return 0.0;
}
double energy_deg_swap(double const rho, double const T, Abundance const& X){
	std::array<double, 2> const wne = findW(rho,T,X);
	if      (wne[0] >10.0) return energy_deg_zero(rho,T,X);
	else if (wne[0] >= 0.1) {
		double eta = solveEta<&fermi_sagar>(rho,T,X);
		return fd_energy<&fermi_sagar>(rho,T,X);
	}
	else if (wne[0] <  0.1) return 1.5*wne[1]*boltzmann_k*T;
	return 0.0;
}
double UpartialRho_deg_swap(double const rho, double const T, Abundance const& X){
	std::array<double, 2> const wne = findW(rho,T,X);
	if      (wne[0] >10.0) return UpartialRho_deg_zero(rho,T,X);
	else if (wne[0] >= 0.1) {
		double eta = solveEta<&fermi_sagar>(rho,T,X);
		return fd_u_del_rho<&fermi_sagar>(rho,T,X);
	}
	else if (wne[0] <  0.1) return 1.5*boltzmann_k*T/proton.mass_CGS/X.mu_e();
	return 0.0;
}
double UpartialT_deg_swap(double const rho, double const T, Abundance const& X){
	std::array<double, 2> const wne = findW(rho,T,X);
	if      (wne[0] >10.0) return UpartialT_deg_zero(rho,T,X);
	else if (wne[0] >= 0.1) {
		double eta = solveEta<&fermi_sagar>(rho,T,X);
		return fd_u_del_T<&fermi_sagar>(rho,T,X);
	}
	else if (wne[0] <  0.1) return 1.5*wne[1]*boltzmann_k;
	return 0.0;
}
double UpartialX_deg_swap(chemical::elem i, double const rho, double const T, Abundance const& X){
	std::array<double, 2> const wne = findW(rho,T,X);
	if      (wne[0] >10.0) return UpartialX_deg_zero(i,rho,T,X);
	else if (wne[0] >= 0.1) {
		double eta = solveEta<&fermi_sagar>(rho,T,X);
		return fd_u_del_X<&fermi_sagar>(i,eta,T,X);
	}
	else if (wne[0] <  0.1) return 1.5*rho*boltzmann_k*T*pow(X.mu_e(),-2)*chemical::partial_mu_e(i,X);
  return 0.0;
}


//pressure of electron gas, where degeneracy can be considered large
// following Pichon 1989, eq 27-29, based on Sommerfeld approximation
// can also be seen in Cox & Giuli, sec 24.2
// *  this partial pressure is less stable, due to high inaccuracies for very high eta
// *  does not smoothly limit to complete degeneracy case
double pressure_deg_large(double const rho, double const T, Abundance const& X){
	using namespace Chandrasekhar;
	using namespace FermiDirac;
	double const beta = boltzmann_k*T/(electron.mass_CGS*C_CGS*C_CGS);
	//approximate eta as in Cox and Giuli eq 24.29', pg 793 (but using 24.327b)
	double const eta_coeff = pow(8.*B0*B0, -1./3.);
	double const rho_coeff = 3.*sqrt(2.)*B0*pow(beta,1.5);
	double const eta1 = eta_coeff*pow(rho, 2./3.)/beta;
	//  finds eta by inverting rho = rho0*( F1/2(eta) + beta*F3/2(eta))
	std::function<double(double)> densDiff = [rho, rho_coeff, beta](double eta)->double {
		return rho_coeff*(FermiDiracEtaLarge1Half(eta,beta) + beta*FermiDiracEtaLarge3Half(eta,beta)) - rho;
	};
	double eta2 = eta1, eta_min = eta1*0.9, eta_max = eta1*1.1;
	rootfind::bisection_find_brackets_newton(densDiff, eta1, eta_min, eta_max);
	rootfind::bisection_search(densDiff, eta2, eta_min, eta_max);

	//return the pressure, as in Cox and Giuli eq 24.326b, pg 850
	double P =  16.*sqrt(2.)*A0*pow(beta,2.5)*(
		FermiDiracEtaLarge3Half(eta2,beta) + 0.5*beta*FermiDiracEtaLarge5Half(eta2,beta));
	return P;
}

//pressure of electon gas with small degeneracy...
//  approximated as Pdeg if larger than Pideal, ignored otherwise
//  not a great method, really, but combined with deg_zero does lead to decent models
double pressure_deg_switch(double const rho, double const T, Abundance const& X){
	double eta, X0;
	double Pdeg = pressure_deg_zero(rho,T,X);
	double Pideal = pressure_ideal(rho,T,X);
	if(Pdeg > Pideal) return Pdeg;
	else return 0.0; 
}

#endif

