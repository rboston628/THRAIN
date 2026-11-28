#ifndef CHANDRACONSTCPP
#define CHANDRACONSTCPP

#include "chandra.h"
#include "./mathlib.h"
#include <cmath>

double Chandrasekhar::factor_f(double x){
	return x*(2.*x*x-3.)*sqrt(1.+x*x) + 3.*asinh(x);
}
	
double Chandrasekhar::factor_g(double x){
	return 8.*x*x*x*(sqrt(1.+x*x) - 1.0) - factor_f(x);
}
	
double Chandrasekhar::factor_h(double x, double sigma, double mue=1.0){
 	return mue*x*x*x + sigma*factor_g(x);
 }


namespace {
// composite Gauss-Legendre integration of F_k, valid for any eta
//   - the fully-occupied region x < eta-30 is smooth and uses a few wide panels in t
//   - the Fermi edge and tail, eta-30 < x < eta+60, uses panels of width 3 in x
//  this uses at most 34 panels (272 function evaluations), with relative error ~1e-7
// Similar methods may be found in
//   - Aparicio (1998, ApJS 117, 637)
//   - Gong, Zejda, Dappen & Aparicio (2001, Computer Physics Communications)
double panel_integrate(double const eta, std::function<double(double)> const &func) {
	double const edge = (eta > 0.0 ? eta : 0.0);
	double const a = (edge > 30.0 ? edge - 30.0 : 0.0);
	double const b = edge + 60.0;
	double const DX = 3.0;
	double sum = 0.0;
	if (a > 0.0) {
		double const ta = sqrt(a);
		int const NPANEL = 4;
		for(int i=0; i<NPANEL; i++)
			sum += math::gauss_legendre_8_point(ta*i/NPANEL, ta*(i+1)/NPANEL, func);
	}
	for(double x = a; x < b; x += DX){
		double const x2 = (x + DX < b ? x + DX : b);
		sum += math::gauss_legendre_8_point(sqrt(x), sqrt(x2), func);
	}
	return sum;
}

double fermi_occupation(double const e) {
	if (e > 0.0) {
		double const ex = exp(-e);
		return ex/(1.0+ex);
	}
	return 1.0/(1.0+exp(e));
}

typedef double (*G)(double x, double const s, double const f, double const e);
typedef double (*Xfactor)(double const t, double const k);

double g(double x, double const s, double const f, double const e) {
	return s * f;
}
double g_delEta(double x, double const s, double const f, double const e) {
	return s * f * fermi_occupation(-e); // d/deta f = f(1-f)
}
double g_delBeta(double x, double const s, double const f, double const e) {
	return 0.25 * x / s * f; // d/dbeta s = x/(4s)
}

double xfactor(double const t, double const k) {
	return 2.0 * pow(t, 2.0*k + 1.0);
}
double xfactor_1Half(double const t, double const k) {
	return 2.0 * t * t;
}
double xfactor_3Half(double const t, double const k) {
	return 2.0 * pow(t, 4);
}
double xfactor_5Half(double const t, double const k) {
	return 2.0 * pow(t, 6);
}

// The Fermi-Diract integrand, f(x) = x^k sqrt(1+beta*x/2) / (1+exp(x-eta))
// here beta  = kB*T/(m c^2)
// see Tooper 1969, ApJ 156
// valid for any eta (including eta <= 0)
template <G g_func, Xfactor xfactor_func>
double fermi_dirac_integrand(double const t, double const k, double const eta, double const beta) {
	double const x = t*t;
	double const s = sqrt(1.0 + 0.5*beta*x);
	double const f = fermi_occupation(x-eta);
	return xfactor_func(t, k) * g_func(x, s, f, x-eta);
}
} // anonymous namespace

/** GENERIC FERMI-DIRAC FUNCTIONS
 * These are the general Fermi-Dirac functions and integrands, using no approximations, special cases, or simplifications
 * FermiDirac = F_k(eta,beta) = int_0^infty x^k sqrt(1+beta*x/2) / (1+exp(x-eta)) dx
 */

//numerically integrate the Fermi-Dirac functon F_k(eta)
double FermiDirac::FermiDirac(double const k, double const eta) {
	auto const &integrand = [k, eta](double t) {
		return fermi_dirac_integrand<&g, &xfactor>(t, k, eta, 0.0);
	};
	return panel_integrate(eta, integrand);
}

//numerically integrate the Fermi-Dirac functon F_k(eta, beta)
double FermiDirac::FermiDirac(double const k, double const eta, double const beta) {
	auto const &integrand = [k, eta, beta](double t) {
		return fermi_dirac_integrand<&g, &xfactor>(t, k, eta, beta);
	};
	return panel_integrate(eta, integrand);
}
//partial derivatives
double FermiDirac::FermiDiracDelEta(double const k, double const eta, double const beta) {
	auto const &integrand = [k, eta, beta](double t) {
		return fermi_dirac_integrand<&g_delEta, &xfactor>(t, k, eta, beta);
	};
	return panel_integrate(eta, integrand);
}
double FermiDirac::FermiDiracDelBeta(double const k, double const eta, double const beta) {
	auto const &integrand = [k, eta, beta](double t) {
		return fermi_dirac_integrand<&g_delBeta, &xfactor>(t, k, eta, beta);
	};
	return panel_integrate(eta, integrand);
}

// special case 1/2
double FermiDirac::FermiDirac1Half(double const eta, double const beta) {
	auto const &integrand = [eta, beta](double t) {
		return fermi_dirac_integrand<&g, &xfactor_1Half>(t, 0.5, eta, beta);
	};
	return panel_integrate(eta, integrand);
}
double FermiDirac::FermiDirac1HalfDelEta(double const eta, double const beta) {
	auto const &integrand = [eta, beta](double t) {
		return fermi_dirac_integrand<&g_delEta, &xfactor_1Half>(t, 0.5, eta, beta);
	};
	return panel_integrate(eta, integrand);
}
double FermiDirac::FermiDirac1HalfDelBeta(double const eta, double const beta) {
	auto const &integrand = [eta, beta](double t) {
		return fermi_dirac_integrand<&g_delBeta, &xfactor_1Half>(t, 0.5, eta, beta);
	};
	return panel_integrate(eta, integrand);
}
// special case 3/2
double FermiDirac::FermiDirac3Half(double const eta, double const beta) {
	auto const &integrand = [eta, beta](double t) {
		return fermi_dirac_integrand<&g, &xfactor_3Half>(t, 1.5, eta, beta);
	};
	return panel_integrate(eta, integrand);
}
double FermiDirac::FermiDirac3HalfDelEta(double const eta, double const beta) {
	auto const &integrand = [eta, beta](double t) {
		return fermi_dirac_integrand<&g_delEta, &xfactor_3Half>(t, 1.5, eta, beta);
	};
	return panel_integrate(eta, integrand);
}
double FermiDirac::FermiDirac3HalfDelBeta(double const eta, double const beta) {
	auto const &integrand = [eta, beta](double t) {
		return fermi_dirac_integrand<&g_delBeta, &xfactor_3Half>(t, 1.5, eta, beta);
	};
	return panel_integrate(eta, integrand);
}
// special case 5/2
double FermiDirac::FermiDirac5Half(double const eta, double const beta) {
	auto const &integrand = [eta, beta](double t) {
		return fermi_dirac_integrand<&g, &xfactor_5Half>(t, 2.5, eta, beta);
	};
	return panel_integrate(eta, integrand);
}
double FermiDirac::FermiDirac5HalfDelEta(double const eta, double const beta) {
	auto const &integrand = [eta, beta](double t) {
		return fermi_dirac_integrand<&g_delEta, &xfactor_5Half>(t, 2.5, eta, beta);
	};
	return panel_integrate(eta, integrand);
}
double FermiDirac::FermiDirac5HalfDelBeta(double const eta, double const beta) {
	auto const &integrand = [eta, beta](double t) {
		return fermi_dirac_integrand<&g_delBeta, &xfactor_5Half>(t, 2.5, eta, beta);
	};
	return panel_integrate(eta, integrand);
}	

/** SAGAR 1991 QUADRATURE METHODS
 * These are special cases of the Fermi-Dirac functions using special quadrature method of Sagar (1991)
 * Because they use only 24 points, they are much faster inside their area of accuracy
 * These methods are reported by Sagar as acurate for -1 < eta < 10,
 * They have been independently verified for -100 < eta < 5 up to ~2e-5, with error on order 1e-4 at eta = 10
 * For regions eta > 10, Sagar suggests using asymptotic approximation
 * Outside accuracy range, the direct Fermi-Dirac integration is used instead
 * See papers:
 * - Robin P Sagar, Computer Physics Communications vol 66, 1991
 * - Robin P Sagar, ApJ 376, 1991
*/

namespace Sagar {
// abscissae (x) and weights (w) for 24-point Gauss quadrature method due to Sagar 1991
double constexpr xFermiDirac1Half[] = {
	1.065765457622244303e-01, 4.230728669515992804e-01, 9.419850084605768359e-01, 
	1.657367156890049269e+00, 2.570213091710842203e+00, 3.687909918000632621e+00, 
	5.019454072639460885e+00, 6.572842000947019917e+00, 8.355768265288956853e+00, 
	1.037702320147209848e+01, 1.264732820630387383e+01, 1.517981042473563459e+01, 
	1.799049833506082646e+01, 2.109902949629965248e+01, 2.452968403197510593e+01, 
	2.831290225833495608e+01, 3.248757170669172601e+01, 3.710461553062605055e+01, 
	4.223293075368851131e+01, 4.796992232165839131e+01, 5.446202083707815918e+01, 
	6.195032240271386363e+01, 7.089542237979413787e+01, 8.247940908146821374e+01
};
double constexpr wFermiDirac1Half[] = {
	3.285175403357927574e-02, 1.078791300687585500e-01, 1.681544217004087412e-01, 
	1.676417853951022623e-01, 1.155215546646083391e-01, 5.734703122687110474e-02, 
	2.119698568082612045e-02, 5.964995125395940070e-03, 1.290895483684744857e-03, 
	2.150426470522999493e-04, 2.743849918141434060e-05, 2.657145441059904433e-06, 
	1.927749118726048725e-07, 1.029994216793772317e-08, 3.964542390084040303e-10, 
	1.068615276190246333e-11, 1.944313861292386839e-13, 2.274868576285609567e-15, 
	1.602346574894849982e-17, 6.189182148247560184e-20, 1.139569672899732326e-22, 
	7.949390158581847580e-26, 1.359843844790543011e-29, 1.905654458662088918e-34
};
double constexpr xFermiDirac3Half[] = {
	2.129716612278677430e-01, 6.237268702986968508e-01, 1.229403044783669869e+00, 
	2.026961451899495619e+00, 3.020610609609100410e+00, 4.218734197248799002e+00, 
	5.629623967479334623e+00, 7.260592688638248677e+00, 9.119217021840780670e+00, 
	1.121447293177990263e+01, 1.355732211119446176e+01, 1.616111273495549703e+01, 
	1.904207813911802493e+01, 2.222005924434693543e+01, 2.571955002541875886e+01, 
	2.957122636257431200e+01, 3.381424706901884974e+01, 3.849986297339433161e+01, 
	4.369738885284278659e+01, 4.950479735755637819e+01, 5.606935454831911818e+01, 
	6.363352566233777885e+01, 7.266054935122159627e+01, 8.433947217510982230e+01
};
double constexpr wFermiDirac3Half[] = {
	1.365937958895298203e-02, 8.753560954084596397e-02, 2.163878832948606529e-01,
	3.004318314223437617e-01, 2.671800670589929453e-01, 1.637863466848732099e-01,
	7.261338055187011266e-02, 2.395111081211760715e-02, 5.957917543163391603e-03,
	1.122153727092106363e-03, 1.596764308326631934e-04, 1.704826903772844090e-05,
	1.350657256618920433e-06, 7.817363319325411844e-08, 3.237319698904724569e-09,
	9.333854978405723236e-11, 1.807632052729200003e-12, 2.241763132816884723e-14,
	1.667878325117133144e-16, 6.785615480105566540e-19, 1.313183740630591210e-21,
	9.616279005133126237e-25, 1.727051369706408515e-28, 2.549456047069159218e-33
};
double constexpr xFermiDirac5Half[] = {
	3.422156821319662829e-01, 8.434644632641465284e-01, 1.532396544460638270e+00, 
	2.409974354200374472e+00, 3.482882831159288010e+00, 4.759639243863683690e+00, 
	6.247759918994630441e+00, 7.954146727611163367e+00, 9.886421301541225215e+00, 
	1.205377252717992426e+01, 1.446738372295618009e+01, 1.714080389169575934e+01, 
	2.009045683883463495e+01, 2.333637615550461493e+01, 2.690326073526996817e+01, 
	3.082201391097157481e+01, 3.513205798444946900e+01, 3.988496357316911226e+01, 
	4.515045472834709800e+01, 5.102706194618288853e+01, 5.766287403453427061e+01, 
	6.530170651996465864e+01, 7.440940721694892776e+01, 8.618189363225279673e+01
};
double constexpr wFermiDirac5Half[] = {
	1.156875083599599910e-02, 1.170293449794717130e-01, 4.040945556823653779e-01, 
	7.236312791373911478e-01, 7.917373471518581726e-01, 5.799992638904753491e-01, 
	3.007177501667138216e-01, 1.138861544833988782e-01, 3.201210622537985194e-02, 
	6.721349241397831141e-03, 1.054098130731385803e-03, 1.228570478281638447e-04, 
	1.053984817271198724e-05, 6.560497130547129842e-07, 2.904716849877691045e-08, 
	8.909289940219604541e-10, 1.827655298772768217e-11, 2.392185495248132500e-13, 
	1.872698757697392744e-15, 7.996932286473866299e-18, 1.621431452583757569e-20, 
	1.242727249852022253e-23, 2.336588459348587113e-27, 3.622936475852263175e-32
};

// independently tested range of eta where the Sagar 1991 quadrature
// accurate to better than ~2e-5
double const SAGAR_ETA_MIN = -100.0;
double const SAGAR_ETA_MAX =    5.0;
bool is_valid(double const eta){
	return (SAGAR_ETA_MIN <= eta) && (eta <= SAGAR_ETA_MAX);
}
} // namespace Sagar

double FermiDirac::FermiDiracSagar1Half(double const eta, double const beta){
	if (!Sagar::is_valid(eta)) return FermiDirac1Half(eta, beta);
	double f;
	double x, w;
	double sum = 0.0;
	for(int i=0; i<24; i++){
		x = Sagar::xFermiDirac1Half[i];
		w = Sagar::wFermiDirac1Half[i];
		f = (1.+exp(x))/(1.+exp(x-eta))*sqrt(1.+0.5*beta*x);
		sum += w*f;
	}
	return sum;
}
double FermiDirac::FermiDiracSagar3Half(double const eta, double const beta){
	if (!Sagar::is_valid(eta)) return FermiDirac3Half(eta, beta);
	double f;
	double x, w;
	double sum = 0.0;
	for(int i=0; i<24; i++){
		x = Sagar::xFermiDirac3Half[i];
		w = Sagar::wFermiDirac3Half[i];
		f = (1.+exp(x))/(1.+exp(x-eta))*sqrt(1.+0.5*beta*x);
		sum += w*f;
	}
	return sum;
}
double FermiDirac::FermiDiracSagar5Half(double const eta, double const beta) {
	if (!Sagar::is_valid(eta)) return FermiDirac5Half(eta, beta);
	double f;
	double x, w;
	double sum = 0.0;
	for(int i=0; i<24; i++){
		x = Sagar::xFermiDirac5Half[i];
		w = Sagar::wFermiDirac5Half[i];
		f = fermi_occupation(x-eta) * (1.+exp(x)) * sqrt(1.+0.5*beta*x);
		sum += w*f;
	}
	return sum;
}

double FermiDirac::FermiDiracSagar1HalfDelEta(double const eta, double const beta) {
	if (!Sagar::is_valid(eta)) return FermiDirac1HalfDelEta(eta, beta);
	double f;
	double x, w;
	double sum = 0.0;
	for(int i=0; i<24; i++){
		x = Sagar::xFermiDirac1Half[i];
		w = Sagar::wFermiDirac1Half[i];
		f = (1.+exp(x))/pow(1.+exp(x-eta),2)*sqrt(1.+0.5*beta*x)*exp(x-eta);
		sum += w*f;
	}
	return sum;
}
double FermiDirac::FermiDiracSagar3HalfDelEta(double const eta, double const beta) {
	if (!Sagar::is_valid(eta)) return FermiDirac3HalfDelEta(eta, beta);
	double f;
	double x, w;
	double sum = 0.0;
	for(int i=0; i<24; i++){
		x = Sagar::xFermiDirac3Half[i];
		w = Sagar::wFermiDirac3Half[i];
		f = (1.+exp(x))/pow(1.+exp(x-eta),2)*sqrt(1.+0.5*beta*x)*exp(x-eta);
		sum += w*f;
	}
	return sum;
}
double FermiDirac::FermiDiracSagar5HalfDelEta(double const eta, double const beta) {
	if (!Sagar::is_valid(eta)) return FermiDirac5HalfDelEta(eta, beta);
	double f;
	double x, w;
	double sum = 0.0;
	for(int i=0; i<24; i++){
		x = Sagar::xFermiDirac5Half[i];
		w = Sagar::wFermiDirac5Half[i];
		f = (1.+exp(x))/pow(1.+exp(x-eta),2)*sqrt(1.+0.5*beta*x)*exp(x-eta);
		sum += w*f;
	}
	return sum;
}
// note d/dbeta sqrt(1+beta*x/2) = 0.25*x/sqrt(1+beta*x/2)
double FermiDirac::FermiDiracSagar1HalfDelBeta(double const eta, double const beta) {
	if (!Sagar::is_valid(eta)) return FermiDirac1HalfDelBeta(eta, beta);
	double f;
	double x, w;
	double sum = 0.0;
	for(int i=0; i<24; i++){
		x = Sagar::xFermiDirac1Half[i];
		w = Sagar::wFermiDirac1Half[i];
		f = (1.+exp(x))/(1.+exp(x-eta))/sqrt(1.+0.5*beta*x)*0.25*x;
		sum += w*f;
	}
	return sum;
}
double FermiDirac::FermiDiracSagar3HalfDelBeta(double const eta, double const beta) {
	if (!Sagar::is_valid(eta)) return FermiDirac3HalfDelBeta(eta, beta);
	double f;
	double x, w;
	double sum = 0.0;
	for(int i=0; i<24; i++){
		x = Sagar::xFermiDirac3Half[i];
		w = Sagar::wFermiDirac3Half[i];
		f = (1.+exp(x))/(1.+exp(x-eta))/sqrt(1.+0.5*beta*x)*0.25*x;
		sum += w*f;
	}
	return sum;
}
double FermiDirac::FermiDiracSagar5HalfDelBeta(double const eta, double const beta) {
	if (!Sagar::is_valid(eta)) return FermiDirac5HalfDelBeta(eta, beta);
	double f;
	double x, w;
	double sum = 0.0;
	for(int i=0; i<24; i++){
		x = Sagar::xFermiDirac5Half[i];
		w = Sagar::wFermiDirac5Half[i];
		f = (1.+exp(x))/(1.+exp(x-eta))/sqrt(1.+0.5*beta*x)*0.25*x;
		sum += w*f;
	}
	return sum;
}

/** ASYMPTOTIC APPROXIMATIONS
 * For very large eta values, the following asymptotic methods will work.
 * See Cox and Giuli 1980, discussion on Sommerfeld's asymptotic expansion
 * See also Pichon 1989, eq 27-28
 */

double FermiDirac::FermiDiracEtaLarge1Half(double const eta) {
	double const e2 = math::PI*math::PI / (eta * eta);
	double const e4 = e2 * e2;
	return (2./3.)*pow(eta, 1.5) * (1. + e2/8. + 7.*e4/640.);
}	
double FermiDirac::FermiDiracEtaLarge3Half(double const eta) {
	double const e2 = math::PI*math::PI / (eta * eta);
	double const e4 = e2 * e2;
	return (2./5.)*pow(eta, 2.5) * (1. + 5.*e2/8. - 7.*e4/384.);
}	
double FermiDirac::FermiDiracEtaLarge5Half(double const eta) {
	double const e2 = math::PI*math::PI / (eta * eta);
	double const e4 = e2 * e2;
	return (2./7.)*pow(eta, 3.5) * (1. + 35.*e2/24. + 49.*e4/384.);
}	

//see Pichon 1989, eq 27-28
// Sommerfeld expansion in kinetic energy eps = beta*x, total energy E = 1+eps, momentum p = sqrt(E^2-1):
//   F_k sqrt(2 beta^{2k+2}) = Phi + (pi^2 beta^2/6) phi' + (7 pi^4 beta^4/360) phi''' + (31 pi^6 beta^6/15120) phi^(5)
//   with phi(E) = (E-1)^{k-1/2} p, Phi = int_1^E phi dE, all evaluated at E = 1+u, u = eta*beta, y = p(E)
//  C1, C2, C3 contain the prefactors and powers of 1/y; the remaining polynomials in u come from phi^(n)
double FermiDirac::FermiDiracEtaLarge1Half(double const eta, double const beta) {
	if (beta == 0.0) return FermiDiracEtaLarge1Half(eta);
	double const y = sqrt( pow(1. + eta*beta, 2.) - 1.);
	double const C1 = pow(math::PI*beta,2)/6./y;
	double const C2 = pow(math::PI*beta/y/y,2)*C1*7./20.;
	double const C3 = pow(math::PI*beta/y/y,4)*C1*31./168.;
	return (FermiDiracLittleF1Half(y) 
			+ (1.+eta*beta)*(C1 + C2 + C3*(4.*y*y+7.)))/sqrt(2.*beta*beta*beta);
}
double FermiDirac::FermiDiracEtaLarge3Half(double const eta, double const beta) {
	if (beta == 0.0) return FermiDiracEtaLarge3Half(eta);
	double const y = sqrt( pow(1. + eta*beta, 2.) - 1.);
	double const C1 = pow(math::PI*beta,2)/6./y;
	double const C2 = pow(math::PI*beta/y/y,2)*C1*7./20.;
	double const C3 = pow(math::PI*beta/y/y,4)*C1*31./168.;
	double const u = eta*beta;
	return (FermiDiracLittleF3Half(y)
			+ u*(C1*(3.+2.*u) - C2 - C3*((4.*u+6.)*u + 3.)))/sqrt(2.*pow(beta,5));
}
double FermiDirac::FermiDiracEtaLarge5Half(double const eta, double const beta) {
	if (beta == 0.0) return FermiDiracEtaLarge5Half(eta);
	double const y = sqrt( pow(1. + eta*beta, 2.) - 1.);
	double const C1 = pow(math::PI*beta,2)/6./y;
	double const C2 = pow(math::PI*beta/y/y,2)*C1*7./20.;
	double const C3 = pow(math::PI*beta/y/y,4)*C1*31./168.;
	double const u = eta*beta;
	return (FermiDiracLittleF5Half(y)
			+ u*u*(C1*(5.+3.*u) + C2*(((2.*u+10.)*u + 15.)*u + 5.) + C3*(3.+5.*u)))/sqrt(2.*pow(beta,7));
}

// these functions appear in asymptotic expansion (Pichon 1989, Cox&Giuli 1980)
double FermiDirac::FermiDiracLittleF1Half(double const y) {
	return 0.5*(y*sqrt(1.+y*y) - asinh(y));
}
double FermiDirac::FermiDiracLittleF3Half(double const y) {
	return 1./3.*y*y*y - FermiDiracLittleF1Half(y);
}
double FermiDirac::FermiDiracLittleF5Half(double const y) {
	return 5./8.*y*(1.+0.4*y*y)*sqrt(1.+y*y) - 5./8.*asinh(y) - 2./3.*y*y*y;
}

#endif