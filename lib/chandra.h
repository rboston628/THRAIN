#ifndef CHANDRACONSTH
#define CHANDRACONSTH

namespace Chandrasekhar {
	// these values are defined to be consistent with NIST CODATA 2022 values, as in constants.h
	double constexpr A0 = 6.002332e22;
	double constexpr B0 = 9.8101895e5;

	// these values taken from Chandrasekhar 1932, pg 416 in Dover reprint
	double constexpr A01939 = 6.01e22;
	double constexpr B01939 = 9.82e5;

	double factor_f(double x);
	double factor_g(double x);
	double factor_h(double x, double, double);
}

namespace FermiDirac {
	// general formulae
	double FermiDirac(double const k, double const eta);
	double FermiDirac(double const k, double const eta, double const beta);
	double FermiDiracDelEta(double const k, double const eta, double const beta);
	double FermiDiracDelBeta(double const k, double const eta, double const beta);
	
	//special cases
	double FermiDirac1Half(double const eta, double const beta);
	double FermiDirac3Half(double const eta, double const beta);
	double FermiDirac5Half(double const eta, double const beta);
	
	//partial derivatives of the special cases
	double FermiDirac1HalfDelEta(double const eta, double const beta);
	double FermiDirac3HalfDelEta(double const eta, double const beta);
	double FermiDirac5HalfDelEta(double const eta, double const beta);
	double FermiDirac1HalfDelBeta(double const eta, double const beta);
	double FermiDirac3HalfDelBeta(double const eta, double const beta);
	double FermiDirac5HalfDelBeta(double const eta, double const beta);
	
	//asymptotic limit for large degeneracy
	double FermiDiracEtaLarge1Half(double const eta);
	double FermiDiracEtaLarge3Half(double const eta);
	double FermiDiracEtaLarge5Half(double const eta);
	double FermiDiracEtaLarge1Half(double const eta, double const beta);
	double FermiDiracEtaLarge3Half(double const eta, double const beta);
	double FermiDiracEtaLarge5Half(double const eta, double const beta);
	//the little-f functions of the asymptotic limit
	// see Cox and Giuli, or Pichon 1989, eq 27-29
	double FermiDiracLittleF1Half(double const y);
	double FermiDiracLittleF3Half(double const y);
	double FermiDiracLittleF5Half(double const y);
}

#endif