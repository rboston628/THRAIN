#include "../lib/stellar.h"
#include "../src/constants.h"
#include "doctest.h"

extern PartialPressure ideal, rad_gas;
extern ElectronPressure deg_zero, deg_sagar, deg_swap, deg_integrate;

namespace {
// generalized Fermi-Dirac integral F_k(eta, beta) partial derivatives, by central difference
double FDDelEtaNumeric(double k, double eta, double beta){
    double const h = 1.e-4*std::max(1.0, fabs(eta));
    return (FermiDirac::FermiDirac(k, eta+h, beta) - FermiDirac::FermiDirac(k, eta-h, beta))/(2.*h);
}
double FDDelBetaNumeric(double k, double eta, double beta){
    double const h = 1.e-4*beta;
    return (FermiDirac::FermiDirac(k, eta, beta+h) - FermiDirac::FermiDirac(k, eta, beta-h))/(2.*h);
}
double betaOf(double T){
    return boltzmann_k*T/(electron.mass_CGS*C_CGS*C_CGS);
}
double const DENSITIES[] = {1.e-8, 1.e-4, 1.e-1, 1.e2, 1.e4, 1.e6, 1.e8};
double const TEMPERATURES[] = {1.e4, 1.e5, 1.e6, 1.e7, 1.e8};
} // namespace

TEST_SUITE("EOS [unit]") {

TEST_CASE("eos: fermi_dirac_limits") {
    /** the numerical integral must work for any eta, including eta <= 0 **/
    // non-degenerate limit, F_k(eta) -> Gamma(k+1) exp(eta)
    for(double eta : {-50.0, -20.0, -10.0}){
        CHECK_DELTA( FermiDirac::FermiDirac(0.5, eta), tgamma(1.5)*exp(eta), 1.e-4 );
        CHECK_DELTA( FermiDirac::FermiDirac(1.5, eta), tgamma(2.5)*exp(eta), 1.e-4 );
        CHECK_DELTA( FermiDirac::FermiDirac(2.5, eta), tgamma(3.5)*exp(eta), 1.e-4 );
    }
    // eta = 0 has an exact value F_1/2(0) = Gamma(3/2) (1 - 2^{-1/2}) zeta(3/2)
    CHECK_DELTA( FermiDirac::FermiDirac(0.5, 0.0), tgamma(1.5)*(1.-1./sqrt(2.))*2.6123753486854883, 1.e-6 );
    // degenerate limit, Sommerfeld expansion
    //   F_1/2(eta) = (2/3) eta^{3/2} (1 + pi^2/(8 eta^2) + 7 pi^4/(640 eta^4) + ...)
    for(double eta : {100.0, 1000.0, 1.e5}){
        CHECK_DELTA( FermiDirac::FermiDirac(0.5, eta),
            (2./3.)*pow(eta,1.5)*(1. + m_pi*m_pi/8./eta/eta + 7.*pow(m_pi,4)/640./pow(eta,4)), 1.e-8 );
    }
}

TEST_CASE("eos: fermi_dirac_sagar_matches_integral") {
    /** the Sagar quadrature, with its fallback, must agree with direct integration at any eta **/
    for(double beta : {1.e-4, 1.e-2, 1.e-1}){
        for(double eta : {-30.0, -5.0, 0.0, 3.0, 5.0, 10.0, 50.0, 300.0, 1.e4}){
            CAPTURE(eta); CAPTURE(beta);
            CHECK_DELTA( FermiDirac::FermiDirac1Half(eta,beta), FermiDirac::FermiDirac(0.5,eta,beta), 1.e-4 );
            CHECK_DELTA( FermiDirac::FermiDirac3Half(eta,beta), FermiDirac::FermiDirac(1.5,eta,beta), 1.e-4 );
            CHECK_DELTA( FermiDirac::FermiDirac5Half(eta,beta), FermiDirac::FermiDirac(2.5,eta,beta), 1.e-4 );
        }
    }
}

TEST_CASE("eos: fermi_dirac_asymptotic_matches_integral") {
    /** the large-eta (Pichon/Sommerfeld) expansions must agree with direct integration **/
    for(double beta : {0.0, 1.e-4, 1.e-2, 1.0}){
        for(double eta : {30.0, 100.0, 1.e3}){
            CAPTURE(eta); CAPTURE(beta);
            CHECK_DELTA( FermiDirac::FermiDiracEtaLarge1Half(eta,beta), FermiDirac::FermiDirac(0.5,eta,beta), 1.e-6 );
            CHECK_DELTA( FermiDirac::FermiDiracEtaLarge3Half(eta,beta), FermiDirac::FermiDirac(1.5,eta,beta), 1.e-6 );
            CHECK_DELTA( FermiDirac::FermiDiracEtaLarge5Half(eta,beta), FermiDirac::FermiDirac(2.5,eta,beta), 1.e-6 );
        }
    }
}

TEST_CASE("eos: fermi_dirac_derivatives") {
    /** the analytic partial derivatives must agree with finite differences **/
    for(double beta : {1.e-3, 1.e-1}){
        for(double eta : {-10.0, 0.0, 3.0, 10.0, 300.0}){
            CAPTURE(eta); CAPTURE(beta);
            for(double k : {0.5, 1.5, 2.5}){
                CAPTURE(k);
                CHECK_DELTA( FermiDirac::FermiDiracDelEta(k,eta,beta),  FDDelEtaNumeric(k,eta,beta),  1.e-5 );
                CHECK_DELTA( FermiDirac::FermiDiracDelBeta(k,eta,beta), FDDelBetaNumeric(k,eta,beta), 1.e-5 );
            }
            CHECK_DELTA( FermiDirac::FermiDirac1HalfDelEta(eta,beta),  FDDelEtaNumeric(0.5,eta,beta),  1.e-4 );
            CHECK_DELTA( FermiDirac::FermiDirac3HalfDelEta(eta,beta),  FDDelEtaNumeric(1.5,eta,beta),  1.e-4 );
            CHECK_DELTA( FermiDirac::FermiDirac5HalfDelEta(eta,beta),  FDDelEtaNumeric(2.5,eta,beta),  1.e-4 );
            CHECK_DELTA( FermiDirac::FermiDirac1HalfDelBeta(eta,beta), FDDelBetaNumeric(0.5,eta,beta), 1.e-4 );
            CHECK_DELTA( FermiDirac::FermiDirac3HalfDelBeta(eta,beta), FDDelBetaNumeric(1.5,eta,beta), 1.e-4 );
            CHECK_DELTA( FermiDirac::FermiDirac5HalfDelBeta(eta,beta), FDDelBetaNumeric(2.5,eta,beta), 1.e-4 );
        }
    }
}

TEST_CASE("eos: find_eta_all_regimes") {
    /** eta must be found for degenerate (eta>0) and non-degenerate (eta<0) electrons
    *   previously, this hung for eta<0 (both) and eta > ~40 (deg_finite) **/
    Abundance const X(0.0, 1.0, 0.0, 0.0);
    for(double rho : DENSITIES){
        for(double T : TEMPERATURES){
            CAPTURE(rho); CAPTURE(T);
            double const eta_int = deg_integrate.state(rho, T, X);
            double const eta_sag = deg_sagar.state(rho, T, X);
            CHECK( std::isfinite(eta_int) );
            CHECK( std::isfinite(eta_sag) );
            CHECK_DELTA( deg_integrate.density(eta_int, T, X), rho, 1.e-10 );
            CHECK_DELTA( deg_sagar.density(eta_sag, T, X), rho, 1.e-10 );
            CHECK( eta_int == doctest::Approx(eta_sag).epsilon(1.e-4).scale(1.0) );
        }
    }
}

TEST_CASE("eos: degeneracy_pressure_limits") {
    Abundance const X(0.0, 1.0, 0.0, 0.0);
    for(double rho : DENSITIES){
        for(double T : TEMPERATURES){
            CAPTURE(rho); CAPTURE(T);
            double const Pfin  = deg_sagar.P(rho, T, X);
            double const Ptrap = deg_integrate.P(rho, T, X);
            CHECK( std::isfinite(Pfin) );
            CHECK( std::isfinite(Ptrap) );
            CHECK_DELTA( Pfin, Ptrap, 1.e-4 );
        }
    }
    // strongly degenerate: approaches the T=0 Chandrasekhar pressure
    for(double rho : {1.e5, 1.e6, 1.e8}){
        CAPTURE(rho);
        CHECK_DELTA( deg_sagar.P(deg_sagar.state(rho, 1.e4, X), 1.e4, X), deg_zero.P(rho, 1.e4, X), 1.e-4 );
        CHECK_DELTA( deg_integrate.P(deg_integrate.state(rho, 1.e4, X), 1.e4, X),   deg_zero.P(rho, 1.e4, X), 1.e-4 );
    }
    // non-degenerate: approaches the ideal electron gas, n_e k T
    //   NOTE: tolerance is limited to ~1.3e-3 because the Chandrasekhar constants A0, B0
    //   are not exactly consistent with m_e c^2 / m_p used for n_e k T  (8 A0/B0 != m_e c^2/m_p)
    for(double rho : {1.e-8, 1.e-4}){
        double const T = 1.e7;
        double const Pideal = rho/proton.mass_CGS/X.mu_e()*boltzmann_k*T;
        CAPTURE(rho);
        CHECK_DELTA( deg_sagar.P(deg_sagar.state(rho, T, X), T, X), Pideal, 2.e-3 );
        CHECK_DELTA( deg_integrate.P(deg_integrate.state(rho, T, X), T, X),   Pideal, 2.e-3 );
    }
}

TEST_CASE("eos: invert_with_finite_degeneracy") {
    /** inverting the EOS for density must terminate and recover the density **/
    Abundance const X(0.0, 1.0, 0.0, 0.0);
    for(ElectronPressure const& deg : {deg_integrate, deg_sagar, deg_swap}){
        EOS eos;
        eos.push_back(deg);
        eos.push_back(ideal);
        eos.push_back(rad_gas);
        for(double rho : {1.e-4, 1.e2, 1.e6}){
            for(double T : {1.e4, 1.e6, 1.e7}){
                CAPTURE(rho); CAPTURE(T);
                double const P = eos(rho, T, X);
                double const rho_found = eos.invert(1.3*rho, P, T, X);
                CHECK_DELTA( rho_found, rho, 1.e-8 );
            }
        }
    }
}

} // end TEST_SUITE
