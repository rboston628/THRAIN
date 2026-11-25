// *******************************************************************
// SIMPLE WD TESTS
// characterization and smoke tests for the SimpleWD model
// the golden values below were captured from the code as of PR #30,
//   before replacing the internal root-finders with lib/rootfind methods;
//   tolerances are loose enough to absorb convergence-tightness differences
// *******************************************************************

#include "doctest.h"
#include "../src/STARS/SimpleWD.h"

namespace {

struct SWDLog {
    SWDLog() {
        ThrainLogger::setOutputFile("tests/artifacts/simplewd_test_log.txt", "w");
    }
    ~SWDLog() {
        ThrainLogger::unsetOutputFile();
    }
};

} // namespace

TEST_SUITE("SimpleWD" * doctest::description("characterization and smoke tests [slow]")) {

TEST_CASE_FIXTURE(SWDLog, "characterization M=0.6 Teff=12000 N=2000"
        * doctest::timeout(30)) {
    SimpleWD star(0.6, 12000., 2000);
    CHECK(star.Mass() == doctest::Approx(0.6 * MSOLAR).epsilon(1.e-10));
    CHECK(star.Radius() == doctest::Approx(8.6785873670350862e+08).epsilon(1.e-4));
    CHECK(star.rho(0) == doctest::Approx(3.9182744702800866e+06).epsilon(1.e-4));
    CHECK(star.P(0) == doctest::Approx(2.0501636603158945e+23).epsilon(1.e-4));
    CHECK(star.SSR() < 1.e-15);
}

TEST_CASE_FIXTURE(SWDLog, "characterization M=0.8 Teff=12000 N=2000"
        * doctest::timeout(30)) {
    SimpleWD star(0.8, 12000., 2000);
    CHECK(star.Mass() == doctest::Approx(0.8 * MSOLAR).epsilon(1.e-10));
    CHECK(star.Radius() == doctest::Approx(6.9278284850234330e+08).epsilon(1.e-4));
    CHECK(star.rho(0) == doctest::Approx(1.1791646791430665e+07).epsilon(1.e-4));
    CHECK(star.P(0) == doctest::Approx(1.0369832065759450e+24).epsilon(1.e-4));
    CHECK(star.SSR() < 1.e-15);
}

TEST_CASE_FIXTURE(SWDLog, "construction smoke test" * doctest::timeout(60)) {
    double const masses[] = {0.4, 1.0};
    double const teffs[] = {10000., 20000.};
    std::size_t const lengths[] = {1000, 2001};
    for(int i = 0; i < 2; i++) {
        CAPTURE(masses[i]);
        SimpleWD star(masses[i], teffs[i], lengths[i]);
        // the constructor rounds up to odd, then expands to the half-grid
        CHECK(star.length() % 2 == 1);
        CHECK(star.length() >= lengths[i]);
        CHECK(star.Radius() > 0.0);
        CHECK(std::isfinite(star.SSR()));
        std::size_t const mid = star.length() / 2;
        CHECK(star.Gamma1(1) > 0.0);
        CHECK(star.Gamma1(mid) > 0.0);
        CHECK(star.Gamma1(star.length() - 2) > 0.0);
        // mass and radius must increase monotonically outward at spot checks
        CHECK(star.mr(mid) < star.mr(star.length() - 1));
        CHECK(star.rad(mid) < star.rad(star.length() - 1));
    }
}

TEST_CASE("parse partial pressure list") {
    std::vector<PartialPressure> const core = parsePartialPressureList("rad\tideal coul deg_zero");
    REQUIRE(core.size() == 4);
    CHECK(core[0].P == rad_gas.P);
    CHECK(core[1].P == ideal.P);
    CHECK(core[2].P == coul.P);
    CHECK(core[3].P == deg_zero.P);
    // unrecognized tokens are logged and skipped
    std::vector<PartialPressure> const partial = parsePartialPressureList("rad bogus ideal");
    REQUIRE(partial.size() == 2);
    CHECK(partial[0].P == rad_gas.P);
    CHECK(partial[1].P == ideal.P);
    CHECK(parsePartialPressureList("").empty());
}

TEST_CASE_FIXTURE(SWDLog, "construction with non-default params"
        * doctest::timeout(30)) {
    // mirrors sampleinput6.txt; goldens from a pre-refactor run of that input
    SimpleWDParams params;
    params.core_pressures = parsePartialPressureList("rad ideal coul deg_zero");
    params.atm_pressures = parsePartialPressureList("rad ideal");
    params.zy = 10.0; params.by = 3.0; params.my = 1.0;
    params.zc = 3.0;  params.bc = 4.0; params.mc = 1.0;
    params.zo = 2.0;  params.bo = 2.0; params.mo = 0.7;
    SimpleWD star(0.608, 12000., 2000, params);
    CHECK(star.Radius() == doctest::Approx(8.5989e+08).epsilon(1.e-3));
    CHECK(star.SSR() < 1.e-15);
    // the modified chemical profile must be reflected in the layer masses
    CHECK(star.Xmass.He4 == doctest::Approx(5.523e-02).epsilon(1.e-2));
    CHECK(star.Xmass.C12 == doctest::Approx(3.808e-01).epsilon(1.e-2));
    CHECK(star.Xmass.O16 == doctest::Approx(5.639e-01).epsilon(1.e-2));
}

} // TEST_SUITE
#include "../src/STARS/SimpleWD.h"
#include "doctest.h"
#include <stdio.h>
#include <cmath>

double center_expand_local(double const x, double const *const ac, int const maxPow){
  double value = ac[0];
  // ac[] here are coefficients for powers x^(2*i) stored in sequence [0,2,4] => we treat indices as 0..maxPow/2
  for(int i=1; i<=maxPow/2; ++i){
    value += ac[i]*pow(x,2*i);
  }
  return value;
}

double surface_expand_inverse_local(double const t, double const *const ac, int const maxPow){
  const int O = 1; // same convention used in SimpleWD
  double value = ac[O-1]/t;
  for(int i=0; i<=maxPow; ++i){
    value += ac[O+i]*pow(t,i);
  }
  return value;
}

class SimpleWDTest : public CxxTest::TestSuite {
public:
  static SimpleWDTest *createSuite (){
    system("mkdir -p tests/artifacts");
    ThrainLogger::setOutputFile("tests/tests/simplewd_log.txt", "w");
    printf("\n## SIMPLE WD TESTS  -- some of these are slow ##\n");
    return new SimpleWDTest();
  }
  static void destroySuite(SimpleWDTest *suite) {
    ThrainLogger::unsetOutputFile();
    printf("\n################");
    delete suite;
  }

  void test_construction_variety_and_SSR(){
    // construct several SimpleWD models and check SSR within a reasonable bound
    const std::size_t LEN = 1001;
    double masses[] = {0.2, 0.6, 1.0};
    double teffs[]  = {9000.0, 12000.0};
    const double expected_error = 10.0 * pow(1.0/LEN, 4); // same scaling heuristic used elsewhere

    for(double M : masses){
      for(double T : teffs){
        SimpleWD *s = new SimpleWD(M, T, LEN);
        double residual = pow(10.0, ceil(log10(s->SSR())));
        TS_ASSERT_LESS_THAN(residual, expected_error);
        delete s;
      }
    }
  }

  void xtest_SSR_scales_with_N(){
    std::size_t LENs[] = {251, 501, 1001};
    const int Ntest = 3;
    double SSR[Ntest];

    for(int i=0;i<Ntest;++i){
      SimpleWD *s = new SimpleWD(0.6, 12000.0, LENs[i]);
      SSR[i] = s->SSR();
      delete s;
    }
    // SSR should decrease with increasing resolution
    for(int i=1;i<Ntest;++i){
        TS_ASSERT_LESS_THAN(SSR[i], fmax(SSR[i-1], 1.e-16));
    }
  }

  void xtest_half_grid_point_behavior(){
    // Verify that constructor enforces odd internal grid length (expands even -> odd)
    std::size_t evenLen = 1000;
    std::size_t oddLen  = 1001;
    SimpleWD *s_even = new SimpleWD(0.6, 12000.0, evenLen);
    SimpleWD *s_odd  = new SimpleWD(0.6, 12000.0, oddLen);

    // final length() should equal the (possibly adjusted) "true" grid size used internally
    // SimpleWD makes even -> odd by incrementing once, so expected length equals original if odd, or +1 if even
    std::size_t len_even = s_even->length();
    std::size_t len_odd  = s_odd->length();
    TS_ASSERT_EQUALS(len_even, evenLen + 1);
    TS_ASSERT_EQUALS(len_odd,  oddLen);

    delete s_even;
    delete s_odd;
  }

  void xtest_center_and_surface_expansions_small(){
    // Ensure get*Center / get*Surface coefficients reproduce the function near center/surface
    const std::size_t LEN = 1001;
    SimpleWD *s = new SimpleWD(0.608, 12000.0, LEN);

    // Center coefficients
    constexpr int MAXPOW_C = 4;
    double Ac[3] = {0.0}, Vc[3] = {0.0}, Uc[3] = {0.0}, cc[3] = {0.0};
    int maxC = MAXPOW_C;
    s->getAstarCenter(Ac, maxC, 5./3.);
    s->getVgCenter(Vc, maxC, 5./3.);
    s->getUCenter(Uc, maxC);
    s->getC1Center(cc, maxC);

    // Test for first few points near center
    const double tol_center = 1e-4;
    std::size_t Ncheck = 4;
    double R = s->Radius();
    for(std::size_t i=0;i<Ncheck;++i){
      double x = s->rad(i)/R;
      double a_fit = center_expand_local(x, Ac, MAXPOW_C);
      double v_fit = center_expand_local(x, Vc, MAXPOW_C);
      double u_fit = center_expand_local(x, Uc, MAXPOW_C);
      double c_fit = center_expand_local(x, cc, MAXPOW_C);

      TS_ASSERT_DELTA(a_fit, s->getAstar(i, 5./3.), tol_center);
      TS_ASSERT_DELTA(v_fit, s->getVg(i, 5./3.), tol_center);
      TS_ASSERT_DELTA(u_fit, s->getU(i), tol_center);
      TS_ASSERT_DELTA(c_fit, s->getC(i), tol_center);
    }

    // Surface coefficients (inverse expansion)
    constexpr int MAXPOW_S = 3;
    int maxS = MAXPOW_S;
    double As[6] = {0.0}, Vs[6] = {0.0}, Us[5] = {0.0}, cs[5] = {0.0};
    s->getAstarSurface(As, maxS, 5./3.);
    s->getVgSurface(Vs, maxS, 5./3.);
    s->getUSurface(Us, maxS);
    s->getC1Surface(cs, maxS);

    const double tol_surf = 1e-3;
    // test a few points near surface (avoid the very last index if degenerate)
    std::size_t L = s->length();
    for(std::size_t idx = L-2; idx > L-6; --idx){
      double t = 1.0 - s->rad(idx)/s->Radius();
      double v_fit = surface_expand_inverse_local(t, Vs, MAXPOW_S);
      double a_fit = surface_expand_inverse_local(t, As, MAXPOW_S);
      TS_ASSERT_DELTA(v_fit, s->getVg(idx, 5./3.), tol_surf);
      // A* can be large near surface => compare relatively
      double refV = s->getVg(idx, 5./3.);
      if(std::abs(refV) > 0.0) {
        TS_ASSERT_LESS_THAN(fabs((a_fit - s->getAstar(idx, 5./3.))/s->getAstar(idx, 5./3.)), 1e-1);
      } else {
        TS_ASSERT_DELTA(a_fit, s->getAstar(idx, 5./3.), tol_surf);
      }
    }
    delete s;
  }
};
