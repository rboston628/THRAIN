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
