// *******************************************************************
// STAR TESTS
// this will run tests of several models of star
// - isopycnic: uniform density model, simplest test case
// - polytrope: useful for many numerical tests
// - Chandrasekhar WD: 
// *******************************************************************

#include "../src/STARS/Star.h"
#include "../src/STARS/Polytrope.h"
#include "../src/STARS/ChandrasekharWD++.h"
#include "../src/MODES/ModeDriver.h"
#include "../src/MODES/NonradialModeDriver.h"
#include "../src/MODES/CowlingModeDriver.h"
#include "../src/MODES/Mode.h"
#include "../src/ThrainMode.h"   // contains Pekeris formula
// classes defined only for testing
#include "test_modes/SineModeDriver.h"
#include "test_modes/DummyMode.h" // contains SineMode
#include "test_stars/Isopycnic.h"
#include "test_stars/DummyStar.h"
#include "doctest.h"
#include <random>
#include <stdio.h>
#include <memory>

namespace {
inline void ensure_artifacts_dir() {
    std::system("mkdir -p tests/artifacts");
}
} // namespace

TEST_SUITE("Mode") {

/***** basic tests of the error functions *****/

TEST_CASE("make_any_mode") {
    std::size_t constexpr LEN(2001);
    auto uniform_star = std::make_unique<Isopycnic>(LEN);
    auto driver_on_uniform = std::make_unique<NonradialModeDriver>(uniform_star.get(), 5./3.);
    auto nonradialMode = std::make_unique<Mode<4UL>>(2, 2, 0, driver_on_uniform.get());
    CHECK_LT( 0.0, nonradialMode->getOmega2() );
    CHECK_LT( 0.0, nonradialMode->getFreq() );
    CHECK_LT( 0.0, nonradialMode->getPeriod() );
}

TEST_CASE("basic_mode_properties") {
    /** tests the following:
    *   - modes can be found with decent SSR
    *   - frequency close to known values
    *   - period, frequency, dimensionelss frequency correctly related
    *   - mode order K and moder numbers (K, L, M) correctly related
    *   - modes can be initialized with initial guess, or with brackets
    */
    std::size_t constexpr LEN(2001);
    double constexpr ACCEPTABLE_ERROR(1.e-6);
    double constexpr FREQ_ERROR(1.e-4);

    auto testStar = std::make_unique<DummyStar>(LEN);
    auto testDriver = std::make_unique<SineModeDriver>(testStar.get());

    double fourpi2 = 4.0 * m_pi * m_pi;
    double sig2omeg = testStar->Gee() * testStar->Mass() * pow(testStar->Radius(), -3);

    for(int i=1; i<10; i++){
        // create a "sine" mode -- just sine, cosine wrapped by Mode object
        auto refMode = std::make_unique<SineMode>(i, testDriver->length());
        
        // create a solved mode using the sine mode driver
        double const w2_guess = fourpi2 * i * i;
        auto testMode = std::make_unique<Mode<2UL>>(w2_guess, 2, 0, testDriver.get());
        double const ssr  = testMode->SSR();
        double const norm = testDriver->innerproduct(testMode.get(), testMode.get());
        double const overlap   = testDriver->innerproduct(testMode.get(), refMode.get());
        
        CHECK_LT(    ssr,          ACCEPTABLE_ERROR);
        CHECK_DELTA( norm,    1.0, ACCEPTABLE_ERROR);
        CHECK_DELTA( overlap, 1.0, ACCEPTABLE_ERROR);
        CHECK_DELTA(testMode->getFreq(), refMode->getFreq(), FREQ_ERROR);
        
        // verify correct relation of omega, f, and P
        CHECK_DELTA(1.0/testMode->getPeriod(), double(i), ACCEPTABLE_ERROR);
        CHECK_EQ( testMode->getPeriod(), 2.*m_pi / testMode->getFreq() );
        CHECK_EQ( sqrt(sig2omeg * testMode->getOmega2()), testMode->getFreq() );
        
        // verify correct relationship of K, L, M
        int K=0, L=0, M=0;
        testMode->modeNumbers(K, L, M);
        CHECK_EQ( testMode->modeOrder() , K );

        // reset the mode
        testMode.release();

        // test bracket init
        double const w2min = 36.0 * i * i;
        double const w2max = 41.0 * i * i;
        testMode = std::make_unique<Mode<2UL>>(w2min, w2max, 2, 0, testDriver.get());
        double const w2 = testMode->getOmega2();

        // ensure w2 is within brackets
        CHECK_LT( w2min, w2 );
        CHECK_LT( w2, w2max );

        // verify reasonable bounds on errors
        double const ssr2 = testMode->SSR();
        double const overlap2 = testDriver->innerproduct(testMode.get(), refMode.get());
        CHECK_LT(    ssr2,          ACCEPTABLE_ERROR);
        CHECK_DELTA( overlap2, 1.0, ACCEPTABLE_ERROR);
        CHECK_DELTA( testMode->getFreq(),refMode->getFreq(), FREQ_ERROR);
        
        // verify correct relation of omega, f, and P
        CHECK_DELTA(1.0/testMode->getPeriod(), double(i), ACCEPTABLE_ERROR);
        CHECK_EQ( testMode->getPeriod(), 2.*m_pi / testMode->getFreq() );
        CHECK_EQ( sqrt(sig2omeg * testMode->getOmega2()), testMode->getFreq() );
        
        // verify correct relationship of K, L, M
        testMode->modeNumbers(K, L, M);
        CHECK_EQ( testMode->modeOrder(), K );
    }
}

TEST_CASE("bad_SSR") {
    /* Tests that low number of grid points
    *  will produce a bad SSR */

    // first check an NAN is raised with fewer than 14 grid points
    std::size_t constexpr LEN1{20};
    std::size_t constexpr MIN_LEN{14};
    int const L{3}, K{5};

    auto testStar = std::make_unique<Isopycnic>(LEN1);
    auto testDriver = std::make_unique<NonradialModeDriver>(testStar.get(), 5./3.);
    auto testMode = std::make_unique<Mode<4UL>>(K, L, 0, testDriver.get());
    CHECK_LT( testDriver->length(), MIN_LEN );
    CHECK_ISNAN(testMode->SSR());
    testStar.release();
    testDriver.release();
    testMode.release();

    // make a slightly larger star and ensure it has a bad SSR
    double const BAD_ERR = 1.e-4;
    std::size_t const LEN2{41};
    testStar = std::make_unique<Isopycnic>(LEN2);
    testDriver = std::make_unique<NonradialModeDriver>(testStar.get(), 5./3.);
    testMode = std::make_unique<Mode<4UL>>(K, L, 0, testDriver.get());
    // ensure the SSR is above the threshold
    double const SSR = testMode->SSR();
    CHECK_LE( MIN_LEN, testDriver->length() );
    CHECK_LT( BAD_ERR, SSR );
    
    // ensure the c0 coefficient is above the threshold
    // must make the f-mode, verify it is f-mode, and verify not test mode
    auto fMode = std::make_unique<Mode<4UL>>(0, L, 0, testDriver.get());
    CHECK_EQ( fMode->modeOrder(), 0 );
    CHECK_NE( testMode->modeOrder(), 0 );

    // now check overlap
    double c0 = fabs(testDriver->innerproduct(testMode.get(), fMode.get()));
    CHECK_LT( BAD_ERR, c0 );
}

TEST_CASE("SSR_scale") {}

TEST_CASE("w2_freq_period") {
    std::size_t const LEN(1000);
    double constexpr GAM1(5./3.);
    double constexpr INDEX(1.5);
    double constexpr BIGM(1.9884099426e33);
    double constexpr BIGR(6.9599e10);
    double constexpr omeg2freq = G_CGS * BIGM / (BIGR * BIGR * BIGR);
    auto star = std::make_unique<Polytrope>(BIGM, BIGR, INDEX, LEN);
    auto driver = std::make_unique<CowlingModeDriver>(star.get(), GAM1);
    auto mode = std::make_unique<Mode<2UL>>(mode::JCD1_5[0][0], 1, 0, driver.get());
    // compare the dimensionless scaling factor to that from JCD-DJM
    // must convert natural frequency to angular freuency
    CHECK_DELTA( sqrt(omeg2freq) * 1.e6, 2. * m_pi * nug, 1.e-5 );
    // ensure dimensionless frequency, frequency, and period are properly related
    CHECK_EQ( sqrt(omeg2freq * mode->getOmega2()), mode->getFreq() );
    CHECK_EQ( mode->getFreq(), 2 * m_pi / mode->getPeriod() );
}

/***** TESTS OF NONRADIAL MODES *****/

TEST_CASE("pekeris_frequency") {
    double constexpr LEN{6001};
    double constexpr Gam1(5./3.);
    double constexpr ACCEPTABLE_ERROR(1.e-6);
    
    auto uniform_star = std::make_unique<Isopycnic>(LEN);
    auto driver_on_uniform = std::make_unique<NonradialModeDriver>(uniform_star.get(), Gam1);
    for(int L=1; L<4; L++){
        for(int K=1; K<10; K++){
            double const w2min = mode::calculate_Pekeris(L, K-1, Gam1);
            double const w2max = mode::calculate_Pekeris(L, K+1, Gam1);
            auto testMode = std::make_unique<Mode<4UL>>(w2min, w2max, L, 0, driver_on_uniform.get());
            int const k = testMode->modeOrder();
            CHECK_EQ( k, K );

            double const w2 = testMode->getOmega2();
            double const wPek2 = mode::calculate_Pekeris(L, k, Gam1);
            CHECK_LT( 2.0 * abs(w2-wPek2)/(w2+wPek2), ACCEPTABLE_ERROR );
            CHECK_LT( testMode->SSR(), ACCEPTABLE_ERROR );
        }
    }
}

/*
Compare to tables in Christensen-Dalsgaard and Mullan,  
Mon. Not. R. Astron. Soc. 270, 921-935 (1994)
*/

TEST_CASE("compare_JCD_DJM_1_5") {
    // Table 1, Christensen-Dalsgaard and Mullan 1994
    double const INDEX(1.5);
    std::size_t const LEN(10000);
    double const GAM1(5./3.);
    // see paragraph after eqn. 3.2 for these values
    // they give an M, but based on bad G; use GM instead
    // GM is given as GM = 1.327124448e26 dyne-cm/g
    // this code uses G = 6.67430e-8, so solving for M
    // M is given as M = 1.9884099426e33
    // R is given as R = 6.9599e10 cm
    double const BIGM(1.9884099426e33);
    double const BIGR(6.9599e10);
    // to convert from angular rad/s to microhertz
    double const CONV(1.e6 / (2.*m_pi));
    auto polytrope = std::make_unique<Polytrope>(BIGM, BIGR, INDEX, LEN);
    auto driver = std::make_unique<NonradialModeDriver>(polytrope.get(), GAM1);
    for (int L=1; L<=3; L++){
        for (int K=2; K<=14; K++){
            auto mode = std::make_unique<Mode<4UL>>(mode::omega2_JCD(INDEX, L, K-1), mode::omega2_JCD(INDEX, L, K+1), L, 0, driver.get());
            CHECK_EQ( K, mode->modeOrder() );
            CHECK_DELTA( mode->getFreq() * CONV, mode::freq_JCD(INDEX, L, K), 1.e-3 );
            CHECK_LE( mode::compare_JCD(INDEX, L, K, sqrt(mode->getOmega2())), 1.01e-4 );
            CHECK_LT( mode->SSR(), 1.e-7 );
        }
    }
}

TEST_CASE("compare_JCD_DJM_3_0") {
    // Table 2, Christensen-Dalsgaard and Mullan 1994
    double const INDEX(3.0);
    std::size_t const LEN(10000);
    double const GAM1(5./3.);
    auto polytrope = std::make_unique<Polytrope>(INDEX, LEN);
    auto driver = std::make_unique<NonradialModeDriver>(polytrope.get(), GAM1);
    for (int L=1; L<=3; L++){
        for (int K=2; K<=14; K++){
            auto mode = std::make_unique<Mode<4UL>>(mode::omega2_JCD(INDEX, L, K-1), mode::omega2_JCD(INDEX, L, K+1), L, 0, driver.get());
            CHECK_EQ( K, mode->modeOrder() );
            CHECK_LE( mode::compare_JCD(INDEX, L, K, sqrt(mode->getOmega2())), 1.01e-4 );
            CHECK_LT( mode->SSR(), 1.e-7 );
        }
    }
}

TEST_CASE("compare_JCD_DJM_4_0") {
    // Table 3, Christensen-Dalsgaard and Mullan 1994
    double const INDEX(4.0);
    std::size_t const LEN(10000);
    double const GAM1(5./3.);
    auto polytrope = std::make_unique<Polytrope>(INDEX, LEN);
    auto driver = std::make_unique<NonradialModeDriver>(polytrope.get(), GAM1);
    Mode<4UL> *mode;
    for (int L=1; L<=3; L++){
        // NOTE Osaki-Scuflaire ordering breaks down for n=4 low-order modes
        // the low-order modes in JCD-DJM use a different classification scheme for these
        for (int K=7; K<=14; K++){
            auto mode = std::make_unique<Mode<4UL>>(mode::omega2_JCD(INDEX, L, K-1), mode::omega2_JCD(INDEX, L, K+1), L, 0, driver.get());
            CHECK_EQ( K, mode->modeOrder() );
            CHECK_LT( mode::compare_JCD(INDEX, L, K, sqrt(mode->getOmega2())), 1.01e-4 );
            CHECK_LT( mode->SSR(), 1.e-7 );
        }
    }
}

/***** TESTS OF COWLING MODES *****/

// there are not as many tests for Cowling modes,
// as they are an approximation and so people 
// don't bother printing highly accurate tables of them.

TEST_CASE("same_coefficients_star") {
    std::size_t const LEN (1001);
    auto testStar = std::make_unique<Polytrope>(2.3, LEN);
    auto cowlingDriver = std::make_unique<CowlingModeDriver>(testStar.get(), 0.0);
    auto nonradialDriver= std::make_unique<NonradialModeDriver>(testStar.get(), 0.0);

    // the mode drivers sould have same radius
    // the upper 2x2 submatrix of the nonradial coefficients array
    // should be the same as the Cowling mode coefficients array
    double fakeW2 = 4.0;
    int L = 2;
    double cowlingMatrix[2UL][2UL];
    double nonradialMatrix[4UL][4UL];
    std::size_t const mode_len = cowlingDriver->length();
    CHECK_EQ( mode_len, nonradialDriver->length());
    for(std::size_t i=0; i<mode_len; i++){
        CHECK_EQ(cowlingDriver->rad(i), nonradialDriver->rad(i));
        // get the coefficient matrices and compared upper-left 2x2 block
        cowlingDriver->getCoeff(&cowlingMatrix[0][0], i, 0, fakeW2, L);
        nonradialDriver->getCoeff(&nonradialMatrix[0][0], i, 0, fakeW2, L);
        for(std::size_t j=0; j<2UL; j++ ){
            for(std::size_t k=0; k<2UL; k++){
                CHECK_EQ(cowlingMatrix[j][k], nonradialMatrix[j][k]);
            }
        }
    }
}

TEST_CASE("same_on_dummy_star") {
    /*
    It can be shown that for the DummyStar with P=rho=const,
    that there is only a single mode per L.
    In the Cowling approximation this mode has a frequency of sqrt(L).
    In the 4th-order wave form, U=3 at the surface sources a potential
    perturbation from the surface distortion, lowering the frequency to
    the Kelvin value 2L(L-1)/(2L+1); for L=1 this is zero (a neutral
    translation), which the mode search cannot converge to.
    */
    double const MACHINE_PRECISION = 1.e-15;
    std::size_t const LEN (1001);
    double const Gam1 (4./3.);
    auto testStar = std::make_unique<DummyStar>(LEN);
    auto cowlingDriver = std::make_unique<CowlingModeDriver>(testStar.get(), Gam1);
    auto nonradialDriver= std::make_unique<NonradialModeDriver>(testStar.get(), Gam1);

    // create a Cowling mode and Nonradial mode on the same star
    for (int L=1; L<10; L++){
        auto cowlingMode   = std::make_unique<Mode<2UL>>(0, L, 0, cowlingDriver.get());
        CHECK_DELTA( cowlingMode->getOmega2(), double(L), MACHINE_PRECISION );
    }
    double const KELVIN_PRECISION = 1.e-12;
    for (int L=2; L<10; L++){
        auto nonradialMode = std::make_unique<Mode<4UL>>(0, L, 0, nonradialDriver.get());
        CHECK_DELTA( nonradialMode->getOmega2(), double(2*L*(L-1))/double(2*L+1), KELVIN_PRECISION );
    }
}

} // end TEST_SUITE