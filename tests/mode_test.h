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
#include <cxxtest/TestSuite.h>
#include <random>
#include <stdio.h>

class ModeTest : public CxxTest::TestSuite {

Isopycnic *uniform_star;
NonradialModeDriver *driver_on_uniform;

public:

static ModeTest *createSuite (){
    system( "mkdir -p tests/artifacts" );
    printf("\n## MODE TESTS  -- some of these are slow ##\n");
    return new ModeTest();
}
static void destroySuite(ModeTest *suite) { 
    printf("\n################");
    delete suite; 
}

void setUp() {
    freopen("tests/artifacts/startest_log.txt", "w", stdout);
    printf("BEGIN MODE TESTS\n");
    constexpr std::size_t LEN = 6001;
    uniform_star = new Isopycnic(LEN);
    driver_on_uniform = new NonradialModeDriver(uniform_star, 5./3.);
}

void tearDown() {
    freopen("/dev/tty", "w", stdout);
}

/***** basic tests of the error functions *****/

void test_make_any_mode(){
    Mode<4UL> *nonradialMode = new Mode<4UL>(2, 2, 0, driver_on_uniform);
    TS_ASSERT_LESS_THAN( 0.0, nonradialMode->getOmega2() );
    TS_ASSERT_LESS_THAN( 0.0, nonradialMode->getFreq() );
    TS_ASSERT_LESS_THAN( 0.0, nonradialMode->getPeriod() );
}

void test_basic_mode_properties(){
    /** tests the following:
    *   - modes can be found with decent SSR
    *   - frequency close to known values
    *   - period, frequency, dimensionelss frequency correctly related
    *   - mode order K and moder numbers (K, L, M) correctly related
    *   - modes can be initialized with initial guess, or with brackets
    */
    std::size_t const LEN(2001);
    double const MACHINE_PRECISION(1.e-15);
    double const ACCEPTABLE_ERROR(1.e-6);
    double const FREQ_ERROR(1.e-4);
    Star *testStar = new DummyStar(LEN);
    ModeDriver *testDriver = new SineModeDriver(testStar);
    Mode<2UL> *testMode;
    SineMode *refMode;
    double fourpi2 = 4.0 * m_pi * m_pi;
    double w2, w2min, w2max;
    int K, L, M;
    double sig2omeg = testStar->Gee() * testStar->Mass() * pow(testStar->Radius(), -3);
    double ssr, norm, overlap;
    for(int i=1; i<10; i++){
        // create a "sine" mode -- just sine, cosine wrapped by Mode object
        refMode = new SineMode(i, testDriver->length());
        // create a solved mode using the sine mode driver
        w2 = fourpi2 * i * i;
        testMode = new Mode<2UL>(w2, 2, 0, testDriver);
        ssr  = testMode->SSR();
        norm = testDriver->innerproduct(testMode, testMode);
        overlap   = testDriver->innerproduct(testMode, refMode);
        TS_ASSERT_LESS_THAN( ssr,      ACCEPTABLE_ERROR);
        TS_ASSERT_DELTA( norm,    1.0, ACCEPTABLE_ERROR);
        TS_ASSERT_DELTA( overlap, 1.0, ACCEPTABLE_ERROR);
        TS_ASSERT_DELTA(testMode->getFreq(), refMode->getFreq(), FREQ_ERROR);
        // verify correct relation of omega, f, and P
        TS_ASSERT_DELTA(1.0/testMode->getPeriod(), double(i), ACCEPTABLE_ERROR);
        TS_ASSERT_EQUALS( testMode->getPeriod(), 2.*M_PI / testMode->getFreq() );
        TS_ASSERT_EQUALS( sqrt(sig2omeg * testMode->getOmega2()), testMode->getFreq() );
        // verify correct relationship of K, L, M
        testMode->modeNumbers(K, L, M);
        TS_ASSERT_EQUALS( testMode->modeOrder() , K );
        delete testMode;
        // test bracket init
        w2min = 36.0 * i * i;
        w2max = 41.0 * i * i;
        testMode = new Mode<2UL>(w2min, w2max, 2, 0, testDriver);
        w2 = testMode->getOmega2();
        // ensure w2 is within brackets
        TS_ASSERT_LESS_THAN( w2min, w2 );
        TS_ASSERT_LESS_THAN( w2, w2max );
        // verify reasonable bounds on errors
        ssr = testMode->SSR();
        overlap = testDriver->innerproduct(testMode, refMode);
        TS_ASSERT_LESS_THAN( ssr,      ACCEPTABLE_ERROR);
        TS_ASSERT_DELTA( overlap, 1.0, ACCEPTABLE_ERROR);
        TS_ASSERT_DELTA( testMode->getFreq(),refMode->getFreq(), FREQ_ERROR);
        // verify correct relation of omega, f, and P
        TS_ASSERT_DELTA(1.0/testMode->getPeriod(), double(i), ACCEPTABLE_ERROR);
        TS_ASSERT_EQUALS( testMode->getPeriod(), 2.*M_PI / testMode->getFreq() );
        TS_ASSERT_EQUALS( sqrt(sig2omeg * testMode->getOmega2()), testMode->getFreq() );
        // verify correct relationship of K, L, M
        testMode->modeNumbers(K, L, M);
        TS_ASSERT_EQUALS( testMode->modeOrder(), K );
        delete testMode;
        delete refMode;
    }
    delete testDriver;
    delete testStar;
}

void test_bad_SSR(){
    /* Tests that low number of grid points
    *  will produce a bad SSR */

    // first check an NAN is raised with fewer than 14 grid points
    std::size_t const LEN1{20};
    std::size_t const MIN_LEN{14};
    int const L{3}, K{5};
    Star *testStar = new Isopycnic(LEN1);
    ModeDriver *testDriver = new NonradialModeDriver(testStar, 5./3.);
    Mode<4UL> *testMode = new Mode<4UL>(K, L, 0, testDriver);
    TS_ASSERT_LESS_THAN( testDriver->length(), MIN_LEN );
    TS_ASSERT( std::isnan(testMode->SSR()) );
    delete testStar;
    delete testDriver;
    delete testMode;

    // make a slightly larger star and ensure it has a bad SSR
    double const BAD_ERR = 1.e-4;
    std::size_t const LEN2{41};
    testStar = new Isopycnic(LEN2);
    testDriver = new NonradialModeDriver(testStar, 5./3.);
    testMode = new Mode<4UL>(K, L, 0, testDriver);
   // ensure the SSR is above the threshold
    double SSR = testMode->SSR();
    TS_ASSERT_LESS_THAN_EQUALS( MIN_LEN, testDriver->length() );
    TS_ASSERT_LESS_THAN( BAD_ERR, SSR );
    // ensure the c0 coefficient is above the threshold
    // must make the f-mode, verify it is f-mode, and verify not test mode
    Mode<4UL> *fMode = new Mode<4UL>(0, L, 0, testDriver);
    TS_ASSERT_EQUALS( fMode->modeOrder(), 0 );
    TS_ASSERT_DIFFERS( testMode->modeOrder(), 0 );
    // now check overlap
    double c0 = fabs(testDriver->innerproduct(testMode, fMode));
    TS_ASSERT_LESS_THAN( BAD_ERR, c0 );
    delete testStar;
    delete testDriver;
    delete testMode;
    delete fMode;
}

void test_SSR_scale(){}

void test_w2_freq_period(){
    std::size_t const LEN(1000);
    double const GAM1(5./3.);
    double const INDEX(1.5);
    double const BIGM(1.9884099426e33);
    double const BIGR(6.9599e10);
    double const omeg2freq = G_CGS * BIGM * pow(BIGR, -3);
    Polytrope *star = new Polytrope(BIGM, BIGR, INDEX, LEN);
    CowlingModeDriver *driver = new CowlingModeDriver(star, GAM1);
    Mode<2UL> *mode = new Mode<2UL>(mode::JCD1_5[0][0], 1, 0, driver);
    // compare the dimensionless scaling factor to that from JCD-DJM
    // must convert natural frequency to angular freuency
    TS_ASSERT_DELTA( sqrt(omeg2freq) * 1.e6, 2. * M_PI * nug, 1.e-5 );
    // ensure dimensionless frequency, frequency, and period are properly related
    TS_ASSERT_EQUALS( sqrt(omeg2freq * mode->getOmega2()), mode->getFreq() );
    TS_ASSERT_EQUALS( mode->getFreq(), 2 * M_PI / mode->getPeriod() )
}

/***** TESTS OF NONRADIAL MODES *****/

void test_pekeris_frequency(){
    fprintf(stderr, "\nPEKERIS TEST: ");
    double const Gam1(5./3.);
    double const ACCEPTABLE_ERROR(1.e-6);
    Mode<4UL> *testMode;
    int k;
    double w2, wPek2, w2min, w2max;
    for(int L=1; L<4; L++){
        for(int K=1; K<10; K++){
            w2min = mode::calculate_Pekeris(L, K-1, Gam1);
            w2max = mode::calculate_Pekeris(L, K+1, Gam1);
            wPek2 = mode::calculate_Pekeris(L, K  , Gam1);
            testMode = new Mode<4UL>(w2min, w2max, L, 0, driver_on_uniform);
            k = testMode->modeOrder();
            TS_ASSERT_EQUALS( k, K );
            w2 = testMode->getOmega2();
            wPek2 = mode::calculate_Pekeris(L, k, Gam1);
            TS_ASSERT_LESS_THAN( 2.0 * abs(w2-wPek2)/(w2+wPek2), ACCEPTABLE_ERROR );
            TS_ASSERT_LESS_THAN( testMode->SSR(), ACCEPTABLE_ERROR );
            delete testMode;
        }
    }
}

/*
Compare to tables in Christensen-Dalsgaard and Mullan,  
Mon. Not. R. Astron. Soc. 270, 921-935 (1994)
*/

void test_compare_JCD_DJM_1_5(){
    // Table 1, Christensen-Dalsgaard and Mullan 1994
    
    fprintf(stderr, "\nCOMPARING TO JCD-DJM for n=1.5: ");
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
    double const CONV(1.e6 / (2.*M_PI));
    Polytrope *polytrope = new Polytrope(BIGM, BIGR, INDEX, LEN);
    NonradialModeDriver *driver = new NonradialModeDriver(polytrope, GAM1);
    Mode<4UL> *mode;
    for (int L=1; L<=3; L++){
        for (int K=2; K<=14; K++){
            mode = new Mode<4UL>(mode::omega2_JCD(INDEX, L, K-1), mode::omega2_JCD(INDEX, L, K+1), L, 0, driver);
            TS_ASSERT_EQUALS( K, mode->modeOrder() );
            TS_ASSERT_DELTA( mode->getFreq() * CONV, mode::freq_JCD(INDEX, L, K), 1.e-3 );
            TS_ASSERT_LESS_THAN_EQUALS( mode::compare_JCD(INDEX, L, K, sqrt(mode->getOmega2())), 1.01e-4 );
            TS_ASSERT_LESS_THAN( mode->SSR(), 1.e-7 );
            delete mode;
        }
    }
    delete driver;
    delete polytrope;
}

void test_compare_JCD_DJM_3_0(){
    // Table 2, Christensen-Dalsgaard and Mullan 1994
    
    fprintf(stderr, "\nCOMPARING TO JCD-DJM for n=3.0: ");
    double const INDEX(3.0);
    std::size_t const LEN(10000);
    double const GAM1(5./3.);
    Polytrope *polytrope = new Polytrope(INDEX, LEN);
    NonradialModeDriver *driver = new NonradialModeDriver(polytrope, GAM1);
    Mode<4UL> *mode;
    for (int L=1; L<=3; L++){
        for (int K=2; K<=14; K++){
            mode = new Mode<4UL>(mode::omega2_JCD(INDEX, L, K-1), mode::omega2_JCD(INDEX, L, K+1), L, 0, driver);
            TS_ASSERT_EQUALS( K, mode->modeOrder() );
            TS_ASSERT_LESS_THAN_EQUALS( mode::compare_JCD(INDEX, L, K, sqrt(mode->getOmega2())), 1.01e-4 );
            TS_ASSERT_LESS_THAN( mode->SSR(), 1.e-7 );
            delete mode;
        }
    }
    delete driver;
    delete polytrope;
}

void test_compare_JCD_DJM_4_0(){
    // Table 3, Christensen-Dalsgaard and Mullan 1994
    
    fprintf(stderr, "\nCOMPARING TO JCD-DJM for n=4.0: ");
    double const INDEX(4.0);
    std::size_t const LEN(10000);
    double const GAM1(5./3.);
    Polytrope *polytrope = new Polytrope(INDEX, LEN);
    NonradialModeDriver *driver = new NonradialModeDriver(polytrope, GAM1);
    Mode<4UL> *mode;
    for (int L=1; L<=3; L++){
        // NOTE Osaki-Scuflaire ordering breaks down for n=4 low-order modes
        // the low-order modes in JCD-DJM use a different classification scheme for these
        for (int K=7; K<=14; K++){
            mode = new Mode<4UL>(mode::omega2_JCD(INDEX, L, K-1), mode::omega2_JCD(INDEX, L, K+1), L, 0, driver);
            TS_ASSERT_EQUALS( K, mode->modeOrder() );
            TS_ASSERT_LESS_THAN( mode::compare_JCD(INDEX, L, K, sqrt(mode->getOmega2())), 1.01e-4 );
            TS_ASSERT_LESS_THAN( mode->SSR(), 1.e-7 );
            delete mode;
        }
    }
    delete driver;
    delete polytrope;
    fprintf(stderr, "\n");
}


/***** TESTS OF COWLING MODES *****/

// there are not as many tests for Cowling modes,
// as they are an approximation and so people 
// don't bother printing highly accurate tables of them.

void test_same_coefficients_star(){
    std::size_t const LEN (1001);
    Polytrope *testStar = new Polytrope(2.3, LEN);
    CowlingModeDriver *cowlingDriver = new CowlingModeDriver(testStar, 0.0);
    NonradialModeDriver *nonradialDriver= new NonradialModeDriver(testStar, 0.0);

    // the mode drivers sould have same radius
    // the upper 2x2 submatrix of the nonradial coefficients array
    // should be the same as the Cowling mode coefficients array
    double fakeW2 = 4.0;
    int L = 2;
    double cowlingMatrix[2UL][2UL];
    double nonradialMatrix[4UL][4UL];
    std::size_t const mode_len = cowlingDriver->length();
    TS_ASSERT_EQUALS( mode_len, nonradialDriver->length());
    for(std::size_t i=0; i<mode_len; i++){
        TS_ASSERT_EQUALS(cowlingDriver->rad(i), nonradialDriver->rad(i));
        // get the coefficient matrices and compared upper-left 2x2 block
        cowlingDriver->getCoeff(&cowlingMatrix[0][0], i, 0, fakeW2, L);
        nonradialDriver->getCoeff(&nonradialMatrix[0][0], i, 0, fakeW2, L);
        for(std::size_t j=0; j<2UL; j++ ){
            for(std::size_t k=0; k<2UL; k++){
                TS_ASSERT_EQUALS(cowlingMatrix[j][k], nonradialMatrix[j][k]);
            }
        }
    }
    delete testStar;
    delete cowlingDriver;
    delete nonradialDriver;
}

void test_same_on_dummy_star(){
    /*
    It can be shown that for the DummyStar with P=rho=const, 
    that there is only a single mode per L,
    that this mode has a frequency of sqrt(L),
    and that it is the same in both Cowling and 4th-order wave forms.
    */
   double const MACHINE_PRECISION = 1.e-15;
    std::size_t const LEN (1001);
    double const Gam1 (4./3.);
    DummyStar *testStar = new DummyStar(LEN);
    CowlingModeDriver *cowlingDriver = new CowlingModeDriver(testStar, Gam1);
    NonradialModeDriver *nonradialDriver= new NonradialModeDriver(testStar, Gam1);
    Mode<4UL> *nonradialMode;
    Mode<2UL> *cowlingMode;

    // create a Cowling mode and Nonradial mode on the same star
    for (int L=1; L<10; L++){
        nonradialMode = new Mode<4UL>(0,L,0, nonradialDriver);
        cowlingMode   = new Mode<2UL>(0,L,0, cowlingDriver);
        TS_ASSERT_DELTA( cowlingMode->getOmega2(), nonradialMode->getOmega2(), MACHINE_PRECISION );
        TS_ASSERT_DELTA( cowlingMode->getOmega2(), double(L), MACHINE_PRECISION );
    }
}


}; // end ModeTest