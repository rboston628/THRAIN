// *******************************************************************
// STAR TESTS
// this will run tests of several models of star
// - isopycnic: uniform density model, simplest test case
// - polytrope: useful for many numerical tests
// - Chandrasekhar WD: 
// *******************************************************************

#include "../src/STARS/Star.h"
#include "../src/STARS/Isopycnic.h"
#include "../src/STARS/Polytrope.h"
#include "../src/STARS/ChandrasekharWD++.h"
#include "../src/MODES/ModeDriver.h"
#include "../src/MODES/NonradialModeDriver.h"
#include "../src/MODES/CowlingModeDriver.h"
#include "../src/MODES/Mode.h"
#include "../src/ThrainMode.h"
#include <cxxtest/TestSuite.h>
#include <random>
#include <stdio.h>

#define ASSERT_REL_DIFFERENCE(x, y, err) { TS_ASSERT_LESS_THAN(relative_difference(x,y), err) }

// double relative_difference(double x, double y){
//     double denom = (x*y==0.0 ? 1.0 : std::abs(x+y));
//     return 2.0 * std::abs(x-y) / denom;
// }

class DummyStar : public Isopycnic {
    // this is a "star" that has constant density and pressure
    // for testing modes.  For this purpose, we only need to 
    // specify that Astar, Vg are zero everywhere.
    // This acts to decouple y1,y2 from y3,y4
private:
    std::size_t len;
    std::size_t indexFit;
public:
    DummyStar(std::size_t L) : len(L), indexFit(L/2), Isopycnic(L) {};
    double rad(std::size_t X) override {return double(X)/double(len-1);};
    std::size_t length() override {return len;};
    
    // throughout star, Astar, Vg are zero everywhere
    double getAstar(std::size_t, double GamPert=0.0) override {return 0.0;};
    double getVg(   std::size_t, double GamPert=0.0) override {return 0.0;};
    // at the boundaries, Astar, Vg coefficients are zero to all orders
    void getAstarCenter(double *Ac, int& maxPow, double g=0) override {
        for (int i=0; i<=maxPow; i++) Ac[i] = 0.0;
    };
    void getAstarSurface(double *As, int& maxPow, double g=0) override{
        getAstarCenter(As, maxPow, g);
    };
    void getVgCenter(double *Vc, int& maxPow, double g=0) override{
        getAstarCenter(Vc, maxPow, g);
    };
    void getVgSurface(double *Vs, int& maxPow, double g=0) override{
        getAstarCenter(Vs, maxPow, g);
    };
};

class ModeTest : public CxxTest::TestSuite {

Isopycnic *uniform_star;
NonradialModeDriver *driver_on_uniform;

public:

static ModeTest *createSuite (){
    system( "mkdir -p tests/artifacts" );
    printf("\n## MODE TESTS ##\n");
    return new ModeTest();
}
static void destroySuite(ModeTest *suite) { 
    printf("\n################");
    delete suite; 
}

void setUp() {
    freopen("tests/artifacts/startest_log.txt", "w", stdout);
    printf("BEGIN MODE TESTS\n");
    constexpr std::size_t LEN = 1001;
    uniform_star = new Isopycnic(LEN);
    driver_on_uniform = new NonradialModeDriver(uniform_star, 5./3.);
}

void tearDown() {
    freopen("/dev/tty", "w", stdout);
}

/***** basic tests of the error functions *****/

void test_make_anything(){
    fprintf(stderr, "START\n");
    Mode<4UL> *nonradialMode = new Mode<4UL>(2, 2, 0, driver_on_uniform);
    fprintf(stderr, "END\n");
}

void test_bad_SSR(){
    /* Tests that low number of grid points
    *  will produce a bad SSR*/

    // first check an NAN is reaised with fewer than 14 grid points
    std::size_t const LEN1{20};
    Star *testStar = new Isopycnic(LEN1);
    ModeDriver *testDriver = new NonradialModeDriver(testStar, 0.0);
    Mode<4UL> *testMode = new Mode<4UL>(2, 2, 0, testDriver);
    TS_ASSERT_LESS_THAN( testDriver->length(), 14 );
    TS_ASSERT_IS_NAN( testMode->SSR() );
    delete testStar;
    delete testDriver;
    delete testMode;

    // make a slihtly larger star and ensure it has a bad SSR
    double const BAD_ERR = 1.e-2;
    std::size_t const LEN2{40};
    testStar = new Isopycnic(LEN2);
    testDriver = new NonradialModeDriver(testStar, 0.0);
    testMode = new Mode<4UL>(2, 2, 0, testDriver);
   // ensure the SSR is above the threshold
    double SSR = testMode->SSR();
    TS_ASSERT_LESS_THAN_EQUALS(14, testDriver->length() );
    TS_ASSERT_LESS_THAN( BAD_ERR, SSR );
    // ensure the c0 coefficient is above the threshold
    // must make the f-mode, verify it is f-mode, and verify not test mode
    Mode<4UL> *fMode = new Mode<4UL>(0, 2, 0, testDriver);
    TS_ASSERT_EQUALS( fMode->modeOrder(), 0 );
    TS_ASSERT_DIFFERS( testMode->modeOrder(), 0 );
    double c0 = testDriver->innerproduct(testMode, fMode);
    TS_ASSERT_LESS_THAN( BAD_ERR, c0 );
    delete testStar;
    delete testDriver;
    delete testMode;
}

void test_SSR_scale(){}

/***** TESTS OF NONRADIAL MODES *****/

void test_pekeris_frequency(){
    double const Gam1 = 5./3.;
    Mode<4UL> *testMode;
    for(int L=1; L<4; L++){
        for(int K=1; K<15; K++){
            double w2min = mode::calculate_Pekeris(L, K-1, Gam1);
            double w2max = mode::calculate_Pekeris(L, K+1, Gam1);
            double wPek2 = mode::calculate_Pekeris(L, K  , Gam1);
            TS_ASSERT_LESS_THAN(w2min, wPek2);
            TS_ASSERT_LESS_THAN(wPek2, w2max);
            testMode = new Mode<4UL>(w2min, w2max, L, 0, driver_on_uniform);
            int k = testMode->modeOrder();
            TS_ASSERT_LESS_THAN(K-1, k);
            TS_ASSERT_LESS_THAN(k, K+1)
            TS_ASSERT_EQUALS( k, K );
            double w2 = testMode->getOmega2();
            TS_ASSERT_LESS_THAN(w2min, w2);
            TS_ASSERT_LESS_THAN(w2, w2max);
            // double wPek2 = mode::calculate_Pekeris(L, k, Gam1);
            // fprintf(stderr, "freq %d %le %le\n", k, w2, mode::compare_Pekeris(sqrt(w2), L, k, Gam1));
            TS_ASSERT_DELTA(w2, mode::calculate_Pekeris(L, k, Gam1), 1.e-5);
            TS_ASSERT_LESS_THAN(testMode->SSR(), 1.e-6);
            delete testMode;
        }
    }
}


/***** TESTS OF COWLING MODES *****/

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
        // get the coefficient matrices andc compared upper-left 2x2 block
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

// void test_same_on_dummy_star(){
//     std::size_t const LEN (1001);
//     double const Gam1 (5./3.);
//     Isopycnic *testStar = new Isopycnic(LEN);
//     CowlingModeDriver *cowlingDriver = new CowlingModeDriver(testStar, Gam1);
//     NonradialModeDriver *nonradialDriver= new NonradialModeDriver(testStar, Gam1);

//     // create a Cowling mode and Nonradial mode on the same star
//     Mode<4UL> *nonradialMode = new Mode<4UL>(1,1,0, nonradialDriver);
//     Mode<2UL> *cowlingMode = new Mode<2UL>(1,1,0, cowlingDriver);
//     fprintf(stderr, "MODES %d %le : %d %le\n", 
//         cowlingMode->modeOrder(), cowlingMode->getOmega2(),
//         nonradialMode->modeOrder(), nonradialMode->getOmega2()
//     );
//     // make sure the modes are comparable
//     TS_ASSERT_EQUALS( cowlingMode->modeOrder(), nonradialMode->modeOrder() );

// }


}; // end ModeTest'[]