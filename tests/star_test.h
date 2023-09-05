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
#include <cxxtest/TestSuite.h>
#include <random>

class StarTest : public CxxTest::TestSuite {
public:

void test_bad_SSR(){
    /* Tests that low number of grid points
    *  will produce a bad SSR*/
    std::size_t const LEN = 10;
    Star *testStar = new Polytrope(1.0,LEN);
    // test SSR
    TS_ASSERT_LESS_THAN( 1.e-6, testStar->SSR() );
    
    delete testStar;
}

void test_scale_SSR(){
    /* Test that SSR decreases up to minimum as gridsize increases*/
    std::size_t const LEN[] = {1000, 2000, 4000, 8000};
    const int Ntest = sizeof(LEN)/sizeof(std::size_t);
    double SSR[Ntest];
    Star* testStar;
    for(int i=0; i<4; i++){
        testStar = new Polytrope(1.0, LEN[i]);
        SSR[i] = testStar->SSR();
        delete testStar;
    }
    
    // test SSRs grow smaller with grid size
    for(int i=1; i<Ntest; i++){
        TS_ASSERT_LESS_THAN(SSR[i], fmax(SSR[i-1], 1.e-14));
    }
    
    // further test that it SCALES like LEN^{-4}
    // if it scaled exactly like LEN^{-4}, this would equal 2^4=16
    // because it reaches limits of machine precision, will be smaller
    double scale_factor = (SSR[1]-SSR[2])/(SSR[2]-SSR[3]);
    TS_ASSERT_LESS_THAN( scale_factor , 16.0);
    // make sure it at least scales this well
    TS_ASSERT_LESS_THAN( 6.0, scale_factor );
}

void test_uniform_star(){
    /* This test really only ensures that the SSR will be 
    ** approximately zero in the case of an exactly known star.  
    ** Because it happens  that g scales linearly with radius, it 
    ** happens to work out that the numerical derivatives in the 
    ** SSR are exactly correct*/
    printf("STAR TESTS ISOPYCNIC\n");
    const std::size_t LEN = 100;
    Star *testStar = new Isopycnic(LEN);

    const double dx = sqrt(6.0)/double(LEN-1);
    const double Rn = sqrt( 1.0/(4.0*m_pi) );

    // test SSR -- this needs to be near-zero
    TS_ASSERT_LESS_THAN( testStar->SSR(), 1.e-15 );

    // test all throughout star
    // this verifies that the getter methods work as intended
    double x=0;
    for(std::size_t i=0; i<LEN; i++){
        x = dx*i;
        TS_ASSERT_EQUALS( testStar->rho(i), 1.0 );
        TS_ASSERT_DELTA( testStar->rad(i), Rn*x, 1.e-15 );
        TS_ASSERT_DELTA( testStar->P(i), 1.0 - x*x/6.0, 1.e-15 );
        TS_ASSERT_DELTA( testStar->mr(i), Rn*(x*x*x)/3.0, 1.e-15 );
    }

    // test that mass and radius equal to known quantities
    TS_ASSERT_DELTA( testStar->Radius(), Rn*sqrt(6.0), 1.e-15 );
    TS_ASSERT_DELTA( testStar->Mass(), Rn*2.0*sqrt(6.0), 1.e-15 )
    delete testStar;
}

void test_uniform_boundaries(){
    /* this will test that the center and surface expansions
    ** in this known case act as good approxmations to the functions*/

    constexpr std::size_t LEN = 1000;
    constexpr double GAM1 = 5./3.;
    Star *testStar = new Isopycnic(LEN);

    /* test center expansions */
    constexpr std::size_t NC=10;
    constexpr int Ncent = 3;
    int maxPow = 2*(Ncent-1); // the integer passed needs to be mutable
    double Ac[Ncent], Vc[Ncent], Uc[Ncent], cc[Ncent];
    
    testStar->getAstarCenter(Ac, maxPow, GAM1);
    testStar->getVgCenter(Vc, maxPow, GAM1);
    testStar->getUCenter(Uc, maxPow);
    testStar->getC1Center(cc, maxPow);

    // simple lambda to evaluate the central expansion
    std::function<double(double,double[Ncent])> center = [Ncent](double x, double a[Ncent])->double {
        double value = a[0];
        for(int i=1; i<Ncent; i++){
            value += a[i]*pow(x,2*i);
        }
        return value;
    };

    double x=0.0, R=testStar->Radius();
    for(std::size_t i=0; i<NC; i++){
        x = testStar->rad(i)/R;
        TS_ASSERT_EQUALS( center(x, Uc), testStar->getU(i) );
        TS_ASSERT_EQUALS( center(x, cc), testStar->getC(i) );
        TS_ASSERT_DELTA( center(x, Vc), testStar->getVg(i,GAM1), 1.e-12 );
        TS_ASSERT_DELTA( center(x, Ac), testStar->getAstar(i,GAM1), 1.e-12 );
    }

    /* test surface expansions */
    constexpr std::size_t NS=10;
    constexpr int Nsurf = 6;
    maxPow = 4; // the integer passed needs to be mutable
    double As[Nsurf], Vs[Nsurf], Us[Nsurf], cs[Nsurf];
    testStar->getAstarSurface(As, maxPow, GAM1);
    testStar->getVgSurface(Vs, maxPow, GAM1);
    testStar->getUSurface(Us, maxPow);
    testStar->getC1Surface(cs, maxPow);

    // simple lambda to evaluate the central expansion
    std::function<double(double,double[Nsurf])> surface = [Nsurf](double t, double a[Nsurf])->double {
        const int O=1;
        double value = a[O-1]/t;
        value += a[O];
        for(int i=1; i<5; i++){
            value += a[i+O]*pow(t,i);
        }
        return value;
    };

    double t=0.0;
    // TODO why are these so much less accurate?
    for(std::size_t i=LEN-2; i>=LEN-1-NS; i--){
        t = 1. - testStar->rad(i)/R;
        TS_ASSERT_DELTA( surface(t, Vs), testStar->getVg(i,GAM1), 1.e-12 );
        TS_ASSERT_DELTA( surface(t, As),  testStar->getAstar(i,GAM1), 1.e-12 );
    }

    delete testStar;
}


void test_polytrope_n0(){
    printf("STAR TESTS n=0.0 POLYTROPE\n");
    const double INDEX = 0.0;
    const std::size_t LEN = 1000;
    Star *uniform = new Isopycnic(LEN);
    Star *testStar = new Polytrope(INDEX,LEN);

    // use the SSR of the polytrope to quantify acceptable error
    const double residual = ceil(testStar->SSR());

    // test they have same mass, radius, up to error
    TS_ASSERT_DELTA( uniform->Radius(), testStar->Radius(), residual );
    TS_ASSERT_DELTA( uniform->Mass(),   testStar->Mass(),   residual );

    // test key features identical inside, up to error
    // the n=0 is constant density,
    double const density = testStar->rho(0);
    for(std::size_t i=0; i<LEN; i++){
        TS_ASSERT_DELTA( uniform->rad(i), testStar->rad(i), residual );
        TS_ASSERT_DELTA( uniform->rho(i), testStar->rho(i), residual );
        TS_ASSERT_DELTA( density,         testStar->rho(i), residual );
        TS_ASSERT_DELTA( uniform->mr(i),  testStar->mr(i),  residual );
        TS_ASSERT_DELTA( uniform->P(i),   testStar->P(i),   residual );
    }

    delete testStar;
    delete uniform;
}

void test_polytrope_n1(){
    printf("STAR TESTS n=1.0 POLYTROPE\n");
    const double INDEX = 1.0;
    const std::size_t LEN = 100;
    Star *testStar = new Polytrope(INDEX,LEN);

    const double residual = ceil(testStar->SSR());

    const double Rn = sqrt( (INDEX+1.0)/(4.0*m_pi) );

    TS_ASSERT_DELTA( testStar->Radius(), Rn*m_pi, residual );
    TS_ASSERT_DELTA( testStar->Mass(), Rn*(INDEX+1)*m_pi, residual );

    double x;
    for(std::size_t i=1; i<LEN; i++){
        x = testStar->rad(i)/Rn;
        TS_ASSERT_DELTA( testStar->P(i),   pow(sin(x)/x,2), residual );
        TS_ASSERT_DELTA( testStar->rho(i), sin(x)/x, residual );
    }

    delete testStar;
}



};