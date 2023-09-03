#include "../src/STARS/Star.h"
#include "../src/STARS/Isopycnic.h"
#include "../src/STARS/Polytrope.h"
#include <cxxtest/TestSuite.h>
#include <random>

class StarTest : public CxxTest::TestSuite {
public:

void test_bad_SSR(){
    const std::size_t LEN = 10;
    Star *testStar = new Polytrope(1.0,LEN);
    // test SSR
    TS_ASSERT_LESS_THAN( 1.e-6, testStar->SSR() );
    
    delete testStar;
}

void test_uniform_star(){
    printf("STAR TESTS ISOPYCNIC\n");
    const std::size_t LEN = 100;
    Star *testStar = new Isopycnic(LEN);

    const double dx = sqrt(6.0)/double(LEN-1);
    const double Rn = sqrt( 1.0/(4.0*m_pi) );
    const double Rn3 = pow(Rn,3);

    // center
    TS_ASSERT_EQUALS(testStar->rad(0), 0.0);
    TS_ASSERT_EQUALS(testStar->P(0), 1.0);
    TS_ASSERT_EQUALS(testStar->rho(0), 1.0);
    // surface
    TS_ASSERT_EQUALS(testStar->rad(LEN-1), Rn*sqrt(6.0));
    TS_ASSERT_EQUALS(testStar->P(LEN-1), 0.0);
    TS_ASSERT_EQUALS(testStar->rho(LEN-1), 1.0);

    // test SSR
    TS_ASSERT_LESS_THAN( testStar->SSR(), 1.e-15 );

    // test all throughout star
    for(std::size_t i=0; i<LEN; i++){
        TS_ASSERT_DELTA( testStar->rad(i), Rn*double(i)*dx, 1.e-15 );
        TS_ASSERT_DELTA( testStar->P(i), 1.0 - pow(dx*i,2)/6.0, 1.e-15 );
        TS_ASSERT_DELTA( testStar->mr(i), 4.0*m_pi*pow(dx*i,3)/3.0*Rn3, 1.e-15 );
    }

    TS_ASSERT_DELTA( testStar->Radius(), Rn*sqrt(6.0), 1.e-15 );
    TS_ASSERT_DELTA( testStar->Mass(), 8.0*m_pi*sqrt(6.0)*Rn3, 1.e-15 )
    
    delete testStar;
}

void test_polytrope_n0(){
    printf("STAR TESTS n=0.0 POLYTROPE\n");
    const double INDEX = 0.0;
    const std::size_t LEN = 100;
    Star *uniform = new Isopycnic(LEN);
    Star *testStar = new Polytrope(INDEX,LEN);

    const double residual = ceil(testStar->SSR());

    TS_ASSERT_DELTA( uniform->Radius(), testStar->Radius(), residual );
    TS_ASSERT_DELTA( uniform->Mass(), testStar->Mass(), residual );

    for(std::size_t i=0; i<LEN; i++){
        TS_ASSERT_DELTA( uniform->rad(i), testStar->rad(i), residual );
        TS_ASSERT_DELTA( uniform->mr(i),  testStar->mr(i),  residual );
        TS_ASSERT_DELTA( uniform->P(i),   testStar->P(i),   residual );
        TS_ASSERT_DELTA( uniform->rho(i), testStar->rho(i), residual );
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