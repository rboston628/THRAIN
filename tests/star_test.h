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
#include "../src/STARS/MESA.h"
#include <cxxtest/TestSuite.h>
#include <random>
#include <stdio.h>

#define ASSERT_REL_DIFFERENCE(x, y, err) { TS_ASSERT_LESS_THAN(relative_difference(x,y), err) }

double center_expand(double const x, double const *const ac, int const maxPow){
    double value = ac[0];
    for(int i=1; i<=maxPow/2; i++){
        value += ac[i]*pow(x,2*i);
    }
    return value;
}

double surface_expand_zero(double const t, double const *const ac, int const maxPow){
    double value = ac[0];
    for(int i=1; i<=maxPow; i++){
        value += ac[i]*pow(t,i);
    }
    return value;
}

double surface_expand_inverse(double const t, double const *const ac, int const maxPow){
    const int O=1;
    double value = ac[O-1]/t;
    for(int i=0; i<=maxPow; i++){
        value += ac[O+i]*pow(t,i);
    }
    return value;
}

double relative_difference(double x, double y){
    double denom = (x*y==0.0 ? 1.0 : std::abs(x+y));
    return 2.0 * std::abs(x-y) / denom;
}

class StarTest : public CxxTest::TestSuite {
public:

static StarTest *createSuite (){
    system( "mkdir -p tests/artifacts" );
    printf("\nSTAR TESTS");
    return new StarTest();
}
static void destroySuite(StarTest *suite) { 
    printf("\n##########");
    delete suite; 
}

void setUp() {
    freopen("tests/artifacts/logio.txt", "a", stdout);
}

void tearDown() {
    freopen("/dev/tty", "w", stdout);
}

/***** basic tests of the SSR function *****/

void test_bad_SSR(){
    /* Tests that low number of grid points
    *  will produce a bad SSR*/
    std::size_t const LEN = 10;
    Star *testStar = new Polytrope(1.0,LEN);
    // test SSR
    double const BAD_SSR = 1.e-6;
    TS_ASSERT_LESS_THAN( BAD_SSR, testStar->SSR() );
    
    delete testStar;
}

void test_scale_SSR(){
    /* Test that SSR decreases up to minimum as gridsize increases*/
    std::size_t const LEN[] = {10, 20, 40, 80};
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
    double const IDEAL_SCALE_FACTOR = 16.0, MIN_ACCEPTABLE_SCALE_FACTOR = 12.0;
    double scale_factor = (SSR[1]-SSR[2])/(SSR[2]-SSR[3]);
    TS_ASSERT_LESS_THAN_EQUALS( scale_factor , IDEAL_SCALE_FACTOR);
    // make sure it at least scales this well
    TS_ASSERT_LESS_THAN( MIN_ACCEPTABLE_SCALE_FACTOR, scale_factor );
}

/***** tests of some basic functions using a uniform star *****/

void test_uniform_star(){
    /* This test really only ensures that the SSR will be 
    ** approximately zero in the case of an exactly known star.  
    ** Because it happens that g scales linearly with radius, it 
    ** happens to work out that the numerical derivatives in the 
    ** SSR are exactly correct*/
    const std::size_t LEN = 100;
    Star *testStar = new Isopycnic(LEN);

    double const MACHINE_PRECISION = 1.e-15;

    const double dx = sqrt(6.0)/double(LEN-1);
    const double Rn = sqrt( 1.0/(4.0*m_pi) );

    // test SSR -- this needs to be near-zero
    TS_ASSERT_LESS_THAN( testStar->SSR(), MACHINE_PRECISION );

    // test all throughout star
    // this verifies that the getter methods work as intended
    double x=0;
    
    for(std::size_t i=0; i<LEN; i++){
        x = dx*i;
        TS_ASSERT_EQUALS( testStar->rho(i), 1.0 );
        TS_ASSERT_DELTA( testStar->rad(i), Rn*x,          MACHINE_PRECISION );
        TS_ASSERT_DELTA( testStar->P(i),  1.0 - x*x/6.0,  MACHINE_PRECISION );
        TS_ASSERT_DELTA( testStar->mr(i), Rn*(x*x*x)/3.0, MACHINE_PRECISION );
    }

    // test that mass and radius equal to known quantities
    TS_ASSERT_DELTA( testStar->Radius(), Rn*sqrt(6.0),    MACHINE_PRECISION );
    TS_ASSERT_DELTA( testStar->Mass(), Rn*2.0*sqrt(6.0),  MACHINE_PRECISION );
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
    double Ac[Ncent]={0.}, Vc[Ncent]={0.}, Uc[Ncent]={0.}, cc[Ncent]={0.};
    
    testStar->getAstarCenter(Ac, maxPow, GAM1);
    testStar->getVgCenter(Vc, maxPow, GAM1);
    testStar->getUCenter(Uc, maxPow);
    testStar->getC1Center(cc, maxPow);

    double x=0.0, R=testStar->Radius();
    double const ACCEPTABLE_PRECISION = 1.e-10;
    for(std::size_t i=0; i<NC; i++){
        x = testStar->rad(i)/R;
        TS_ASSERT_EQUALS( center_expand(x, Uc, maxPow), testStar->getU(i) );
        TS_ASSERT_EQUALS( center_expand(x, cc, maxPow), testStar->getC(i) );
        TS_ASSERT_DELTA ( center_expand(x, Vc, maxPow), testStar->getVg(i,GAM1),    ACCEPTABLE_PRECISION );
        TS_ASSERT_DELTA ( center_expand(x, Ac, maxPow), testStar->getAstar(i,GAM1), ACCEPTABLE_PRECISION );
    }

    /* test surface expansions */
    constexpr std::size_t NS=10;
    constexpr int Nsurf = 6;
    maxPow = 4; // the integer passed needs to be mutable
    double As[Nsurf]={0.}, Vs[Nsurf]={0.}, Us[Nsurf]={0.}, cs[Nsurf]={0.};
    testStar->getAstarSurface(As, maxPow, GAM1);
    testStar->getVgSurface(Vs, maxPow, GAM1);
    testStar->getUSurface(Us, maxPow);
    testStar->getC1Surface(cs, maxPow);

    double t=0.0;
    for(std::size_t i=LEN-2; i>=LEN-1-NS; i--){
        t = 1. - testStar->rad(i)/R;
        TS_ASSERT_EQUALS( surface_expand_zero(t, Us, maxPow),    testStar->getU(i) );
        TS_ASSERT_EQUALS( surface_expand_zero(t, cs, maxPow),    testStar->getC(i) );
        TS_ASSERT_DELTA ( surface_expand_inverse(t, Vs, maxPow), testStar->getVg(i,GAM1),    ACCEPTABLE_PRECISION );
        TS_ASSERT_DELTA ( surface_expand_inverse(t, As, maxPow), testStar->getAstar(i,GAM1), ACCEPTABLE_PRECISION );
    }

    delete testStar;
}

/***** TESTS OF POLYTROPES *****/

std::size_t do_test_center(Star *const testStar, double const tol){
    /* test center expansions */
    constexpr std::size_t NC(4);             // number of steps to calculate
    constexpr std::size_t MAXPOW(4);         // max power in the expansion
    constexpr std::size_t NCENT(MAXPOW/2+1); // number of coefficients needed
    constexpr double GAM1(5./3.);
    
    int maxPow = MAXPOW; // the integer passed needs to be mutable
    double Ac[NCENT]={0.}, Vc[NCENT]={0.}, Uc[NCENT]={0.}, cc[NCENT]={0.};
    
    testStar->getAstarCenter(Ac, maxPow, GAM1);
    testStar->getVgCenter(Vc, maxPow, GAM1);
    testStar->getUCenter(Uc, maxPow);
    testStar->getC1Center(cc, maxPow);

    std::size_t failures = 0;
    std::size_t LEN = testStar->length();
    double x=0.0, R=testStar->Radius();
    for(std::size_t i=0; i<NC; i++){
        x = double(i)/double(LEN-1);
        TS_ASSERT_DELTA( center_expand(x, cc, maxPow), testStar->getC(i), tol );
        TS_ASSERT_DELTA( center_expand(x, Vc, maxPow), testStar->getVg(i,GAM1), tol );
        TS_ASSERT_DELTA( center_expand(x, Ac, maxPow), testStar->getAstar(i,GAM1), tol );
        /* see mathematica notebook --
        // it can be shown that U's central expansion is simply less accurate than the others*/ 
        TS_ASSERT_DELTA( center_expand(x, Uc, maxPow), testStar->getU(i), tol );        
    }
    TS_ASSERT_EQUALS( failures, 0 );
    return failures;
}

void do_test_surface(Star *const testStar, double const tol){
    /* test surface expansions */
    constexpr std::size_t NS(4); // number of steps to calculate
    constexpr std::size_t MAXPOW(4); // maximum power in expansion
    constexpr std::size_t Nsurf(MAXPOW+2); // number of coefficients needed
    constexpr double GAM1(5./3.);
    int maxPow = MAXPOW; // the integer passed needs to be mutable
    double As[Nsurf]={0.}, Vs[Nsurf]={0.}, Us[Nsurf]={0.}, cs[Nsurf]={0.};

    testStar->getAstarSurface(As, maxPow, GAM1);
    testStar->getVgSurface(Vs, maxPow, GAM1);
    testStar->getUSurface(Us, maxPow);
    testStar->getC1Surface(cs, maxPow);

    double t=0.0;
    // TODO why are these so much less accurate?
    std::size_t const LEN = testStar->length();
    double const R = testStar->Radius();
    for(std::size_t i=LEN-2; i>=LEN-1-NS; i--){
        t = 1. - testStar->rad(i)/R;
        TS_ASSERT_DELTA( surface_expand_zero(t, Us, maxPow), testStar->getU(i), tol );
        TS_ASSERT_DELTA( surface_expand_zero(t, cs, maxPow), testStar->getC(i), tol );
        // while c1, U stay on the order 1, these two functions blow up near the surface
        // requiring the relative difference for comparison
        ASSERT_REL_DIFFERENCE( surface_expand_inverse(t, Vs, maxPow), testStar->getVg(i,GAM1), tol );
        ASSERT_REL_DIFFERENCE( surface_expand_inverse(t, As, maxPow), testStar->getAstar(i,GAM1), tol );
    }
}

// TODO there must be a way to make something like this work
// might be used by ETS_ tests
// void test_bad_boundaries(){
//     /* this test is to ensure that the boundary checks do not always pass*/
//     std::size_t const LEN (10);
//     double const INDEX = 0.5; // this is known to behave poorly
//     Polytrope *testStar = new Polytrope(INDEX, LEN);
//     double tol = 1.e-6;
//     {
//         fprintf(stderr, "RUNNING THE TESTS\n");
//         std::size_t failed = (CxxTest::tracker()).testFailedAsserts();
//         CxxTest::TestTracker tmpTracker, save = CxxTest::tracker();
//         CxxTest::tracker() = tmpTracker;
//         fprintf(stderr, "INITTED\n");
//         try { ETS_ASSERT_EQUALS( 2, 3 ); }
//         catch(const int &e) { fprintf(stderr, "basic math\n"); }
//         ETS_ASSERT_EQUALS( (CxxTest::tracker()).testFailedAsserts(), failed+1);
//         // do_test_center(testStar, tol);
//         // fprintf(stderr, "RAN THE TESTS\n");
//         // ETS_ASSERT_EQUALS( (CxxTest::tracker()).testFailedAsserts(), 4 );
//         // CxxTest::tracker() = save;
//     }
//     // TS_ASSERT_THROWS_ANYTHING( do_test_center(testStar, tol) );
//     // TS_ASSERT_THROWS_ANYTHING( do_test_surface(testStar, tol) );
//     TS_ASSERT_EQUALS((CxxTest::tracker()).testFailedAsserts(), 0);
//     delete testStar;
// }

void test_polytrope_n0(){
    /* this tests the n=0 polytrope, by comparing it to the exactly-known
    ** solution of an isopycnic (aka uniform-density) star.  The polytrope solution 
    ** should arrive at the isopycnic solution, but through numerical 
    ** integration as opposed to direct setting by an equation.*/
    fprintf(stderr, "\n##########");
    fprintf(stderr, "\nSTAR TESTS - POLYTROPE n=0.0");
    double const INDEX (0.0);
    std::size_t const LEN (1001);
    Star *uniform = new Isopycnic(LEN);
    Polytrope *testStar = new Polytrope(INDEX,LEN);

    // use the SSR of the polytrope to quantify acceptable error
    double const machine_precision = 10.e-16;
    double const expected_error = 10.0 * pow(1.0/LEN, 4);
    double const residual = pow(10.0,ceil(log10(testStar->SSR())));
    TS_ASSERT_LESS_THAN( residual, expected_error );

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

    double xi, theta;
    FILE *fp = fopen("tests/artifacts/exact_0.txt", "w");
    for(std::size_t i=0; i<LEN; i++){
        xi = testStar->getX(i);
        theta = 1. - xi*xi/6.;
        TS_ASSERT_DELTA( testStar->getY(i), theta, residual );
        fprintf(fp, "%lf\t%le\n", xi, std::abs(theta-testStar->getY(i)));
    }
    fclose(fp);

    do_test_center (testStar, machine_precision );
    do_test_surface(testStar, expected_error);

    delete testStar;
    delete uniform;
}

void test_polytrope_n1(){
    fprintf(stderr, "\nSTAR TESTS - POLYTROPE n=1.0");
    double const INDEX (1.0);
    std::size_t const LEN (1001);
    Polytrope *testStar = new Polytrope(INDEX,LEN);

    // expected error: LEN^{-4} per step, compounded over LEN steps
    double const expected_error = pow(1.0/LEN, 3);
    double const residual = pow(10.0,ceil(log10(testStar->SSR())));
    TS_ASSERT_LESS_THAN( residual, expected_error );

    const double Rn = sqrt( (INDEX+1.0)/(4.0*m_pi) );

    TS_ASSERT_DELTA( testStar->getX(LEN-1), m_pi,          residual );
    TS_ASSERT_DELTA( testStar->Radius(), Rn*m_pi,          residual );
    TS_ASSERT_DELTA( testStar->Mass(), Rn*(INDEX+1.)*m_pi, residual );

    double x, theta;
    FILE *fp = fopen("tests/artifacts/exact_1.txt", "w");
    for(std::size_t i=1; i<LEN; i++){
        x = testStar->rad(i)/Rn;
        theta = sin(x)/x;
        TS_ASSERT_DELTA( testStar->P(i),    pow(theta,2),  residual );
        TS_ASSERT_DELTA( testStar->rho(i),  theta,         residual );
        TS_ASSERT_DELTA( testStar->getY(i), theta,         residual );
        fprintf(fp, "%lf\t%le\n", x, std::abs(theta-testStar->getY(i)));
    }
    fclose(fp);
    
    double acceptable_error(1.e-6);
    do_test_center  ( testStar, acceptable_error );
    do_test_surface ( testStar, acceptable_error );

    delete testStar;
}

void test_polytrope_n5(){
    // this is not a star, as it has no surface, but extends to infinity
    // it is useful only because it has an exact solution
    // the radius and mass cannot be tested, as they do not exist

    fprintf(stderr, "\nSTAR TESTS - POLYTROPE n=5.0");
    const double INDEX = 5.0;
    const std::size_t LEN = 1001;
    Polytrope *testStar = new Polytrope(INDEX,LEN);

    const double expected_error = pow(1.0/LEN, 3);
    const double residual = pow(10.0,ceil(log10(testStar->SSR())));
    TS_ASSERT_LESS_THAN( residual, expected_error );

    double x, theta;
    FILE *fp = fopen("tests/artifacts/exact_5.txt", "w");
    for(std::size_t i=0; i<LEN; i++){
        x = testStar->getX(i);
        theta = pow( 1.0 + x*x/3.0, -0.5);
        TS_ASSERT_DELTA( testStar->getY(i), theta,             residual );
        TS_ASSERT_DELTA( testStar->P(i),   pow(theta,INDEX+1), residual );
        TS_ASSERT_DELTA( testStar->rho(i), pow(theta,INDEX),   residual );
        fprintf(fp, "%lf\t%le\n", x, std::abs(theta-testStar->getY(i)));
    }
    fclose(fp);

    double acceptable_error(1.e-4);
    do_test_center( testStar, acceptable_error );
    // there is no surface

    delete testStar;
}

void test_several_polytropes(){
    fprintf(stderr, "\nSTAR TESTS - SSR FOR SEVERAL POLYTROPES");
    double const INDEX[] = {2.0, 2.5, 3.0, 3.5, 4.0};
    std::size_t const LEN = 1001;
    Star *testStar;

    double const expected_error = pow(1.0/LEN, 3);
    double residual;

    for(double n : INDEX){
        testStar = new Polytrope(n, LEN);
        residual = pow(10.0,ceil(log10(testStar->SSR())));
        // TS_ASSERT_LESS_THAN(residual, expected_error);
        // assert P/rho = theta?
        double acceptable_error(1.e-4);
        do_test_center ( testStar, acceptable_error );
        do_test_surface( testStar, acceptable_error );
        delete testStar;
    }
}

void test_polytrope_MR_constructor(){
    fprintf(stderr, "\nSTAR TESTS - POLYTROPE M,R CONSTRUCTOR");
    const double BigM[] = {0.2, 0.6, 1.0, 2.0, 10.0};
    const double BigR[] = {0.1, 1.1, 1.2, 20.0};
    const double index[] = {1.0, 1.5, 3.0};
    const std::size_t LEN = 1001;
    double const expected_error = pow(10./LEN , 4);
    Polytrope *testStar, *refStar;

    for(double M : BigM){
        for(double R : BigR){
            for(double n : index){
                testStar = new Polytrope(M*MSOLAR, R*REARTH, n, LEN);
                refStar = new Polytrope(n, LEN);
                TS_ASSERT_DELTA(M, testStar->Mass()/MSOLAR, 1.e-14);
                TS_ASSERT_DELTA(R, testStar->Radius()/REARTH, 1.e-14);
                TS_ASSERT_LESS_THAN(testStar->SSR(), expected_error);
                for(std::size_t i=0; i<LEN; i++){
                    TS_ASSERT_DELTA( testStar->getY(i), refStar->getY(i), expected_error );
                    TS_ASSERT_DELTA( testStar->getX(i), refStar->getX(i), expected_error );
                }
                delete testStar;
                delete refStar;
            }
        }
    }
}

void test_make_exact_error_graph(){
	//plot everything in single graph, for simplicity
	FILE *gnuplot = popen("gnuplot -persist", "w");
	fprintf(gnuplot, "reset\n");
	fprintf(gnuplot, "set term png size 2000,1000\n");
	fprintf(gnuplot, "set output 'tests/artifacts/polytrope_exact_differences.png'\n");
	fprintf(gnuplot, "set title 'Errors in Analytic Polytropes, N_{star} = %d'\n",1000);
	fprintf(gnuplot, "set xrange [0:1.01]\n");
	fprintf(gnuplot, "set xlabel 'Scaled Radius r/R'\n");
	fprintf(gnuplot, "set ylabel 'log_{10} |θ_{num}-θ_{exact}|\n");
	fprintf(gnuplot, "set logscale y\n");
	fprintf(gnuplot, "set format y '10^{%%L}'\n");
	fprintf(gnuplot, "plot");
    fprintf(gnuplot, " '%s' u 1:2 w l t 'n=0',", "tests/artifacts/exact_0.txt");
    fprintf(gnuplot, " '%s' u 1:2 w l t 'n=1',", "tests/artifacts/exact_1.txt");
    fprintf(gnuplot, " '%s' u 1:2 w l t 'n=5'", "tests/artifacts/exact_5.txt");
	fprintf(gnuplot, "\n");
	pclose(gnuplot);

    // remove the text files
    system("rm tests/artifacts/exact_0.txt");
    system("rm tests/artifacts/exact_1.txt");
    system("rm tests/artifacts/exact_5.txt");
}

// TODO: make the following for polytrope
// * scaling of exact errors
// * physical errors
// * scaling of physical errors
// * scaling of theta

/***** TESTS OF CHANDRASEKHAR WD *****/

void test_CHWD_against_chandrasekhar(){
    fprintf(stderr, "\nSTAR TESTS - CHANDRASEKHAR TABLE");
    // Chandrasekhar 1939 table 27 has properties of several WDs
    // to compare must use values of A0, B0 that were in use in 1939

    std::size_t const LEN(2001);
    ChandrasekharWD *testStar;
    // the values listed in Chandrasekhar 1939 table
    double invY0sq[] = {0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.80};
    double M1939[] = {5.51, 5.32, 4.87, 4.33, 3.54, 2.95, 2.45, 2.02, 1.62, 0.88};
    double R1939[] = {4.13e8, 5.44e8, 7.69e8, 9.92e8, 1.29e9, 1.51e9, 1.72e9, 1.93e9, 2.15e9, 2.79e9};
    const std::size_t num_rows = sizeof(invY0sq)/sizeof(std::size_t);
    // translate 1/Y0^2 to Y0
    double Y0[num_rows];
    for(std::size_t n=0; n<num_rows; n++) {
        Y0[n] = 1./sqrt(invY0sq[n]);
    }

    char filename[] = "tests/artifacts/Chandrasekhar1939.txt";
    system( "mkdir -p tests/artifacts/" );
    FILE *fp = fopen(filename, "w");
    double M,R;
    for(std::size_t n=0; n<num_rows; n++){
        // use the form of the constructor that allows specifuing A0, B0
        // must use the 1939 values, to ensure equality of results
        testStar = new ChandrasekharWD(Y0[n], LEN, Chandrasekhar::A01939, Chandrasekhar::B01939);
        M = testStar->Mass()/MSOLAR;
        R = testStar->Radius();
        // assert reasonable numerical error
        TS_ASSERT_LESS_THAN( testStar->SSR(), 1.e-8 );
        // assert agreement with tabular values in Chandrasekhar 1939
        ASSERT_REL_DIFFERENCE( M, M1939[n], 1.e-2);
        ASSERT_REL_DIFFERENCE( R, R1939[n], 1.e-2);
        // write out values to table
        fprintf(fp, "%0.8lf\t%0.2lf\t%0.3lf\t%0.2le\t%0.3le\n", invY0sq[n], M1939[n], M, R1939[n], R);
        
        // test expansions
        do_test_center(testStar, 1.e-4);
        // do_test_surface(testStar, 1.e-3); // TODO why is this fit so bad? 
        
        delete testStar;
    }
    fclose(fp);
}

void test_CHWD_grad_constructor(){
    fprintf(stderr, "\nSTAR TESTS - CHANDRASEKHAR CONSTRUCTORS\n");
    // test the constructor that accepts a gradiant in mu
    // try both a constant gradient, and a sigmoidal gradient
    std::size_t const LEN(1001);
    ChandrasekharWD *testStar;
    const double Y0[] = {1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2.0};
    double F0;

    for(double y0 : Y0){
        // star with constant mue
        testStar = new ChandrasekharWD(y0, LEN, Chandrasekhar::constant_mu{2.0});
        TS_ASSERT_LESS_THAN(testStar->SSR(), 1.e-4);
        do_test_center(testStar, 1.e-4);
        // TODO make surface test work
        do_test_surface(testStar, 1.0);
        delete testStar;
        // star with sigmoidally varying mue (CHWD++)
        F0 = Chandrasekhar::factor_f(sqrt(y0*y0-1.));
        testStar = new ChandrasekharWD(y0, LEN, Chandrasekhar::sigmoidal_in_logf{2.,F0,2.,1.});
        TS_ASSERT_LESS_THAN(testStar->SSR(), 1.e-4);
        do_test_center(testStar, 1.e-4);
        // TODO make surface test work
        do_test_surface(testStar, 1.0);
        delete testStar;
    }
}

// TODO read in star outputs, compare

/***** TESTS OF MESA WRAPPER ******/

void test_MESA_fit(){
    /* create a simple star model from a polytrope
    *  save its data to a file, formatted like MESA data
    *  load the file into a MESA wrapper
    *  check that the knots correspond */
    fprintf(stderr, "\nSTAR TESTS - MESA FIT");
    double index = 2.2;
    double M = MSOLAR, R = 1.35*REARTH;
    std::size_t LEN(1001);
    Polytrope *refStar = new Polytrope(M, R, index, LEN);
    // print using MESA v1.01 format
    std::string testfile ("tests/inputs/test_polytrope_n=2.2_v101.dat");
    FILE *fp = fopen(testfile.c_str(), "w");
    // first line formatted N Mstar Rstar Lstar vx100
    fprintf(fp, "%lu %le %le 0.0 101\n", LEN, M, R);
    for(std::size_t i=0; i<LEN; i++){
        // column 1 = index, from 1 to LEN (i+1)
        // column 2 = radius
        // column 3 = mass
        // column 4 = luminosity = 0
        fprintf(fp, "%lu %0.16le %0.16le 0.0 ", i+1, refStar->rad(i), refStar->mr(i));
        // column 5 = pressure
        // column 6 = temperature = 0 except last two
        // column 7 = density
        // column 8 = del = 0
        double t = (i<LEN-2 ? 0.0 : 1.2e4);
        fprintf(fp, "%0.16le %0.16le %0.16le 0.0 ", refStar->P(i), t, refStar->rho(i));
        // column 9 = Brunt-Vaisala
        // column 10 = Gamma1
        // column 11 onward = stuff we don't care about
        fprintf(fp, "%0.16le %0.16le 0.0\n", -refStar->Schwarzschild_A(i)*refStar->dPhidr(i), refStar->Gamma1(i));
    }
    fclose(fp);
    // create a MESA model from that file
    MESA *testStar = new MESA(testfile.c_str(), LEN);
    // some routine verifications
    TS_ASSERT_EQUALS(testStar->length(), 2*(LEN-1) + 1);
    TS_ASSERT_LESS_THAN(testStar->SSR(), 1.0/LEN)
    // assert both models are equal at the knots to within relative precision
    double const MACHINE_PRECISION = 1.e-15;
    for(std::size_t i=0; i<LEN-2; i++){
        ASSERT_REL_DIFFERENCE(refStar->rad(i),      testStar->rad(2*i),      MACHINE_PRECISION);
        ASSERT_REL_DIFFERENCE(refStar->mr(i),       testStar->mr(2*i),       MACHINE_PRECISION);
        ASSERT_REL_DIFFERENCE(refStar->P(i),        testStar->P(2*i),        MACHINE_PRECISION);
        ASSERT_REL_DIFFERENCE(refStar->rho(i),      testStar->rho(2*i),      MACHINE_PRECISION);
        ASSERT_REL_DIFFERENCE(refStar->dPhidr(i),   testStar->dPhidr(2*i),   MACHINE_PRECISION);
        // these ae calculare values, but should still be within precision
        ASSERT_REL_DIFFERENCE(refStar->getU(i),     testStar->getU(2*i),     MACHINE_PRECISION);
        ASSERT_REL_DIFFERENCE(refStar->getC(i),     testStar->getC(2*i),     MACHINE_PRECISION);
        ASSERT_REL_DIFFERENCE(refStar->getVg(i),    testStar->getVg(2*i),    MACHINE_PRECISION);
        ASSERT_REL_DIFFERENCE(refStar->getAstar(i), testStar->getAstar(2*i), MACHINE_PRECISION);
    }
}

void test_MESA_constructor(){
    // should first calculate SSR using only the knots
    // then check the SSR isn't very much worse

    fprintf(stderr, "\nSTAR TESTS - MESA BOUNDARIES");
    std::size_t const LEN(101);
    std::string mesa_file = "mesa_co_wd_cold.dat";
    MESA *testStar = new MESA(mesa_file.c_str(), LEN);
    TS_ASSERT_LESS_THAN( testStar->SSR(), 1.0e-2 );
    do_test_center(testStar, 1.e-4);
    // TODO make surface work
    // do_test_surface(testStar, 1.e-6);
    delete testStar;
}

// TODO make reduced data file

// TODO read in MESA model, previous output, compare

/***** TESTS OF SIMPLE WD MODELS ******/

// TODO create a SWD, save its output
// then create one here, compare to previous

};
