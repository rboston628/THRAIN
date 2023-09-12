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
#include <cxxtest/TestSuite.h>
#include <random>

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


class StarTest : public CxxTest::TestSuite {
public:

StarTest() {
    system( "mkdir -p tests/artifacts" );
}

/***** basic tests of the SSR function *****/

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
    printf("%le %le %le\n", SSR[1], SSR[2], SSR[3]);
    double scale_factor = (SSR[1]-SSR[2])/(SSR[2]-SSR[3]);
    TS_ASSERT_LESS_THAN( scale_factor , 16.0);
    // make sure it at least scales this well
    TS_ASSERT_LESS_THAN( 6.0, scale_factor );
}

/***** tests of some basic functions using a uniform star *****/

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
    double Ac[Ncent]={0.}, Vc[Ncent]={0.}, Uc[Ncent]={0.}, cc[Ncent]={0.};
    
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
        TS_ASSERT_EQUALS( center_expand(x, Uc, maxPow), testStar->getU(i) );
        TS_ASSERT_EQUALS( center_expand(x, cc, maxPow), testStar->getC(i) );
        TS_ASSERT_DELTA( center_expand(x, Vc, maxPow), testStar->getVg(i,GAM1), 1.e-12 );
        TS_ASSERT_DELTA( center_expand(x, Ac, maxPow), testStar->getAstar(i,GAM1), 1.e-12 );
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
    // TODO why are these so much less accurate?
    for(std::size_t i=LEN-2; i>=LEN-1-NS; i--){
        t = 1. - testStar->rad(i)/R;
        TS_ASSERT_EQUALS( surface_expand_zero(t, Us, maxPow), testStar->getU(i) );
        TS_ASSERT_EQUALS( surface_expand_zero(t, cs, maxPow), testStar->getC(i) );
        TS_ASSERT_DELTA( surface_expand_inverse(t, Vs, maxPow), testStar->getVg(i,GAM1), 1.e-1 );
        TS_ASSERT_DELTA( surface_expand_inverse(t, As, maxPow), testStar->getAstar(i,GAM1), 1.e-1 );
    }

    delete testStar;
}

/***** TESTS OF POLYTROPES *****/

void do_test_center(Star* testStar, double const tol){
    /* test center expansions */
    constexpr std::size_t NC=10;
    constexpr double GAM1 = 5./3.;
    constexpr int Ncent = 3;
    int maxPow = 2*(Ncent-1); // the integer passed needs to be mutable
    double Ac[Ncent]={0.}, Vc[Ncent]={0.}, Uc[Ncent]={0.}, cc[Ncent]={0.};
    
    testStar->getAstarCenter(Ac, maxPow, GAM1);
    testStar->getVgCenter(Vc, maxPow, GAM1);
    testStar->getUCenter(Uc, maxPow);
    testStar->getC1Center(cc, maxPow);

    double x=0.0, R=testStar->Radius();
    for(std::size_t i=0; i<NC; i++){
        x = testStar->rad(i)/R;
        TS_ASSERT_DELTA( center_expand(x, Uc, maxPow), testStar->getU(i), tol );
        TS_ASSERT_DELTA( center_expand(x, cc, maxPow), testStar->getC(i), tol );
        TS_ASSERT_DELTA( center_expand(x, Vc, maxPow), testStar->getVg(i,GAM1), tol );
        TS_ASSERT_DELTA( center_expand(x, Ac, maxPow), testStar->getAstar(i,GAM1), tol );
    }
}

void do_test_surface(Star* testStar, double const tol){
    /* test surface expansions */
    constexpr std::size_t NS=10;
    constexpr double GAM1 = 5./3.;
    constexpr int Nsurf = 6;
    int maxPow = 4; // the integer passed needs to be mutable
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
        TS_ASSERT_DELTA( surface_expand_inverse(t, Vs, maxPow), testStar->getVg(i,GAM1), tol );
        TS_ASSERT_DELTA( surface_expand_inverse(t, As, maxPow), testStar->getAstar(i,GAM1), tol );
    }
}

void test_polytrope_n0(){
    /* this tests the n=0 polytrope, by comparing it to the exactly-known
    ** solution of an isopycnic (aka uniform-density) star.  The polytrope solution 
    ** should arrive exactly at the isopycnic solution, but through numerical 
    ** integration as opposed to direct setting by an equation.*/
    printf("STAR TESTS n=0.0 POLYTROPE\n");
    const double INDEX = 0.0;
    const std::size_t LEN = 1000;
    Star *uniform = new Isopycnic(LEN);
    Polytrope *testStar = new Polytrope(INDEX,LEN);

    // use the SSR of the polytrope to quantify acceptable error
    const double residual = pow(10.0,ceil(log10(testStar->SSR())));
    TS_ASSERT_LESS_THAN( residual, 1.e-8 );

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
        fprintf(fp, "%lf\t%le\n", xi, std::abs(theta-testStar->getY(i)));
    }
    fclose(fp);

    do_test_center(testStar, 1.e-12);
    do_test_surface(testStar, 1.e-4);

    delete testStar;
    delete uniform;
}

void test_polytrope_n1(){
    printf("STAR TESTS n=1.0 POLYTROPE\n");
    const double INDEX = 1.0;
    const std::size_t LEN = 1000;
    Polytrope *testStar = new Polytrope(INDEX,LEN);

    const double residual = pow(10.0,ceil(log10(testStar->SSR())));
    TS_ASSERT_LESS_THAN( residual, 1.e-8 );

    const double Rn = sqrt( (INDEX+1.0)/(4.0*m_pi) );

    TS_ASSERT_DELTA( testStar->Radius(), Rn*m_pi, residual );
    TS_ASSERT_DELTA( testStar->Mass(), Rn*(INDEX+1)*m_pi, residual );

    double x, theta;
    FILE *fp = fopen("tests/artifacts/exact_1.txt", "w");
    for(std::size_t i=1; i<LEN; i++){
        x = testStar->rad(i)/Rn;
        theta = sin(x)/x;
        TS_ASSERT_DELTA( testStar->P(i),   pow(theta,2), residual );
        TS_ASSERT_DELTA( testStar->rho(i), theta, residual );
        fprintf(fp, "%lf\t%le\n", x, std::abs(theta-testStar->getY(i)));
    }
    fclose(fp);

    do_test_center(testStar, 1.e-4);
    do_test_surface(testStar, 1.e-4);

    delete testStar;
}

void test_polytrope_n5(){
    // this is not a star, as it has no surface, but extends to infinity
    // it is useful only because it has an exact solution
    // the radius and mass cannot be tested, as they do not exist

    printf("STAR TESTS n=5.0 POLYTROPE\n");
    const double INDEX = 5.0;
    const std::size_t LEN = 1000;
    Polytrope *testStar = new Polytrope(INDEX,LEN);

    const double residual = pow(10.0,ceil(log10(testStar->SSR())));
    TS_ASSERT_LESS_THAN( residual, 1.e-8 );

    double x, theta;
    FILE *fp = fopen("tests/artifacts/exact_5.txt", "w");
    for(std::size_t i=0; i<LEN; i++){
        x = testStar->getX(i);
        theta = pow( 1.0 + x*x/3.0, -0.5);
        TS_ASSERT_DELTA( testStar->getY(i), theta, residual );
        TS_ASSERT_DELTA( testStar->P(i),   pow(theta,INDEX+1), residual );
        TS_ASSERT_DELTA( testStar->rho(i), pow(theta,INDEX), residual );
        fprintf(fp, "%lf\t%le\n", x, std::abs(theta-testStar->getY(i)));
    }
    fclose(fp);

    do_test_center(testStar, 1.e-4);
    // there is no surface

    delete testStar;
}

void test_several_polytropes(){
    printf("STAR TESTS SSR FOR SEVERAL POLYTROPES\n");
    double const INDEX[] = {1.5, 2.0, 2.5, 3.0, 3.5, 4.0};
    std::size_t const LEN = 1000;
    Star *testStar;
    for(double n : INDEX){
        testStar = new Polytrope(n, LEN);
        do_test_center(testStar, 1.e-4);
        do_test_center(testStar, 1.e-4);
        TS_ASSERT_LESS_THAN(testStar->SSR(), 1.e-8);
        delete testStar;
    }
}

void test_polytrope_MR_constructor(){
    printf("STAR TESTS POLYTROPE M,R CONSTRUCTOR\n");
    const double BigM[] = {0.2, 0.6, 1.0, 2.0, 10.0};
    const double BigR[] = {0.1, 1.1, 1.2, 20.0};
    const double index[] = {1.0, 1.5, 3.0};
    const std::size_t LEN = 1000;
    Polytrope *testStar;

    for(double M : BigM){
        for(double R : BigR){
            for(double n : index){
                testStar = new Polytrope(M*MSOLAR, R*REARTH, n, LEN);
                TS_ASSERT_DELTA(M, testStar->Mass()/MSOLAR, 1.e-14);
                TS_ASSERT_DELTA(R, testStar->Radius()/REARTH, 1.e-14);
                TS_ASSERT_LESS_THAN(testStar->SSR(), 1.e-8);
                delete testStar;
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
    printf("STAR TEST CHANDRASEKHAR TABLE\n");
    // Chandrasekhar 1939 table 27 has properties of several WDs
    // to compare must use values of A0, B0 that were in use in 1939

    std::size_t const LEN(2000);
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
        testStar = new ChandrasekharWD(Y0[n], LEN, Chandrasekhar::A01939, Chandrasekhar::B01939);
        M = testStar->Mass()/MSOLAR;
        R = testStar->Radius();
        TS_ASSERT_LESS_THAN( testStar->SSR(), 1.e-8 );
        TS_ASSERT_DELTA( M, M1939[n], 1.e-1);
        TS_ASSERT_DELTA( R, R1939[n], 1.e8);
        fprintf(fp, "%0.8lf\t%0.2lf\t%0.3lf\t%0.2le\t%0.3le\n", invY0sq[n], M1939[n], M, R1939[n], R);
        
        // test expansions
        do_test_center(testStar, 1.e-4);
        do_test_surface(testStar, 2.e-0);
        
        delete testStar;
    }
    fclose(fp);
}

void test_CHWD_grad_constructor(){
    printf("STAR TEST CHANDRASEKHAR CONSTRUCTORS\n");
    std::size_t const LEN(1000);
    ChandrasekharWD *testStar;
    const double Y0[] = {1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2.0};
    double F0;

    for(double y0 : Y0){
        testStar = new ChandrasekharWD(y0, LEN, Chandrasekhar::constant_mu{2.0});
        TS_ASSERT_LESS_THAN(testStar->SSR(), 1.e-4);
        delete testStar;
        F0 = Chandrasekhar::factor_f(sqrt(y0*y0-1.));
        testStar = new ChandrasekharWD(y0, LEN, Chandrasekhar::sigmoidal_in_logf{2.,F0,2.,1.});
        TS_ASSERT_LESS_THAN(testStar->SSR(), 1.e-4);
        delete testStar;
    }
}

// TODO read in star outputs, compare

/***** TESTS OF MESA WRAPPER ******/

// TODO read in MESA modle, previous output, compare

/***** TESTS OF SIMPLE WD MODELS ******/

// TODO create a SWD, save its output
// then create one here, compare to previous


};