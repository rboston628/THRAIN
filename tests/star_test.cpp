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
#include "../src/STARS/MESA.h"
#include "../src/STARS/SimpleWD.h"
#include "doctest.h"
//import for testing
#include "test_stars/Isopycnic.h"
#include <random>
#include <stdio.h>
#include <memory>

namespace {
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

std::size_t do_test_center(Star *const testStar, double const tol){
  /* test center expansions */
  constexpr std::size_t NC(4);       // number of steps to calculate
  constexpr std::size_t MAXPOW(4);     // max power in the expansion
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
    CHECK_DELTA( center_expand(x, cc, maxPow), testStar->getC(i), tol );
    CHECK_DELTA( center_expand(x, Vc, maxPow), testStar->getVg(i,GAM1), tol );
    CHECK_DELTA( center_expand(x, Ac, maxPow), testStar->getAstar(i,GAM1), tol );
    /* see mathematica notebook --
    // it can be shown that U's central expansion is simply less accurate than the others*/ 
    CHECK_DELTA( center_expand(x, Uc, maxPow), testStar->getU(i), tol );    
  }
  CHECK_EQ( failures, 0 );
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
    CHECK_DELTA( surface_expand_zero(t, Us, maxPow), testStar->getU(i), tol );
    CHECK_DELTA( surface_expand_zero(t, cs, maxPow), testStar->getC(i), tol );
    // while c1, U stay on the order 1, these two functions blow up near the surface
    // requiring the relative difference for comparison
    CHECK_EPS( surface_expand_inverse(t, Vs, maxPow), testStar->getVg(i,GAM1), tol );
    CHECK_EPS( surface_expand_inverse(t, As, maxPow), testStar->getAstar(i,GAM1), tol );
  }
}
} // namespace


TEST_SUITE("Star [unit]") {

/***** basic tests of the SSR function *****/

TEST_CASE("star: star_bad_SSR") {
  /* Tests that low number of grid points
  *  will produce a bad SSR*/
  std::size_t const LEN = 10;
  auto testStar = std::make_unique<Polytrope>(1.0, LEN);
  // test SSR
  double const BAD_SSR = 1.e-6;
  CHECK_LT( BAD_SSR, testStar->SSR() );
}

TEST_CASE("star: star_scale_SSR") {
  /* Test that SSR decreases up to minimum as gridsize increases*/
  std::size_t const LEN[] = {10, 20, 40, 80};
  const int Ntest = sizeof(LEN)/sizeof(std::size_t);
  double SSR[Ntest];
  for(int i=0; i<4; i++){
    auto testStar = std::make_unique<Polytrope>(1.0, LEN[i]);
    SSR[i] = testStar->SSR();
  }
  
  // test SSRs grow smaller with grid size
  for(int i=1; i<Ntest; i++){
    CHECK_LT(SSR[i], fmax(SSR[i-1], 1.e-14));
  }
  
  // further test that it SCALES like LEN^{-4}
  // if it scaled exactly like LEN^{-4}, this would equal 2^4=16
  // because it reaches limits of machine precision, will be smaller
  double const IDEAL_SCALE_FACTOR = 16.0, MIN_ACCEPTABLE_SCALE_FACTOR = 12.0;
  double scale_factor = (SSR[1]-SSR[2])/(SSR[2]-SSR[3]);
  CHECK_LE( scale_factor , IDEAL_SCALE_FACTOR);
  // make sure it at least scales this well
  CHECK_LT( MIN_ACCEPTABLE_SCALE_FACTOR, scale_factor );
}

/***** tests of some basic functions using a uniform star *****/

TEST_CASE("star: star_uniform_star") {
  /* This test really only ensures that the SSR will be 
  ** approximately zero in the case of an exactly known star.  
  ** Because it happens that g scales linearly with radius, it 
  ** happens to work out that the numerical derivatives in the 
  ** SSR are exactly correct*/
  const std::size_t LEN = 100;
  auto testStar = std::make_unique<Isopycnic>(LEN);

  double const MACHINE_PRECISION = 1.e-15;

  const double dx = sqrt(6.0)/double(LEN-1);
  const double Rn = sqrt( 1.0/(4.0*m_pi) );

  // test SSR -- this needs to be near-zero
  CHECK_LT( testStar->SSR(), MACHINE_PRECISION );

  // test all throughout star
  // this verifies that the getter methods work as intended  
  for(std::size_t i=0; i<LEN; i++){
    double const x = dx*i;
    CHECK_EQ( testStar->rho(i), 1.0 );
    CHECK_DELTA( testStar->rad(i), Rn*x,      MACHINE_PRECISION );
    CHECK_DELTA( testStar->P(i),  1.0 - x*x/6.0,  MACHINE_PRECISION );
    CHECK_DELTA( testStar->mr(i), Rn*(x*x*x)/3.0, MACHINE_PRECISION );
  }

  // test that mass and radius equal to known quantities
  CHECK_DELTA( testStar->Radius(), Rn*sqrt(6.0),  MACHINE_PRECISION );
  CHECK_DELTA( testStar->Mass(), Rn*2.0*sqrt(6.0),  MACHINE_PRECISION );
}

TEST_CASE("star: star_uniform_boundaries") {
  /* this will test that the center and surface expansions
  ** in this known case act as good approxmations to the functions*/

  constexpr std::size_t LEN = 1000;
  constexpr double GAM1 = 5./3.;
  auto testStar = std::make_unique<Isopycnic>(LEN);

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
    CHECK_EQ(    center_expand(x, Uc, maxPow), testStar->getU(i) );
    CHECK_EQ(    center_expand(x, cc, maxPow), testStar->getC(i) );
    CHECK_DELTA( center_expand(x, Vc, maxPow), testStar->getVg(i,GAM1),  ACCEPTABLE_PRECISION );
    CHECK_DELTA( center_expand(x, Ac, maxPow), testStar->getAstar(i,GAM1), ACCEPTABLE_PRECISION );
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
    CHECK_EQ(    surface_expand_zero(   t, Us, maxPow),  testStar->getU(i) );
    CHECK_EQ(    surface_expand_zero(   t, cs, maxPow),  testStar->getC(i) );
    CHECK_DELTA( surface_expand_inverse(t, Vs, maxPow), testStar->getVg(i,GAM1),  ACCEPTABLE_PRECISION );
    CHECK_DELTA( surface_expand_inverse(t, As, maxPow), testStar->getAstar(i,GAM1), ACCEPTABLE_PRECISION );
  }
}

/***** TESTS OF POLYTROPES *****/

// // TODO enable this test
// // by better handling errors in Polytrope.cpp
// void test_fail_polytrope_bad_inputs(){
//   std::size_t const LEN(10);
//   double const low_index(-1.0e-3);
//   double const high_index(5.001);
//   TS_ASSERT_THROWS_ANYTHING(Polytrope(low_index, LEN));
//   TS_ASSERT_THROWS_ANYTHING(Polytrope(low_index, LEN));
// }

// TODO there must be a way to make something like this work
// might be used by ETS_ tests
// void test_bad_boundaries(){
//   /* this test is to ensure that the boundary checks do not always pass*/
//   std::size_t const LEN (10);
//   double const INDEX = 0.5; // this is known to behave poorly
//   Polytrope *testStar = new Polytrope(INDEX, LEN);
//   double tol = 1.e-6;
//   {
//     fprintf(stderr, "RUNNING THE TESTS\n");
//     std::size_t failed = (CxxTest::tracker()).testFailedAsserts();
//     CxxTest::TestTracker tmpTracker, save = CxxTest::tracker();
//     CxxTest::tracker() = tmpTracker;
//     fprintf(stderr, "INITTED\n");
//     try { ECHECK_EQ( 2, 3 ); }
//     catch(const int &e) { fprintf(stderr, "basic math\n"); }
//     ECHECK_EQ( (CxxTest::tracker()).testFailedAsserts(), failed+1);
//     // do_test_center(testStar, tol);
//     // fprintf(stderr, "RAN THE TESTS\n");
//     // ECHECK_EQ( (CxxTest::tracker()).testFailedAsserts(), 4 );
//     // CxxTest::tracker() = save;
//   }
//   // TS_ASSERT_THROWS_ANYTHING( do_test_center(testStar, tol) );
//   // TS_ASSERT_THROWS_ANYTHING( do_test_surface(testStar, tol) );
//   CHECK_EQ((CxxTest::tracker()).testFailedAsserts(), 0);
//   delete testStar;
// }

TEST_CASE("star: star_polytrope_n0") {
  /* this tests the n=0 polytrope, by comparing it to the exactly-known
  ** solution of an isopycnic (aka uniform-density) star.  The polytrope solution 
  ** should arrive at the isopycnic solution, but through numerical 
  ** integration as opposed to direct setting by an equation.*/
  double const INDEX (0.0);
  std::size_t const LEN (1001);
  auto uniform = std::make_unique<Isopycnic>(LEN);
  auto testStar = std::make_unique<Polytrope>(INDEX, LEN);

  // use the SSR of the polytrope to quantify acceptable error
  double const machine_precision = 10.e-16;
  double const expected_error = 10.0 * pow(1.0/LEN, 4);
  double const residual = pow(10.0,ceil(log10(testStar->SSR())));
  CHECK_LT( residual, expected_error );

  // test they have same mass, radius, up to error
  CHECK_DELTA( uniform->Radius(), testStar->Radius(), residual );
  CHECK_DELTA( uniform->Mass(),   testStar->Mass(),   residual );

  // test key features identical inside, up to error
  // the n=0 is constant density,
  double const density = testStar->rho(0);
  for(std::size_t i=0; i<LEN; i++){
    CHECK_DELTA( uniform->rad(i), testStar->rad(i), residual );
    CHECK_DELTA( uniform->rho(i), testStar->rho(i), residual );
    CHECK_DELTA( density,     testStar->rho(i), residual );
    CHECK_DELTA( uniform->mr(i),  testStar->mr(i),  residual );
    CHECK_DELTA( uniform->P(i),   testStar->P(i),   residual );
  }

  double xi, theta;
  FILE *fp = fopen("tests/artifacts/exact_0.txt", "w");
  for(std::size_t i=0; i<LEN; i++){
    xi = testStar->getX(i);
    theta = 1. - xi*xi/6.;
    CHECK_DELTA( testStar->getY(i), theta, residual );
    fprintf(fp, "%lf\t%le\n", xi, std::abs(theta-testStar->getY(i)));
  }
  fclose(fp);

  do_test_center (testStar.get(), machine_precision );
  do_test_surface(testStar.get(), expected_error);
}

TEST_CASE("star: star_polytrope_n1") {
  double const INDEX (1.0);
  std::size_t const LEN (1001);
  auto testStar = std::make_unique<Polytrope>(INDEX, LEN);

  // expected error: LEN^{-4} per step, compounded over LEN steps
  double const expected_error = pow(1.0/LEN, 3);
  double const residual = pow(10.0,ceil(log10(testStar->SSR())));
  CHECK_LT( residual, expected_error );

  const double Rn = sqrt( (INDEX+1.0)/(4.0*m_pi) );

  CHECK_DELTA( testStar->getX(LEN-1), m_pi,      residual );
  CHECK_DELTA( testStar->Radius(), Rn*m_pi,      residual );
  CHECK_DELTA( testStar->Mass(), Rn*(INDEX+1.)*m_pi, residual );

  double x, theta;
  FILE *fp = fopen("tests/artifacts/exact_1.txt", "w");
  for(std::size_t i=1; i<LEN; i++){
    x = testStar->rad(i)/Rn;
    theta = sin(x)/x;
    CHECK_DELTA( testStar->P(i),  pow(theta,2),  residual );
    CHECK_DELTA( testStar->rho(i),  theta,     residual );
    CHECK_DELTA( testStar->getY(i), theta,     residual );
    fprintf(fp, "%lf\t%le\n", x, std::abs(theta-testStar->getY(i)));
  }
  fclose(fp);
  
  double acceptable_error(1.e-6);
  do_test_center  ( testStar.get(), acceptable_error );
  do_test_surface ( testStar.get(), acceptable_error );
}

TEST_CASE("star: star_polytrope_n5") {
  // this is not a star, as it has no surface, but extends to infinity
  // it is useful only because it has an exact solution
  // the radius and mass cannot be tested, as they do not exist
  const double INDEX = 5.0;
  const std::size_t LEN = 1001;
  auto testStar = std::make_unique<Polytrope>(INDEX, LEN);

  const double expected_error = pow(1.0/LEN, 3);
  const double residual = pow(10.0,ceil(log10(testStar->SSR())));
  CHECK_LT( residual, expected_error );

  double x, theta;
  FILE *fp = fopen("tests/artifacts/exact_5.txt", "w");
  for(std::size_t i=0; i<LEN; i++){
    x = testStar->getX(i);
    theta = pow( 1.0 + x*x/3.0, -0.5);
    CHECK_DELTA( testStar->getY(i), theta,       residual );
    CHECK_DELTA( testStar->P(i),   pow(theta,INDEX+1), residual );
    CHECK_DELTA( testStar->rho(i), pow(theta,INDEX),   residual );
    fprintf(fp, "%lf\t%le\n", x, std::abs(theta-testStar->getY(i)));
  }
  fclose(fp);

  double acceptable_error(1.e-4);
  do_test_center( testStar.get(), acceptable_error );
  // there is no surface
}

TEST_CASE("star: star_several_polytropes") {
  double const INDEX[] = {2.0, 2.5, 3.0, 3.5, 4.0};
  std::size_t const LEN = 1001;

  double const expected_error = 3.0 * pow(1.0/LEN, 3);
  double const acceptable_error(1.e-4);
  double residual;

  for(double n : INDEX){
    CAPTURE(n);
    auto testStar = std::make_unique<Polytrope>(n, LEN);
    residual = testStar->SSR();
    CHECK_LT(residual, expected_error);
    // assert P/rho = theta?
    do_test_center ( testStar.get(), acceptable_error );
    do_test_surface( testStar.get(), acceptable_error );
  }
}

TEST_CASE("star: star_polytrope_MR_constructor"){
  const double BigM[] = {0.2, 0.6, 1.0, 2.0, 10.0};
  const double BigR[] = {0.1, 1.1, 1.2, 20.0};
  const double index[] = {1.0, 1.5, 3.0};
  const std::size_t LEN = 1001;
  double const expected_error = pow(10./LEN , 4);

  for(double M : BigM){
    for(double R : BigR){
      for(double n : index){
        auto testStar = std::make_unique<Polytrope>(M*MSOLAR, R*REARTH, n, LEN);
        auto refStar = std::make_unique<Polytrope>(n, LEN);
        CHECK_DELTA(M, testStar->Mass()/MSOLAR, 1.e-14);
        CHECK_DELTA(R, testStar->Radius()/REARTH, 1.e-14);
        CHECK_LT(testStar->SSR(), expected_error);
        for(std::size_t i=0; i<LEN; i++){
          CHECK_EQ( testStar->getY(i), refStar->getY(i) );
          CHECK_EQ( testStar->getX(i), refStar->getX(i) );
        }
      }
    }
  }
}

TEST_CASE("star: star_test_make_exact_error_graph [artifacts]") {
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

TEST_CASE("star: star_CHWD_against_chandrasekhar") {
  // Chandrasekhar 1939 table 27 has properties of several WDs
  // to compare must use values of A0, B0 that were in use in 1939

  std::size_t const LEN(2001);
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
    auto testStar = std::make_unique<ChandrasekharWD>(Y0[n], LEN, Chandrasekhar::A01939, Chandrasekhar::B01939);
    M = testStar->Mass()/MSOLAR;
    R = testStar->Radius();
    // assert reasonable numerical error
    CHECK_LT( testStar->SSR(), 1.e-8 );
    // assert agreement with tabular values in Chandrasekhar 1939
    CHECK_EPS( M, M1939[n], 1.e-2);
    CHECK_EPS( R, R1939[n], 1.e-2);
    // write out values to table
    fprintf(fp, "%0.8lf\t%0.2lf\t%0.3lf\t%0.2le\t%0.3le\n", invY0sq[n], M1939[n], M, R1939[n], R);
    
    // test expansions
    do_test_center(testStar.get(), 1.e-4);
    // do_test_surface(testStar.get(), 1.e-3); // TODO why is this fit so bad? 
  }
  fclose(fp);
}

TEST_CASE("star: star_CHWD_grad_constructor") {
  // test the constructor that accepts a gradiant in mu
  // try both a constant gradient, and a sigmoidal gradient
  std::size_t const LEN(1001);
  ChandrasekharWD *testStar;
  const double Y0[] = {1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 2.0};
  double F0;

  for(double y0 : Y0){
    // star with constant mue
    auto testStar = std::make_unique<ChandrasekharWD>(y0, LEN, Chandrasekhar::constant_mu{2.0});
    CHECK_LT(testStar->SSR(), 1.e-4);
    do_test_center(testStar.get(), 1.e-4);
    // TODO make surface test work
    do_test_surface(testStar.get(), 1.0);
    testStar.release();
    // star with sigmoidally varying mue (CHWD++)
    F0 = Chandrasekhar::factor_f(sqrt(y0*y0-1.));
    testStar = std::make_unique<ChandrasekharWD>(y0, LEN, Chandrasekhar::sigmoidal_in_logf{2.,F0,2.,1.});
    CHECK_LT(testStar->SSR(), 1.e-4);
    do_test_center(testStar.get(), 1.e-4);
    // TODO make surface test work
    do_test_surface(testStar.get(), 1.0);
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
  MESA *testStar = new MESA(testfile, LEN);
  // some routine verifications
  CHECK_EQ(testStar->length(), 2*(LEN-1) + 1);
  CHECK_LT(testStar->SSR(), 1.0/LEN);
  // assert both models are equal at the knots to within relative precision
  double const MACHINE_PRECISION = 1.e-15;
  for(std::size_t i=0; i<LEN-2; i++){
    CHECK_EPS(refStar->rad(i),    testStar->rad(2*i),    MACHINE_PRECISION);
    CHECK_EPS(refStar->mr(i),     testStar->mr(2*i),     MACHINE_PRECISION);
    CHECK_EPS(refStar->P(i),    testStar->P(2*i),    MACHINE_PRECISION);
    CHECK_EPS(refStar->rho(i),    testStar->rho(2*i),    MACHINE_PRECISION);
    CHECK_EPS(refStar->dPhidr(i),   testStar->dPhidr(2*i),   MACHINE_PRECISION);
    // these ae calculare values, but should still be within precision
    CHECK_EPS(refStar->getU(i),   testStar->getU(2*i),   MACHINE_PRECISION);
    CHECK_EPS(refStar->getC(i),   testStar->getC(2*i),   MACHINE_PRECISION);
    CHECK_EPS(refStar->getVg(i),  testStar->getVg(2*i),  MACHINE_PRECISION);
    CHECK_EPS(refStar->getAstar(i), testStar->getAstar(2*i), MACHINE_PRECISION);
  }
}

TEST_CASE("star: star_MESA_constructor") {
  // should first calculate SSR using only the knots
  // then check the SSR isn't very much worse
  std::size_t const LEN(101);
  std::string mesa_file = "mesa_co_wd_cold.dat";
  auto testStar = std::make_unique<MESA>(mesa_file, LEN);
  CHECK_LT( testStar->SSR(), 1.0e-2 );
  do_test_center(testStar.get(), 1.e-4);
  // TODO make surface work
  // do_test_surface(testStar.get(), 1.e-6);
}

// TODO make reduced data file

// TODO read in MESA model, previous output, compare

/***** TESTS OF SIMPLE WD MODELS ******/

// TODO create a SWD, save its output
// then create one here, compare to previous

void test_simple_wd_constructor(){
  /* this tests the n=0 polytrope, by comparing it to the exactly-known
  ** solution of an isopycnic (aka uniform-density) star.  The polytrope solution 
  ** should arrive at the isopycnic solution, but through numerical 
  ** integration as opposed to direct setting by an equation.*/
  fprintf(stderr, "\nSTAR TESTS - SIMPLEWD");
  double const INDEX (0.0);
  std::size_t const LEN (1001);
  auto testStar = std::make_unique<SimpleWD>(0.608, 1.2e4, LEN);

  // use the SSR to quantify acceptable error
  double const machine_precision = 10.e-16;
  double const expected_error = 10.0 * pow(1.0/LEN, 4);
  double const residual = pow(10.0,ceil(log10(testStar->SSR())));
  CHECK_LT( residual, expected_error );

  do_test_center (testStar.get(), expected_error);
  do_test_surface(testStar.get(), expected_error);
}

} // end TEST_SUITE
