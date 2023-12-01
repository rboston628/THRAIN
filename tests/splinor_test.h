#include "../lib/Splinor.h"
#include "../src/constants.h"
#include <cxxtest/TestSuite.h>
#include <random>
#include <stdio.h>

class SplinorTest : public CxxTest::TestSuite {
public:

std::default_random_engine generator;

SplinorTest() {
    // remove old artifacts
    system("rm -r    ./tests/artifacts/spline");
    system("mkdir -p ./tests/artifacts/spline");
    // prepare a random number generator for tests
    generator = std::default_random_engine();
}

void test_slinor_basic(){
    printf("SPLINE BASIC\n");

    // range of test data
    double x0=0.0, x1=1.0;

    // create the data -- sixth-order polynomial of random coefficients
    std::size_t const data_order = 6;
    // generate random coefficients
    double y0=-5.0, y1=7.2;
    std::uniform_real_distribution<double> coefficient_maker(y0,y1);
    double A[data_order];
    for(std::size_t i=0; i<data_order; i++){
        A[i] = coefficient_maker(this->generator);
    }
    std::function<double(double)> polynomial = [A, data_order](double x) -> double {
        double y = A[0];
        for (std::size_t i=1; i < data_order; i++){
            y += A[i] * pow(x, i);
        }
        return y;
    };
    // create the x and y data
    std::size_t const N = 15;
    double delta = (x1-x0)/(N-1);
    double xtest[N], ytest[N];
    for(std::size_t k=0; k<N; k++){
        xtest[k] = k*delta;
        ytest[k] = polynomial(xtest[k]);
    }

    // create each sort of spline for the data
    Splinor spline_natural(xtest,ytest,N, Splinor::bc::NATURAL);
    Splinor spline_quadratic(xtest,ytest,N, Splinor::bc::QUADRATIC);
    double yprime0 = A[0];
    double yprime1 = 0.0;
    for(std::size_t i=1; i<data_order; i++){
        yprime1 += i * A[i] * pow(x1, i-1);
    }
    Splinor spline_clamped(xtest,ytest,N, Splinor::bc::CLAMPED, yprime0, yprime1);
    
    // ensure that the splines fit the data
    for(std::size_t k=0; k<N; k++){
        TS_ASSERT_EQUALS( spline_natural(xtest[k]), ytest[k] );
        TS_ASSERT_EQUALS( spline_quadratic(xtest[k]), ytest[k] );
        TS_ASSERT_EQUALS( spline_clamped(xtest[k]), ytest[k] );
    }

    // ensure the derivatives at the end of the clamped spline match what they should
    TS_ASSERT_DELTA( spline_clamped.deriv(x0), yprime0, 1.e-14);
    TS_ASSERT_DELTA( spline_clamped.deriv(x1), yprime1, 1.e-14);

    // ensure that the third-order term (a) of quadratic spline is zerp
    double a, b, c, d;
    std::tie(a, b, c, d) = spline_quadratic.getCoefficients(0);
    TS_ASSERT_DELTA( a, 0.0, 1.e-15 );
    std::tie(a, b, c, d) = spline_quadratic.getCoefficients(N-2);
    TS_ASSERT_DELTA( a, 0.0, 1.e-15);

    // test that the splines interpolate the data
    std::size_t count=0;
    double dx = delta/7.0;
    double sum_natural = 0., sum_quadratic = 0., sum_clamped = 0.;
    for(double x = x0; x <= x1; x += dx){
        double y = polynomial(x);
        sum_natural   += fabs(y - spline_natural(x))/fabs(y);
        sum_clamped   += fabs(y - spline_clamped(x))/fabs(y);
        sum_quadratic += fabs(y - spline_quadratic(x))/fabs(y);
        count++;
    }
    sum_natural /= count;
    sum_clamped /= count;
    sum_quadratic /= count;
    TS_ASSERT_LESS_THAN( sum_natural, 0.01 );
    TS_ASSERT_LESS_THAN( sum_clamped, 0.01 );
    TS_ASSERT_LESS_THAN( sum_quadratic, 0.01 );
}

void test_splinor_line_natural(){
    printf("SPLINOR LINE NATURAL\n");

    // test distribution, line of slope 1
    std::size_t const num_points = 10;
    double xtest[num_points], ytest[num_points];
    for(std::size_t i=0; i<num_points; i++){
        xtest[i] = ytest[i] = i;
    }    

    // spline fit the line
    Splinor spline(xtest,ytest,num_points);

    // run through the x-axis and check the fit
    // a natural cubic spline can exactly fit a line
    double a,b,c,d;
    double dx=0.01;
    for(double x=0.0; x<=num_points; x+=dx){
        TS_ASSERT_EQUALS(x, spline(x));
        TS_ASSERT_EQUALS(1.0, spline.deriv(x));   
    }

    // ensure the coefficients are correct     
    for (std::size_t i=0; i<num_points-1; i++){
        std::tie(a, b, c, d) = spline.getCoefficients(i);
        TS_ASSERT_EQUALS(a, 0.0);
        TS_ASSERT_EQUALS(b, 0.0);
        TS_ASSERT_EQUALS(c, 1.0);
        TS_ASSERT_EQUALS(d, i);
    }

    // pick random spots on the x-axis, and ensure exact fit
    std::uniform_real_distribution<double> sample_x_axis(0.0,double(num_points));
    double x=0.0;
    for(int I=0; I<100; I++){
        x = sample_x_axis(this->generator);
        TS_ASSERT_EQUALS(x, spline(x));
        TS_ASSERT_EQUALS(1.0, spline.deriv(x));
    }
}

void test_splinor_line_clamped(){
    printf("SPLINOR LINE CLAMPED\n");

    // test distribution, line of slope 1
    std::size_t const num_points = 10;
    std::uniform_real_distribution<double> sample_x_axis(0.0,double(num_points));
    double m = sample_x_axis(this->generator);
    double xtest[num_points], ytest[num_points];
    for(std::size_t i=0; i<num_points; i++){
        xtest[i] = i;
        ytest[i] = i * m;
    }    

    // spline fit the line
    Splinor spline(xtest, ytest, num_points, Splinor::bc::CLAMPED, m, m);

    // run through the x-axis and check the fit
    // a natural cubic spline can exactly fit a line
    double dx=0.01;
    for(double x=0.0; x<=num_points; x+=dx){
        TS_ASSERT_DELTA(m * x, spline(x), 1.e-14);
        TS_ASSERT_DELTA(m, spline.deriv(x), 1.e-14);   
    }

    // ensure the coefficients are correct    
    double a,b,c,d; 
    for (std::size_t i=0; i<num_points-1; i++){
        std::tie(a, b, c, d) = spline.getCoefficients(i);
        TS_ASSERT_DELTA(a, 0.0, 1.e-14);
        TS_ASSERT_DELTA(b, 0.0, 1.e-14);
        TS_ASSERT_DELTA(c, m, 1.e-14);
        TS_ASSERT_DELTA(d, i * m, 1.e-14);
    }

    // pick random spots on the x-axis, and ensure exact fit
    double x=0.0;
    for(int I=0; I<100; I++){
        x = sample_x_axis(this->generator);
        TS_ASSERT_DELTA(m * x, spline(x), 1.e-14);
        TS_ASSERT_DELTA(m, spline.deriv(x), 1.e-14);
    }
}

void test_splinor_quadratic(){
    printf("SPLINOR QUADRATIC\n");

    double x0 = 0.0, x1 = 10.0;
    std::uniform_real_distribution<double> sample_x_axis(x0,x1);

    // number of data points, for testing error scaling
    std::size_t const num_data[4] = {10, 20, 40, 80};

    // the quadratic function to be fit
    double A=1.2, B=0.7, C=2.3;
    std::function<double(double)> quadratic = [A,B,C](double x)->double {return (A*x+B)*x+C;};

    double resid[4], *xtest, *ytest, deltax;
    for(int I=0; I<4; I++){
        // create the fit data
        xtest = new double[num_data[I]+1];
        ytest = new double[num_data[I]+1];
        deltax = (x1-x0)/(num_data[I]);
        for(int i=0; i<num_data[I]+1; i++){
            xtest[i] = i*deltax;
            ytest[i] = quadratic(xtest[i]);
        }
        // create the spline fit -- use quadratic boundary conditions
        Splinor spline(xtest,ytest,num_data[I], Splinor::bc::QUADRATIC);

        // now test the interpolaton
        // create a summed residual between the fit and function
        double dx = deltax/8.0;
        double sum_fit = 0.0;
        std::size_t count = 0;
        for(double x=1.0; x<9.0; x+=dx){
            TS_ASSERT_DELTA( quadratic(x), spline(x), 1.e-12 );
            sum_fit += fabs(quadratic(x)-spline(x))/fabs(quadratic(x));
            count++;
        }
        sum_fit /= count;
        resid[I] = sum_fit;

        // verify the coefficients
        double a,b,c,d;
        for (std::size_t i=2; i<num_data[I]-2; i++) {
            std::tie(a, b, c, d) = spline.getCoefficients(i);
            TS_ASSERT_DELTA(a, 0.0, 1.e-10);
            TS_ASSERT_DELTA(b,   A, 1.e-10);
            TS_ASSERT_DELTA(c,   2.*A*xtest[i]+B, 1.e-10);
        }

        delete xtest;
        delete ytest;
    }

    // assert that the summed residual is decreasing with number of points 
    for(int I=1; I<4; I++){
        TS_ASSERT_LESS_THAN( resid[I], fmax(1.e-15,resid[I-1]) );
    }
}

void test_splinor_cubic(){
    printf("SPLINOR CUBIC\n");
    
    double A=1.2, B=0.5, C=3.0, D=-5.0;
    std::function<double(double)> cubic = [A,B,C,D](double x)->double {return ((A*x+B)*x+C)*x+D;};
    
    int const num_data = 100;
    double xtest[num_data], ytest[num_data];
    double xstart = 0.0, xstop = 10.0;
    double deltax = (xstop-xstart)/(num_data-1);

    for(int i=0; i<num_data; i++){
        xtest[i] = i*deltax;
        ytest[i] = cubic(xtest[i]);
    }

    Splinor spline(xtest,ytest,num_data);

    double dx=deltax/7.0;
    double sum_fit=0.0, sum_deriv=0.0;
    std::size_t count=0;
    for(double x=1.0; x<9.0; x+=dx){
        sum_fit += fabs(cubic(x) - spline(x))/fabs(cubic(x));
        double yprime = (3.*A*x+2.*B)*x+C;
        sum_deriv += fabs(yprime - spline.deriv(x))/fabs(yprime);
        count++;
    }
    sum_fit /= count;
    sum_deriv /= count;
    TS_ASSERT_LESS_THAN(sum_fit, 1.e-8);
    TS_ASSERT_LESS_THAN(sum_deriv, 1.e-8);
}

void test_graphical_fits(){
    /** This test utilizes several test data files to fit cubic splines
     * and creates graphs of the spline fits for visual validation.
    */

   printf("SPLINOR GRAPHS\n");
   std::string tests[] = {"log", "cosine", "runge", "step", "window"};
   std::size_t N;
   double *xdata, *ydata;

    for (std::string test : tests) {
        
        // open the data file and read in the number from the first row
        std::string name = "./tests/inputs/spline/test_" + test + ".dat";
        FILE *data_file = fopen(name.c_str(), "r");
        fscanf(data_file, "%lu %*[^\n]", &N);
        
        // create data arrays to hold the data
        xdata = new double[N];
        ydata = new double[N];

        // now read in the data from the file
        for(std::size_t i=0; i<N; i++){
            fscanf(data_file, "%lf\t%lf\n", &xdata[i], &ydata[i]);
        }
        fclose(data_file);

        std::size_t N_sub = 10;
        std::size_t NN = N_sub*(N-1)+1;

        // create the splinor
        Splinor spline(xdata, ydata, N);

        // split out the interpolated data
        std::string outfile = "./tests/artifacts/spline/compare.txt";
        FILE *out = fopen(outfile.c_str(), "w");
        double dx = (xdata[N-1]-xdata[0])/(NN-2);
        for(double x = xdata[0]; x<= xdata[N-1]; x += dx){
            fprintf(out, "%lf\t%lf\n", x, spline(x));
        }
        fclose(out);
        
    	//now graph the results
        std::string output = "./tests/artifacts/spline/test_" + test + ".png";
        FILE *gnuplot = popen("gnuplot -persist", "w");
        fprintf(gnuplot, "reset\n");
        fprintf(gnuplot, "set term png size 1600,800\n");
        fprintf(gnuplot, "set output '%s'\n", output.c_str());
        fprintf(gnuplot, "set title 'Interpolation For %s'\n", test.c_str());
        fprintf(gnuplot, "set xlabel 'x'\n");
        //plot interpolations with lines
        fprintf(gnuplot, "plot '%s' u 1:2 w l t 'spline'", outfile.c_str());
        //plot the data with points
        fprintf(gnuplot, "   , '%s' u 1:2 w p ls 4 t 'data'", name.c_str());
        fprintf(gnuplot, "\n");	
        pclose(gnuplot);

        // clean up
        system("rm ./tests/artifacts/spline/compare.txt"); 
        delete[] xdata;
        delete[] ydata; 
    }
}


};

