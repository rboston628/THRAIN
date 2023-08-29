#include "../lib/Splinor.h"
#include "../src/constants.h"
#include <cxxtest/TestSuite.h>
#include <random>

class SplinorTest : public CxxTest::TestSuite {
public:

void test_splinor_basic(){
    printf("SPLINOR BASIC\n");
    double xtest[5] = {0, 1, 2, 3, 4};
    double ytest[5] = {0, 1, 2, 3, 4};
    
    Splinor line(xtest,ytest,5);

    double dx=0.01;
    for(double x=0.0; x<4.0; x+=dx){
        TS_ASSERT_EQUALS(x, line(x));
        TS_ASSERT_EQUALS(1.0, line.deriv(x));
    }

    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0,4.0);

    double x=0.0;
    for(int I=0; I<100; I++){
        x = distribution(generator);
        TS_ASSERT_EQUALS(x, line(x));
        TS_ASSERT_EQUALS(1.0, line.deriv(x));
    }
}

void test_splinor_quadratic(){
    // printf("SPLINOR QUADRATIC\n");

    // double A=1.2, B=0.5, C=3.0;
    // std::function<double(double)> quadratic = [A,B,C](double x)->double {return (A*x+B)*x+C;};
    
    // double x0 = 0.0, x1 = 10.0, deltax;

    // size_t const num_data[4] = {10, 100, 1000, 10000};
    // double resid[4], *xtest, *ytest;
    // for(int I=0; I<4; I++){
    //     xtest = new double[num_data[I]];
    //     ytest = new double[num_data[I]];
    //     deltax = (x1-x0)/num_data[I];
    //     for(int i=0; i<num_data[I]; i++){
    //         xtest[i] = i*deltax;
    //         ytest[i] = quadratic(xtest[i]);
    //     }
    //     Splinor spline(xtest,ytest,num_data[I]);
    //     double dx = deltax/7.0;
    //     double sum_fit = 0.0;
    //     size_t count = 0;
    //     for(double x=1.0; x<9.0; x+=dx){
    //         sum_fit += fabs(quadratic(x)-spline(x))/fabs(quadratic(x));
    //         count++;
    //     }
    //     sum_fit /= count;
    //     resid[I] = sum_fit;
    // }

    // for(int I=1; I<4; I++){
    //     TS_ASSERT_LESS_THAN( resid[I], fmax(1.e-15,resid[I-1]) );
    // }
    // TS_ASSERT_LESS_THAN(resid[3], 1.e-12);
    
    // std::default_random_engine generator;
    // std::uniform_real_distribution<double> distribution(0.0,10.0);

    // double x=0.0;
    // for(int I=0; I<100; I++){
    //     x = distribution(generator);
    //     TS_ASSERT_EQUALS(quadratic(x), spline(x));
    //     TS_ASSERT_EQUALS(2.0*A*x+B, spline.deriv(x));
    // }
}

void test_splinor_cubic(){
    // printf("SPLINOR CUBIC\n");
    
    // double A=1.2, B=0.5, C=3.0, D=-5.0;
    // std::function<double(double)> cubic = [A,B,C,D](double x)->double {return ((A*x+B)*x+C)*x+D;};
    
    // int const num_data = 100;
    // double xtest[num_data], ytest[num_data];
    // double xstart = 0.0, xstop = 10.0;
    // double deltax = (xstop-xstart)/num_data;

    // for(int i=0; i<num_data; i++){
    //     xtest[i] = i*deltax;
    //     ytest[i] = cubic(xtest[i]);
    // }

    // Splinor spline(xtest,ytest,num_data);

    // double dx=deltax/7.0;
    // double sum_fit=0.0, sum_deriv=0.0;
    // size_t count=0;
    // for(double x=1.0; x<9.0; x+=dx){
    //     sum_fit += fabs(cubic(x) - spline(x))/fabs(cubic(x));
    //     sum_deriv += fabs((3.*A*x+2.*B)*x+C - spline(x))/fabs((3.*A*x+2.*B)*x+C);
    //     count++;
    // }
    // sum_fit /= count;
    // sum_deriv /= count;
    // TS_ASSERT_LESS_THAN(sum_fit, 1.e-8);
    // TS_ASSERT_LESS_THAN(sum_deriv, 1.e-8);

    // std::default_random_engine generator;
    // std::uniform_real_distribution<double> distribution(0.0,10.0);

    // double x=0.0;
    // for(int I=0; I<100; I++){
    //     x = distribution(generator);
    //     TS_ASSERT_EQUALS(quadratic(x), spline(x));
    //     TS_ASSERT_EQUALS(2.0*A*x+B, spline.deriv(x));
    // }
}

};

