#include "../lib/matrix.h"
#include <cxxtest/TestSuite.h>
#include <functional>
#include <random>
#include <cmath>

// TODO: time-profile solutions
// time profile 100 dets and inversions
// time profile det, inversion of one enormous matrix
// have tests requiring certain time performance
// TODO: check singular, sparse matrices

// create a matrix using some function to fill in elements
template<std::size_t N, std::size_t M, class T>
void make_matrix(T (&m)[N][M], std::function<T(void)> make){
    for(int i=0; i<N; i++){
        for(int j=0; j<M; j++){
            m[i][j] = make();
        }
    }
}

/* these functions provide a source of verification
*  they are conveptually intuitive but slower */

// simple 2x2 determinant using familiar cross-product rule
template<class T>
T det_2x2(T const (&m)[2][2]){
    return m[0][0]*m[1][1]-m[0][1]*m[1][0];
}

// simple 3x3 determinant using familiar cross-product rule
template<class T>
T det_3x3(T const (&m)[3][3]){
    T M00[2][2] = {{m[1][1],m[1][2]},{m[2][1],m[2][2]}};
    T M01[2][2] = {{m[1][0],m[1][2]},{m[2][0],m[2][2]}};
    T M02[2][2] = {{m[1][0],m[1][1]},{m[2][0],m[2][1]}};
    return m[0][0]*det_2x2(M00) - m[0][1]*det_2x2(M01) + m[0][2]*det_2x2(M02);
}

// a specialization of the more general NxN to a 3x3
template<class T>
T det_NxN(T const (&m)[3][3]){
    return det_3x3<T>(m);
}

// recursive process for general NxN determinant
// this is a terrible way to find the determinant for real purposes
// being used here only for a more sound verification
template<std::size_t N, class T>
T det_NxN(T const (&m)[N][N]){
    T det = 0.0;
    T minor[N-1][N-1];
    
    double sign = 1.0;
    for(int k=0; k<N; k++){
        // create the minor matrix
        for(int i=1; i<N; i++){
            for(int j=0, ji=0; j<N && ji<N-1; j++){
                if(j!=k){
                    minor[i-1][ji] = m[i][j];
                    ji++;
                }
            }
        }
        det += sign*m[0][k]*det_NxN(minor);
        sign = -sign;
    }
    return det;
}

// verify that a matrix has been put into upper-triangular form
// the matrix inversion method will leave matrix in row-echelon form
template<std::size_t N, class T>
bool matrix_triangular(T const (&m)[N][N]){
    bool upper = true;
    for(int i=0; i<N && upper; i++){
        for(int j=0; j<i && upper; j++){
            upper &= (m[i][j]==T(0));
        }
        //further ensure the diagonal is not zero
        upper &= (m[i][i]!=T(0));
    }
    return upper;
}

class MatrixTest : public CxxTest::TestSuite {
public:

std::default_random_engine generator;
std::uniform_real_distribution<double> distribution;
std::function<double(void)> rando;
std::function<std::complex<double>(void)> randoC;

MatrixTest() {
    generator = std::default_random_engine();
    distribution = std::uniform_real_distribution<double>(-10.0,10.0);
    rando = [this]()->double {
        return this->distribution(this->generator);
    };
    randoC = [this]()->std::complex<double> {
        return std::complex<double>{this->distribution(this->generator),this->distribution(this->generator)};
    };    
}

void test_det_singular(){
    const double ZERO = 0.0;
    double m[4][4] = {{1,4,6,2},{-7,2,1,4}, {0,0,0,0}, {1,9,8,5}};
    TS_ASSERT_EQUALS(matrix::determinant(m), ZERO);
}

void test_det_scalar(){
    const double VALUE = 13.0;
    double m[1][1] = {{VALUE}};
    TS_ASSERT_EQUALS(matrix::determinant(m), VALUE);
}

void test_det_2x2(){
    //|  1 3 | -7 2 | -7   2   | -7  2   |
    //| -7 2 |  1 3 |  0 3+2/7 |  0 23/7 |
    //| d=1  | d=-1 | d=7      | d=23    |
    double m[2][2] = {{1,3},{-7,2}};
    const double RESULT = det_2x2(m);
    TS_ASSERT_EQUALS( matrix::determinant(m), RESULT);
    // check the elements have been changed as expected
    TS_ASSERT_EQUALS(m[0][0],-7.0);
    TS_ASSERT_EQUALS(m[0][1],2.0);
    // m[1][0] is not changed by the algorithm
    TS_ASSERT_EQUALS(m[1][1],23./7.);
}

void test_det_3x3(){
    double m[3][3] = {
        {1,2,3},
        {3,5,6},
        {8,2,0}
    };
    const double RESULT = det_3x3(m);
    TS_ASSERT_DELTA(matrix::determinant(m), RESULT, 1.0e-10);
    // TODO: test shape of matrix
}

void test_many_2x2(){
    const std::size_t N=2, num_trials = 100;
    
    double sum = 0.0;
    double m[N][N];
    for(int I=0; I<num_trials; I++){
        make_matrix<N,N,double>(m, this->rando);
        double ans = det_2x2(m);
        sum += fabs( ans - matrix::determinant(m) )/fabs(ans);
    }
    sum /= num_trials;
    TS_ASSERT_LESS_THAN(sum, 1.e-12);
}

void test_det_many_3x3(){
    const std::size_t N = 3, num_trials = 100;

    double sum = 0.0;
    double m[N][N];
    for(int I=0; I<num_trials; I++){
        make_matrix<N,N,double>(m, this->rando);
        double ans = det_3x3(m);
        sum += fabs( ans - matrix::determinant(m) )/fabs(ans);
    }
    sum /= num_trials;
    TS_ASSERT_LESS_THAN(sum, 1.e-12);
}

void test_det_many_NxN(){
    const std::size_t N = 8, num_trials = 100;

    double sum = 0.0;
    double m[N][N];
    for(int I=0; I<num_trials; I++){
        make_matrix<N,N,double>(m, this->rando);
        double ans = det_NxN(m);
        sum += fabs( ans - matrix::determinant(m) )/fabs(ans);
    }
    sum /= num_trials;
    TS_ASSERT_LESS_THAN(sum, 1.e-12);
}

void test_det_NxN_complex(){
    const int N=8;
    typedef std::complex<double> C;

    C m[N][N];
    make_matrix<N,N,C>(m, this->randoC);
    C ans = det_NxN(m);
    C res = matrix::determinant(m);
    TS_ASSERT_LESS_THAN( abs(ans - res)/abs(ans), 1.e-12 ); 
}

void test_det_many_NxN_complex(){
    const int N=4, num_trials = 100;
    typedef std::complex<double> C;

    double sum = 0.0;
    C m[N][N];
    for(int I=0; I<num_trials; I++){
        make_matrix<N,N,C>(m, this->randoC);
        C ans = det_NxN(m);
        sum += abs( ans - matrix::determinant(m) )/abs(ans);
    }
    sum /= num_trials;
    TS_ASSERT_LESS_THAN(sum, 1.e-12);
}

void test_invert_2x2(){
    const int N=2;

    // solving A*x = b.
    // specify A, x, solve for b.
    // then invert to solve for x and compare.
    double m[N][N], ans[N], b[N];
    make_matrix<N,N,double>(m, this->rando);
    for(int i=0; i<N; i++) ans[i] = this->rando();
    for(int i=0; i<N; i++){
        b[i] = 0.0;
        for(int j=0; j<N; j++){
            b[i] += m[i][j]*ans[j];
        }
    }
    double res[N];
    matrix::invertMatrix(m,b,res);
    for(int i=0; i<N; i++){
        TS_ASSERT_LESS_THAN( fabs(ans[i] - res[i])/fabs(ans[i]),  1.e-12 ); 
    }
    TS_ASSERT(matrix_triangular(m));
}

void test_invert_many_2x2(){
    const int N=2, num_trials = 100;

    // solving A*x = b.
    // specify A, x, solve for b.
    // then invert to solve for x and compare.
    double sum = 0.0;
    double m[N][N], ans[N], b[N];
    make_matrix<N,N,double>(m, this->rando);
    for(int I=0; I<num_trials; I++){
        for(int i=0; i<N; i++) ans[i] = this->rando();
        for(int i=0; i<N; i++){
            b[i] = 0.0;
            for(int j=0; j<N; j++){
                b[i] += m[i][j]*ans[j];
            }
        }
        double res[N];
        matrix::invertMatrix(m,b,res);
        for(int i=0; i<N; i++){
            sum += fabs(ans[i] - res[i])/fabs(ans[i]); 
        }
        TS_ASSERT(matrix_triangular(m));
    }
    sum /= N;
    sum /= num_trials;
    TS_ASSERT_LESS_THAN( sum, 1.e-12 );
}

void test_invert_many_NxN(){
    const int N=6, num_trials = 100;

    // solving A*x = b.
    // specify A, x, solve for b.
    // then invert to solve for x and compare.
    double sum = 0.0;
    double m[N][N], ans[N], b[N];
    make_matrix<N,N,double>(m, this->rando);
    for(int I=0; I<num_trials; I++){
        for(int i=0; i<N; i++) ans[i] = this->rando();
        for(int i=0; i<N; i++){
            b[i] = 0.0;
            for(int j=0; j<N; j++){
                b[i] += m[i][j]*ans[j];
            }
        }
        double res[N];
        matrix::invertMatrix(m,b,res);
        for(int i=0; i<N; i++){
            sum += fabs(ans[i] - res[i])/fabs(ans[i]); 
        }
        TS_ASSERT(matrix_triangular(m));
    }
    sum /= N;
    sum /= num_trials;
    TS_ASSERT_LESS_THAN( sum, 1.e-12 );
}

void test_invert_many_NxN_complex(){
    const int N=6, num_trials = 10;
    typedef std::complex<double> C;

    // solving A*x = b.
    // specify A, x, calculate b.
    // then invert to solve for x from A,b and compare.
    double sum = 0.0;
    C m[N][N], ans[N], b[N];
    make_matrix<N,N,C>(m, this->randoC);
    for(int I=0; I<num_trials; I++){
        for(int i=0; i<N; i++) ans[i] = this->randoC();
        for(int i=0; i<N; i++){
            b[i] = 0.0;
            for(int j=0; j<N; j++){
                b[i] += m[i][j]*ans[j];
            }
        }
        C res[N];
        matrix::invertMatrix(m,b,res);
        for(int i=0; i<N; i++){
            sum += abs(ans[i] - res[i])/abs(ans[i]); 
        }
        TS_ASSERT(matrix_triangular(m));
    }
    sum /= N;
    sum /= num_trials;
    TS_ASSERT_LESS_THAN( sum, 1.e-12 );
}


void test_invert_tridiagonal(){
    const int N=20, num_trials = 100;

    // solving T*x = b where T is tridiagonal
    // specify each diagonal, x, solve for b.
    // then invert to solve for x and compare.
    double sum = 0.0;
    double left[N], diag[N], rite[N];
    double ans[N], b[N];
    for(int I=0; I<num_trials; I++){
        for(int i=0; i<N; i++) left[i] = this->rando();
        for(int i=0; i<N; i++) rite[i] = this->rando();
        for(int i=0; i<N; i++) diag[i] = 2.0*(left[i]+rite[i]);
        for(int i=0; i<N; i++) ans[i] = this->rando();
        left[0] = rite[N-1] = 0.0;
        b[0] = diag[0]*ans[0] + rite[0]*ans[1];
        for(int i=1; i<N-1; i++){
            b[i] = left[i]*ans[i-1] + diag[i]*ans[i] + rite[i]*ans[i+1];
        }
        b[N-1] = left[N-1]*ans[N-2] + diag[N-1]*ans[N-1];
        double res[N];
        matrix::invertTridiagonal(left,diag,rite,b,res, N);
        sum = 0.0;
        for(int i=0; i<N; i++){
            // TS_ASSERT_EQUALS( ans[i], res[i] );
            sum += fabs(ans[i] - res[i])/fabs(ans[i]); 
        }
        sum /= N;
        TS_ASSERT_LESS_THAN( sum, 1.e-12 );
    }
}
};

