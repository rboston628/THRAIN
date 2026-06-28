#include "../lib/matrix.h"
#include "doctest.h"
#include <functional>
#include <random>
#include <cmath>
#include <complex>

// TODO: time-profile solutions
// time profile 100 dets and inversions
// time profile det, inversion of one enormous matrix
// have tests requiring certain time performance
// TODO: check singular, sparse matrices

namespace {

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

// Doctest fixture for RNG
struct Rngfixture {
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution;
    std::function<double(void)> rando;
    std::function<std::complex<double>(void)> randoC;

    explicit Rngfixture(unsigned seed = 0) : generator(seed), distribution(-10.0, 10.0) {
        rando = [this]()->double { return this->distribution(this->generator); };
        randoC = [this]()->std::complex<double> { return {distribution(generator),distribution(generator)}; };    
    }
};

} // namespace


TEST_SUITE("Matrix") {

TEST_CASE("det_singular"){
    double constexpr ZERO = 0.0;
    double m[4][4] = {{1,4,6,2},{-7,2,1,4}, {0,0,0,0}, {1,9,8,5}};
    CHECK_EQ(matrix::determinant(m), ZERO);
}

TEST_CASE("det_scalar"){
    double constexpr VALUE = 13.0;
    double m[1][1] = {{VALUE}};
    CHECK_EQ(matrix::determinant(m), VALUE);
}

TEST_CASE("det_2x2"){
    //|  1 3 | -7 2 | -7   2   | -7  2   |
    //| -7 2 |  1 3 |  0 3+2/7 |  0 23/7 |
    //| d=1  | d=-1 | d=7      | d=23    |
    double m[2][2] = {{1,3}, {-7,2}};
    double const RESULT = det_2x2(m);
    CHECK_EQ( matrix::determinant(m), RESULT);
    // check the elements have been changed as expected
    CHECK_EQ(m[0][0], -7.0);
    CHECK_EQ(m[0][1], 2.0);
    // m[1][0] is not changed by the algorithm
    CHECK_EQ(m[1][1], 23./7.);
}

TEST_CASE("det_3x3"){
    double m[3][3] = {
        {1,2,3},
        {3,5,6},
        {8,2,0}
    };
    double const RESULT = det_3x3(m);
    CHECK_DELTA(matrix::determinant(m), RESULT, 1.0e-10);
    // TODO: test shape of matrix
}

TEST_CASE("det_identity") {
    std::size_t constexpr N=5;
    double I[N][N];
    for(int i=0; i<N; i++){
        I[i][i] = 1.0;
    }
    CHECK_EQ(matrix::determinant(I), 1.0);
}

TEST_CASE("det_diagonal") {
    std::size_t constexpr N=5;
    double D[N][N];
    for(int i=0; i<N; i++){
        D[i][i] = 1.0 + i;
    }
    double expected = 1.0;
    for(int i=0; i<N; i++){
        expected *= D[i][i];
    }
    CHECK_EQ(matrix::determinant(D), expected);
}

TEST_CASE("det_row_swap_flips_sign") {
    std::size_t constexpr N=3;
    double A[N][N] = {
        {1, 2, 3},
        {4, 0, 6},
        {7, 8, 9}
    };
    double B[N][N], C[N][N], D[N][N];
    for (int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            B[i][j] = A[i][j];
            C[i][j] = A[i][j];
            D[i][j] = A[i][j];
        }
    }
    // all start equal to A
    // B is A with rows 0 and 1 swapped
    matrix::swap_rows(B, 0, 1);
    // C is A with rows 0 and 1 swapped,  then 1 and 2 swapped
    matrix::swap_rows(C, 0, 1);
    matrix::swap_rows(C, 1, 2);
    // D is A with rows 0 and 2 swapped
    matrix::swap_rows(D, 0, 1); // 1, 0, 2
    matrix::swap_rows(D, 1, 2); // 1, 2, 0
    matrix::swap_rows(D, 0, 2); // 0, 2, 1
    double const detA = matrix::determinant(A);
    double const detB = matrix::determinant(B);
    double const detC = matrix::determinant(C);
    double const detD = matrix::determinant(D);
    CHECK_EQ(detB, -detA);
    CHECK_EQ(detC, -detB);
    CHECK_EQ(detC, detA);
    CHECK_EQ(detD, -detC);
    CHECK_EQ(detD, -detA);
}

TEST_CASE("det_duplicate_rows_zero") {
    double m[3][3] = {
        {1, 2, 3},
        {1, 2, 3},
        {2, 4, 6}
    };
    CHECK_EQ(matrix::determinant(m), 0.0);
}

TEST_CASE_FIXTURE(Rngfixture, "det_many_2x2"){
    const std::size_t N=2, num_trials = 100;
    
    double sum = 0.0;
    double m[N][N];
    for(int I=0; I<num_trials; I++){
        make_matrix<N,N,double>(m, this->rando);
        double ans = det_2x2(m);
        sum += fabs( ans - matrix::determinant(m) )/fabs(ans);
    }
    sum /= num_trials;
    CHECK_LT(sum, 1.e-12);
}

TEST_CASE_FIXTURE(Rngfixture, "det_many_3x3"){
    const std::size_t N = 3, num_trials = 100;

    double sum = 0.0;
    double m[N][N];
    for(int I=0; I<num_trials; I++){
        make_matrix<N,N,double>(m, this->rando);
        double ans = det_3x3(m);
        sum += fabs( ans - matrix::determinant(m) )/fabs(ans);
    }
    sum /= num_trials;
    CHECK_LT(sum, 1.e-12);
}

TEST_CASE_FIXTURE(Rngfixture, "det_many_NxN"){
    const std::size_t N = 8, num_trials = 100;

    double sum = 0.0;
    double m[N][N];
    for(int I=0; I<num_trials; I++){
        make_matrix<N,N,double>(m, this->rando);
        double ans = det_NxN(m);
        sum += fabs( ans - matrix::determinant(m) )/fabs(ans);
    }
    sum /= num_trials;
    CHECK_LT(sum, 1.e-12);
}

TEST_CASE_FIXTURE(Rngfixture, "det_NxN_complex"){
    const int N=8;
    typedef std::complex<double> C;

    C m[N][N];
    make_matrix<N,N,C>(m, this->randoC);
    C ans = det_NxN(m);
    C res = matrix::determinant(m);
    CHECK_LT( abs(ans - res)/abs(ans), 1.e-12 ); 
}

TEST_CASE_FIXTURE(Rngfixture, "det_many_NxN_complex"){
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
    CHECK_LT(sum, 1.e-12);
}

TEST_CASE_FIXTURE(Rngfixture,"invert_2x2"){
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
        CHECK_LT( fabs(ans[i] - res[i])/fabs(ans[i]),  1.e-12 ); 
    }
    CHECK(matrix_triangular(m));
}

TEST_CASE_FIXTURE(Rngfixture, "invert_many_2x2"){
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
        CHECK(matrix_triangular(m));
    }
    sum /= N;
    sum /= num_trials;
    CHECK_LT( sum, 1.e-12 );
}

TEST_CASE_FIXTURE(Rngfixture, "invert_many_NxN"){
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
        CHECK(matrix_triangular(m));
    }
    sum /= N;
    sum /= num_trials;
    CHECK_LT( sum, 1.e-12 );
}

TEST_CASE_FIXTURE(Rngfixture, "invert_many_NxN_complex"){
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
        CHECK(matrix_triangular(m));
    }
    sum /= N;
    sum /= num_trials;
    CHECK_LT( sum, 1.e-12 );
}


TEST_CASE_FIXTURE(Rngfixture, "invert_tridiagonal"){
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
            // CHECK_EQ( ans[i], res[i] );
            sum += fabs(ans[i] - res[i])/fabs(ans[i]); 
        }
        sum /= N;
        CHECK_LT( sum, 1.e-12 );
    }
}

TEST_CASE_FIXTURE(Rngfixture, "invert_identity_returns_b"){
    int constexpr N=5, num_trials = 100;
    double sum = 0.0;
    double I[N][N] = {0.0}, res[N], b[N];
    for(int t=0; t<num_trials; t++){
        // re-form the identity, which is destroyed by the algorithm, and make a new random b and ans
        for(int i=0; i<N; i++) I[i][i] = 1.0;
        for(int i=0; i<N; i++) b[i] = this->rando();
        double res[N];
        CHECK_EQ(matrix::invertMatrix(I, b, res), 0);
        for(int i=0; i<N; i++){
            CHECK_EQ(b[i], res[i]);
        }
        CHECK(matrix_triangular(I));
    }
}

TEST_CASE("invertMatrix_singular_should_fail") {
    double A[3][3] = {
        {1,2,3},
        {1,2,3}, // repeat = singular
        {4,5,6}
    };
    double b[3] = {1,2,3};
    double x[3];

    // Singular matrix shoud fail and return 1
    ThrainLogger::setLogLevel(ThrainLogger::LogLevel::MUTE); // suppress error message
    CHECK_EQ(matrix::invertMatrix(A, b, x), 1);
    ThrainLogger::setLogLevel(ThrainLogger::LogLevel::INFO); // restore error message
}

TEST_CASE("invertTridiagonal_small_N") {
    {   // N=1
        const std::size_t N = 1;
        double left[N] = {0};
        double diag[N] = {2};
        double rite[N] = {0};
        double rhs[N]  = {6};
        double x[N]    = {0};
        CHECK_EQ(matrix::invertTridiagonal(left, diag, rite, rhs, x, N), 0);
        CHECK_EQ(x[0], 3.0);
    }
    {   // N=2
        const std::size_t N = 2;
        // [2 -1; -1 2] x = rhs, with x = [1, 2] => rhs = [0, 3]
        double left[N] = {0, -1};
        double diag[N] = {2, 2};
        double rite[N] = {-1, 0};
        double rhs[N]  = {0, 3};
        double x[N]    = {0, 0};

        CHECK_EQ(matrix::invertTridiagonal(left, diag, rite, rhs, x, N), 0);
        CHECK_DELTA(x[0], 1.0, 1e-12);
        CHECK_DELTA(x[1], 2.0, 1e-12);
    }
    {   // N=3
        const std::size_t N = 3;
        // [2 -1 0 ] [1] = [1] 
        // [-1 2 -1] [1] = [0]
        // [0 -1  2] [1] = [1]
        // [2 -1 0; -1 2 -1; 0 -1 2] x = rhs, with x = [1, 1, 1] => rhs = [1, 0, 1]
        double left[N] = {0, -1, -1};
        double diag[N] = {2, 2, 2};
        double rite[N] = {-1, -1, 0};
        double rhs[N]  = {1, 0, 1};
        double x[N]    = {0, 0, 0};

        CHECK_EQ(matrix::invertTridiagonal(left, diag, rite, rhs, x, N), 0);
        CHECK_DELTA(x[0], 1.0, 1e-12);
        CHECK_DELTA(x[1], 1.0, 1e-12);
        CHECK_DELTA(x[2], 1.0, 1e-12);
    }
}

} // test suite

