#include "../lib/rootfind.h"
#include "../src/constants.h"
#include "doctest.h"
#include <random>
#include <set>

namespace {

struct RootFindLineFixture {
    double root = 2.0;
    double over_root = 2.4;
    double under_root = 1.6;
};
struct Line {
    virtual ~Line() = default;
    virtual double operator()(double x) const = 0;
};
struct RootFindLineFlat : RootFindLineFixture, Line {
    double operator()(double x) const override final{ return 1.0; }
} flat;
struct RootFindLineUp : RootFindLineFixture, Line {
    double operator()(double x) const override final{ return (x-root); }
} slopeup;
struct RootFindLineDown : RootFindLineFixture, Line {
    double operator()(double x) const override final{ return (root-x); }
} slopedown;
struct RootFindSinusoid : RootFindLineFixture, Line {
    double operator()(double x) const override final{ return sin( x * m_pi/root ); }
} sinusoid;
struct RootFindParabola : RootFindLineFixture, Line {
    double operator()(double x) const override final{ return (x-root)*(x-root) + 1.0; }
} parabola;

struct RootFind2DFixture {
    double xroot[2] = {1.5, 3.7};
};
struct Field {
    virtual ~Field() = default;
    virtual void operator()(double f[2], double x[2]) const = 0;
};
struct RootFind2DUpUp : RootFind2DFixture, Field {
    void operator()(double f[2], double x[2]) const override final{
        f[0] = (x[0] - xroot[0]);
        f[1] = (x[1] - xroot[1]);
    }
} upup;
struct RootFind2DUpDown : RootFind2DFixture, Field {
    void operator()(double f[2], double x[2]) const override final{
        f[0] = (x[0] - xroot[0]);
        f[1] =-(x[1] - xroot[1]);
    }
} updown;
struct RootFind2DDownUp : RootFind2DFixture, Field {
    void operator()(double f[2], double x[2]) const override final{
        f[0] =-(x[0] - xroot[0]);
        f[1] = (x[1] - xroot[1]);
    }
} downup;
struct RootFind2DDownDown : RootFind2DFixture, Field {
    void operator()(double f[2], double x[2]) const override final{
        f[0] =-(x[0] - xroot[0]);
        f[1] =-(x[1] - xroot[1]);
    }
} downdown;
} // namespace

TEST_SUITE("RootFind [unit]") {

  TEST_CASE("rootfind: pseudo_unif_not_repeated") {
    /** 
     * make sure the pseudo_unif produces a value between 0 and 1
     * also make sure that in 1000 iterations, the same value is not repeated
    */
    std::size_t const LEN(1000);
    double dart;
    std::set<double> darts;
    for (std::size_t i=0; i<LEN; i++){
        dart = rootfind::pseudo_unif();
        CHECK_LE( 0.0, dart);
        CHECK_LE(dart,  1.0);
        CHECK_EQ(darts.count(dart), 0);
        darts.insert(dart);
    }
}

  TEST_CASE("rootfind: pseudo_unif_range") {
    /**
     * make sure the pseudo_unif produces a value between xmin and xmax
     * also make sure that in 1000 iterations, the same value is not repeated
    */
    std::size_t const LEN(1000);
    double xmin = 3.0, xmax = 7.0, dart;
    std::set<double> darts;
    for (std::size_t i=0; i<LEN; i++){
        dart = rootfind::pseudo_unif(xmin, xmax);
        CHECK_LE(xmin, dart);
        CHECK_LE(dart, xmax);
        CHECK_EQ(darts.count(dart), 0);
        darts.insert(dart);
    }
}

TEST_CASE_FIXTURE(RootFindLineFixture, "rootfind: find_brackets_move_line_up"){
    double xmin, xmax, xguess;

    // test situation
    //        /
    // ---[--O-x]----
    //      /
    xguess = over_root;
    rootfind::bisection_find_brackets_move(slopeup, xguess, xmin, xmax);
    CHECK_LT( xmin, xmax );
    CHECK_LT( xmin, root );
    CHECK_LT( root, xmax );
    // specific to side of zero
    CHECK_LT( xmin, xguess );
    CHECK_EQ( xguess, xmax );
    // specific to slope up
    CHECK_LT( slopeup(xmin), slopeup(xmax) );
    CHECK_LT( slopeup(xmin), 0.0 );
    CHECK_LT( 0.0, slopeup(xmax) );

    // test situation
    //        /
    // ---[x-O--]----
    //      /       
    xmin = xmax = 0.0;
    xguess = under_root;
    rootfind::bisection_find_brackets_move(slopeup, xguess, xmin, xmax);
    CHECK_LT( xmin, xmax );
    CHECK_LT( xmin, root );
    CHECK_LT( root, xmax );
    // specific to side of zero
    CHECK_EQ( xmin, xguess );
    CHECK_LT( xguess, xmax);
    // specific to slope up
    CHECK_LT( slopeup(xmin), slopeup(xmax) );
    CHECK_LT( slopeup(xmin), 0.0 );
    CHECK_LT( 0.0, slopeup(xmax));
}

TEST_CASE_FIXTURE(RootFindLineFixture, "rootfind: find_brackets_move_line_down"){
    double xmin, xmax, xguess;
    /* test situation
    //      \        
    // ---[x-O--]----
    //        \    */
    xmin = xmax = 0.0;
    xguess = under_root;
    rootfind::bisection_find_brackets_move(slopedown, xguess, xmin, xmax);
    CHECK_LT( xmin, xmax );
    CHECK_LT( xmin, root );
    CHECK_LT( root, xmax );
    // specific to side of zero
    CHECK_EQ( xmin, xguess );
    CHECK_LT(xguess, xmax);
    // specific to slope up
    CHECK_LT( slopedown(xmax), slopedown(xmin) );
    CHECK_LT( slopedown(xmax), 0.0 );
    CHECK_LT( 0.0, slopedown(xmin));

    /* test situation
    //      \        
    // ---[--O-x]----
    //        \     */      
    xmin = xmax = 0.0;
    xguess = over_root;
    rootfind::bisection_find_brackets_move(slopedown, xguess, xmin, xmax);
    CHECK_LT( xmin, xmax );
    CHECK_LT( xmin, root );
    CHECK_LT( root, xmax );
    // specific to side of zero
    CHECK_LT( xmin, xguess );
    CHECK_EQ(xguess, xmax);
    // specific to slope up
    CHECK_LT( slopedown(xmax), slopedown(xmin) );
    CHECK_LT( slopedown(xmax), 0.0 );
    CHECK_LT( 0.0, slopedown(xmin));
}

TEST_CASE_FIXTURE(RootFindLineFixture, "rootfind: find_brackets_move_sine"){
    double xmin, xmax, xguess;

    /* test situation
    #  /\        /\
    # /  \      /  \
    #/-[x-o--]-/----\- 
    #      \  /      \
    #       \/      */ 
    xguess = under_root;
    rootfind::bisection_find_brackets_move(sinusoid, xguess, xmin, xmax);
    CHECK_LT( xmin, xmax );
    CHECK_LT( xmin, root );
    CHECK_LT( root, xmax );
    // specific to side of zero
    CHECK_EQ( xmin, xguess );
    CHECK_LT( xguess, xmax );
    // specific to slope down
    CHECK_LT( sinusoid(xmax), sinusoid(xmin) );
    CHECK_LT( sinusoid(xmax), 0.0 );
    CHECK_LT( 0.0, sinusoid(xmin));

    /* test situation
    #  /\        /\
    # /  \      /  \
    #/-[--o-x]-/----\- 
    #      \  /      \
    #       \/      */ 
    xguess = over_root;
    rootfind::bisection_find_brackets_move(sinusoid, xguess, xmin, xmax);
    CHECK_LT( xmin, xmax );
    CHECK_LT( xmin, root );
    CHECK_LT( root, xmax );
    // specific to side of zero
    CHECK_EQ( xmax, xguess );
    CHECK_LT( xmin, xguess );
    // specific to slope down
    CHECK_LT( sinusoid(xmax), sinusoid(xmin) );
    CHECK_LT( sinusoid(xmax), 0.0 );
    CHECK_LT( 0.0, sinusoid(xmin));

    /* test situation
    #  /\        /\
    # /  \      /  \
    #/----\-[x-o--]-\- 
    #      \  /      \
    #       \/      */ 
    xguess = 2.*under_root;
    rootfind::bisection_find_brackets_move(sinusoid, xguess, xmin, xmax);
    CHECK_LT( xmin, xmax );
    CHECK_LT( xmin, 2.*root );
    CHECK_LT( 2.*root, xmax );
    // specific to side of zero
    CHECK_EQ( xmin, xguess );
    CHECK_LT(xguess, xmax);
    // specific to slope up
    CHECK_LT( sinusoid(xmin), sinusoid(xmax) );
    CHECK_LT( sinusoid(xmin), 0.0 );
    CHECK_LT( 0.0, sinusoid(xmax));

    /* test situation
    #  /\        /\
    # /  \      /  \
    #/----\-[--o-x]-\- 
    #      \  /      \
    #       \/      */ 
    xguess = 2.*over_root;
    rootfind::bisection_find_brackets_move(sinusoid, xguess, xmin, xmax);
    CHECK_LT( xmin, xmax );
    CHECK_LT( xmin, 2.*root );
    CHECK_LT( 2.*root, xmax );
    // specific to side of zero
    CHECK_EQ( xmax, xguess );
    CHECK_LT(xmin, xguess);
    // specific to slope up
    CHECK_LT( sinusoid(xmin), sinusoid(xmax) );
    CHECK_LT( sinusoid(xmin), 0.0 );
    CHECK_LT( 0.0, sinusoid(xmax));

    /* test situation - concave down
    #  /-\         /-\
    # /   \       /   \
    #/-[x--o-]---/-----\- 
    #       \   /       \
    #        \_/       */ 
    xguess = 1.0;
    rootfind::bisection_find_brackets_move(sinusoid, xguess, xmin, xmax);
    CHECK_LT( xmin, xmax );
    CHECK_LT( xmin, root );
    CHECK_LT( root, xmax );
    // specific to side of zero
    CHECK_EQ( xmin, xguess );
    CHECK_LT(xguess, xmax);
    // specific to concave down (same as slope down)
    CHECK_LT( sinusoid(xmax), sinusoid(xmin) );
    CHECK_LT( sinusoid(xmax), 0.0 );
    CHECK_LT( 0.0, sinusoid(xmin));

    /* test situation - concave up
    #  /-\         /_\
    # /   \       /   \
    #/-----\-[x--o-]---\- 
    #       \   /       \
    #        \_/        */ 
    xguess = 3.0;
    rootfind::bisection_find_brackets_move(sinusoid, xguess, xmin, xmax);
    CHECK_LT( xmin, xmax );
    CHECK_LT( xmin, 2.*root );
    CHECK_LT( 2.*root, xmax );
    // specific to side of zero
    CHECK_EQ( xmin, xguess );
    CHECK_LT(xguess, xmax);
    // specific to concave up (same as slpe up)
    CHECK_LT( sinusoid(xmin), sinusoid(xmax) );
    CHECK_LT( sinusoid(xmin), 0.0 );
    CHECK_LT( 0.0, sinusoid(xmax));
}

TEST_CASE_FIXTURE(RootFindLineFixture, "rootfind: find_brackets_move_failures") {
    double xmin, xmax, xguess;

    /* test situation
    #----------------             
    #----x----------- 
    #              */
    xguess = under_root;
    rootfind::bisection_find_brackets_move(flat, xguess, xmin, xmax);
    CHECK(doctest::IsNaN<double>(xmin));
    CHECK(doctest::IsNaN<double>(xmax));

    /* test situation
    #  \   /
    #   \-/          
    #----x----------- 
    #
    #              */
    xguess = root;
    rootfind::bisection_find_brackets_move(parabola, xguess, xmin, xmax);
    CHECK(doctest::IsNaN<double>(xmin));
    CHECK(doctest::IsNaN<double>(xmax));
}

TEST_CASE_FIXTURE(RootFindLineFixture, "rootfind: find_brackets_newton"){
   double xmin, xmax, xguess;

    /* test situation
    #  /\        /\
    # /  \      /  \
    #/-[x-o--]-/----\- 
    #      \  /      \
    #       \/      */ 
    xguess = under_root;
    rootfind::bisection_find_brackets_newton(sinusoid, xguess, xmin, xmax);
    CHECK_LT( xmin, xmax );
    CHECK_LT( xmin, root );
    CHECK_LT( root, xmax );
    // specific to side of zero
    CHECK_EQ( xmin, xguess );
    CHECK_LT( xguess, xmax );
    // specific to slope down
    CHECK_LT( sinusoid(xmax), sinusoid(xmin) );
    CHECK_LT( sinusoid(xmax), 0.0 );
    CHECK_LT( 0.0, sinusoid(xmin));

    /* test situation
    #  /\        /\
    # /  \      /  \
    #/-[--o-x]-/----\- 
    #      \  /      \
    #       \/      */ 
    xguess = over_root;
    rootfind::bisection_find_brackets_newton(sinusoid, xguess, xmin, xmax);
    CHECK_LT( xmin, xmax );
    CHECK_LT( xmin, root );
    CHECK_LT( root, xmax );
    // specific to side of zero
    CHECK_EQ( xmax, xguess );
    CHECK_LT( xmin, xguess );
    // specific to slope down
    CHECK_LT( sinusoid(xmax), sinusoid(xmin) );
    CHECK_LT( sinusoid(xmax), 0.0 );
    CHECK_LT( 0.0, sinusoid(xmin));

    /* test situation
    #  /\        /\
    # /  \      /  \
    #/----\-[x-o--]-\- 
    #      \  /      \
    #       \/      */ 
    xguess = 2.*under_root;
    rootfind::bisection_find_brackets_newton(sinusoid, xguess, xmin, xmax);
    CHECK_LT( xmin, xmax );
    CHECK_LT( xmin, 2.*root );
    CHECK_LT( 2.*root, xmax );
    // specific to side of zero
    CHECK_EQ( xmin, xguess );
    CHECK_LT(xguess, xmax);
    // specific to slope up
    CHECK_LT( sinusoid(xmin), sinusoid(xmax) );
    CHECK_LT( sinusoid(xmin), 0.0 );
    CHECK_LT( 0.0, sinusoid(xmax));

    /* test situation
    #  /\        /\
    # /  \      /  \
    #/----\-[--o-x]-\- 
    #      \  /      \
    #       \/      */ 
    xguess = 2.*over_root;
    rootfind::bisection_find_brackets_newton(sinusoid, xguess, xmin, xmax);
    CHECK_LT( xmin, xmax );
    CHECK_LT( xmin, 2.*root );
    CHECK_LT( 2.*root, xmax );
    // specific to side of zero
    CHECK_EQ( xmax, xguess );
    CHECK_LT(xmin, xguess);
    // specific to slope up
    CHECK_LT( sinusoid(xmin), sinusoid(xmax) );
    CHECK_LT( sinusoid(xmin), 0.0 );
    CHECK_LT( 0.0, sinusoid(xmax));

    /* test situation - concave down
    #  /-\         /-\
    # /   \       /   \
    #/-[x--o-]---/-----\- 
    #       \   /       \
    #        \_/       */ 
    xguess = 1.0;
    rootfind::bisection_find_brackets_newton(sinusoid, xguess, xmin, xmax);
    CHECK_LT( xmin, xmax );
    CHECK_LT( xmin, root );
    CHECK_LT( root, xmax );
    // specific to side of zero
    CHECK_EQ( xmin, xguess );
    CHECK_LT(xguess, xmax);
    // specific to concave down (same as slope down)
    CHECK_LT( sinusoid(xmax), sinusoid(xmin) );
    CHECK_LT( sinusoid(xmax), 0.0 );
    CHECK_LT( 0.0, sinusoid(xmin));

    /* test situation - concave up
    #  /-\         /_\
    # /   \       /   \
    #/-----\-[x--o-]---\- 
    #       \   /       \
    #        \_/        */ 
    xguess = 3.0;
    rootfind::bisection_find_brackets_newton(sinusoid, xguess, xmin, xmax);
    CHECK_LT( xmin, xmax );
    CHECK_LT( xmin, 2.*root );
    CHECK_LT( 2.*root, xmax );
    // specific to side of zero
    CHECK_EQ( xmin, xguess );
    CHECK_LT(xguess, xmax);
    // specific to concave up (same as slpe up)
    CHECK_LT( sinusoid(xmin), sinusoid(xmax) );
    CHECK_LT( sinusoid(xmin), 0.0 );
    CHECK_LT( 0.0, sinusoid(xmax));
}

TEST_CASE_FIXTURE(RootFindLineFixture, "rootfind: bisection_search"){
    double xmin, xmax, xguess;

    /* test situation
    #  /\        /\
    # /  \      /  \
    #/-[x-o--]-/----\- 
    #      \  /      \
    #       \/      */ 
    xguess = under_root;
    rootfind::bisection_find_brackets_newton(sinusoid, xguess, xmin, xmax);
    rootfind::bisection_search(sinusoid, xguess, xmin, xmax);
    CHECK_DELTA( xguess, root, 1.e-10 );

    /* test situation
    #  /\        /\
    # /  \      /  \
    #/-[--o-x]-/----\- 
    #      \  /      \
    #       \/      */ 
    xguess = over_root;
    rootfind::bisection_find_brackets_newton(sinusoid, xguess, xmin, xmax);
    rootfind::bisection_search(sinusoid, xguess, xmin, xmax);
    CHECK_DELTA( xguess, root, 1.e-10 );

    /* test situation
    #  /\        /\
    # /  \      /  \
    #/----\-[x-o--]-\- 
    #      \  /      \
    #       \/      */ 
    xguess = 2.*under_root;
    rootfind::bisection_find_brackets_newton(sinusoid, xguess, xmin, xmax);
    rootfind::bisection_search(sinusoid, xguess, xmin, xmax);
    CHECK_DELTA( xguess, 2.0*root, 1.e-10 );

    /* test situation
    #  /\        /\
    # /  \      /  \
    #/----\-[--o-x]-\- 
    #      \  /      \
    #       \/      */ 
    xguess = 2.*over_root;
    rootfind::bisection_find_brackets_newton(sinusoid, xguess, xmin, xmax);
    rootfind::bisection_search(sinusoid, xguess, xmin, xmax);
    CHECK_DELTA( xguess, 2.0*root, 1.e-10 );

    /* test situation - concave down
    #  /-\         /-\
    # /   \       /   \
    #/-[x--o-]---/-----\- 
    #       \   /       \
    #        \_/       */ 
    xguess = 1.0;
    rootfind::bisection_find_brackets_newton(sinusoid, xguess, xmin, xmax);
    rootfind::bisection_search(sinusoid, xguess, xmin, xmax);
    CHECK_DELTA( xguess, root, 1.e-10 );

    /* test situation - concave up
    #  /-\         /_\
    # /   \       /   \
    #/-----\-[x--o-]---\- 
    #       \   /       \
    #        \_/        */ 
    xguess = 3.0;
    rootfind::bisection_find_brackets_newton(sinusoid, xguess, xmin, xmax);
    rootfind::bisection_search(sinusoid, xguess, xmin, xmax);
    CHECK_DELTA( xguess, 2.*root, 1.e-10 );
}

TEST_CASE_FIXTURE(RootFindLineFixture, "rootfind: newton_search_1d"){
    double xguess, tol = 1.e-10;
    
    // happy path -- the guess is the root
    xguess = root;
    rootfind::newton_search<double>(sinusoid, xguess, 0.1, tol, 0);
    CHECK_EQ( xguess, root );

    // the root is nearby
    xguess = under_root;
    rootfind::newton_search<double>(sinusoid, xguess, 0.1, tol, 0);
    CHECK_DELTA( xguess, root, tol );

    xguess = over_root;
    rootfind::newton_search<double>(sinusoid, xguess, 0.1, tol, 0);
    CHECK_DELTA( xguess, root, tol );
    // if only allowed on step, will not find root within tolerance
    xguess = under_root;
    rootfind::newton_search<double>(sinusoid, xguess, 0.1, tol, 1);
    CHECK_LT( tol, fabs(xguess-root) );

    // finds root before ten iterations
    xguess = under_root;
    rootfind::newton_search<double>(sinusoid, xguess, 0.1, tol, 10);
    CHECK_DELTA( xguess, root, tol );
}

TEST_CASE_FIXTURE(RootFindLineFixture, "rootfind: newton_search_target_1d"){
    double xguess, tol = 1.e-10;
    double xtarget = 1./3., ytarget=0.5;
    
    // happy path -- the guess is the target
    xguess = xtarget;
    rootfind::newton_search<double>(sinusoid, ytarget, xguess, 0.1, tol, 0);
    CHECK_EQ( xguess, xtarget );
    CHECK_DELTA( sinusoid(xguess), ytarget, tol );

    // the target is nearby
    xguess = 0.1;
    rootfind::newton_search<double>(sinusoid, ytarget, xguess, 0.1, tol, 0);
    CHECK_DELTA( xguess, xtarget, tol );

    // finds target before ten iterations
    xguess = 0.1;
    rootfind::newton_search<double>(sinusoid, ytarget, xguess, 0.1, tol, 10);
    CHECK_DELTA( xguess, xtarget, tol );
}

TEST_CASE_FIXTURE(RootFindLineFixture, "rootfind: newton_search_1d_complex"){
    typedef std::complex<double> C;
    C croot {this->root, 0.5*this->root};
    std::function<C(C)> up_complex = [croot](C x)->C{return (x-croot);};
    C xguess;
    double tol = 1.e-10;
    
    // happy path -- the guess is the root
    xguess = root;
    rootfind::newton_search<C>(up_complex, xguess, C{0.1,0.2}, tol, 0);
    CHECK_LT( std::abs(xguess-croot), 1.e-15 );

    // the root is nearby
    xguess = under_root;
    rootfind::newton_search<C>(up_complex, xguess, C{0.1,0.2}, tol, 0);
    CHECK_LT( std::abs(xguess-croot), tol );

    xguess = over_root;
    rootfind::newton_search<C>(up_complex, xguess, C{0.1,0.2}, tol, 0);
    CHECK_LT( std::abs(xguess-croot), tol );

    // if only allowed one step, will not find root within tolerance
    xguess = under_root;
    rootfind::newton_search<C>(up_complex, xguess, C{0.1,0.2}, tol, 1);
    CHECK_LT( tol, std::abs(xguess-croot) );

    // finds root before ten iterations
    xguess = under_root;
    rootfind::newton_search<C>(up_complex, xguess, C{0.1,0.2}, tol, 10);
    CHECK_LT( std::abs(xguess-croot), tol );

    // test finding a target 
    C target {1.5*this->root, 3.0*this->root};
    C xroot  = croot + target;
    xguess = under_root;
    rootfind::newton_search<C>(up_complex, target, xguess, C{0.2,0.1}, tol, 10);
    CHECK_LT( std::abs(xguess-xroot), tol );
}

TEST_CASE_FIXTURE(RootFind2DFixture, "rootfind: newton_search_2d") {
    double xguess[2];
    double dx[2] = {0.1,0.2};
    double tol = 1.0e-10;

    // happy path -- the guess is the root
    xguess[0] = xroot[0];
    xguess[1] = xroot[1];
    rootfind::newton_search(upup, xguess, dx, tol, 0);
    CHECK_DELTA( xguess[0], xroot[0], 1.e-15 );
    CHECK_DELTA( xguess[1], xroot[1], 1.e-15 );
    xguess[0] = xroot[0];
    xguess[1] = xroot[1];
    rootfind::newton_search(updown, xguess, dx, tol, 0);
    CHECK_DELTA( xguess[0], xroot[0], 1.e-15 );
    CHECK_DELTA( xguess[1], xroot[1], 1.e-15 );
    xguess[0] = xroot[0];
    xguess[1] = xroot[1];
    rootfind::newton_search(downup, xguess, dx, tol, 0);
    CHECK_DELTA( xguess[0], xroot[0], 1.e-15 );
    CHECK_DELTA( xguess[1], xroot[1], 1.e-15 );
    xguess[0] = xroot[0];
    xguess[1] = xroot[1];
    rootfind::newton_search(downdown, xguess, dx, tol, 0);
    CHECK_DELTA( xguess[0], xroot[0], 1.e-15 );
    CHECK_DELTA( xguess[1], xroot[1], 1.e-15 );

    // the root can be found
    xguess[0] = 3.*xroot[0];
    xguess[1] = 5.*xroot[1];
    rootfind::newton_search(upup, xguess, dx, tol, 10);
    CHECK_DELTA( xguess[0], xroot[0], tol );
    CHECK_DELTA( xguess[1], xroot[1], tol );
    xguess[0] = 3.*xroot[0];
    xguess[1] = 5.*xroot[1];
    rootfind::newton_search(updown, xguess, dx, tol, 10);
    CHECK_DELTA( xguess[0], xroot[0], tol );
    CHECK_DELTA( xguess[1], xroot[1], tol );
    xguess[0] = 3.*xroot[0];
    xguess[1] = 5.*xroot[1];
    rootfind::newton_search(downup, xguess, dx, tol, 10);
    CHECK_DELTA( xguess[0], xroot[0], tol );
    CHECK_DELTA( xguess[1], xroot[1], tol );
    xguess[0] = 3.*xroot[0];
    xguess[1] = 5.*xroot[1];
    rootfind::newton_search(downdown, xguess, dx, tol, 10);
    CHECK_DELTA( xguess[0], xroot[0], tol );
    CHECK_DELTA( xguess[1], xroot[1], tol );
}

} // end TEST_SUITE
