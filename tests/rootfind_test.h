#include "../lib/rootfind.h"
#include "../src/constants.h"
#include <cxxtest/TestSuite.h>
#include <random>
#include <set>

class RootFindTest : public CxxTest::TestSuite {
public:

    double root, over_root, under_root ;
    typedef std::function<double(double)> Line;
    Line slopeup, slopedown, sinusoid, flat, parabola;
    
    double xroot[2];
    typedef std::function<void(double[2],double[2])> Field2D;
    Field2D updown, downup, upup, downdown;

RootFindTest(){
    root = 2.0;
    over_root = 2.4;
    under_root = 1.6;
    slopeup = [this](double x)->double { return (x-this->root); };
    slopedown = [this](double x)->double { return (this->root-x); };
    sinusoid = [this](double x)->double { return sin( x * m_pi/this->root ); };
    flat = [](double x)->double {return 1.0; };
    parabola = [this](double x)->double {return (x-this->root)*(x-this->root) + 1.0;};

    xroot[0] = 1.5;
    xroot[1] = 3.7;
    upup = [this](double f[2],double x[2]){
        f[0] = (x[0] - this->xroot[0]);
        f[1] = (x[1] - this->xroot[1]);
    };
    updown = [this](double f[2],double x[2]){
        f[0] = (x[0] - this->xroot[0]);
        f[1] =-(x[1] - this->xroot[1]);
    };
    downup = [this](double f[2],double x[2]){
        f[0] =-(x[0] - this->xroot[0]);
        f[1] = (x[1] - this->xroot[1]);
    };
    downdown = [this](double f[2],double x[2]){
        f[0] =-(x[0] - this->xroot[0]);
        f[1] =-(x[1] - this->xroot[1]);
    };
}

void test_pseudo_unif(){
    /** 
     * make sure the pseudo_unif produces a value between 0 and 1
     * also make sure that in 1000 iterations, the same value is not repeated
    */
    std::size_t const LEN(1000);
    double dart;
    std::set<double> darts;
    for (std::size_t i=0; i<LEN; i++){
        dart = rootfind::pseudo_unif();
        TS_ASSERT_LESS_THAN_EQUALS( 0.0, dart);
        TS_ASSERT_LESS_THAN_EQUALS(dart,  1.0);
        TS_ASSERT_EQUALS(darts.count(dart), 0);
        darts.insert(dart);
    }
}

void test_pseudo_unif_range(){
    /**
     * make sure the pseudo_unif produces a value between xmin and xmax
     * also make sure that in 1000 iterations, the same value is not repeated
    */
    std::size_t const LEN(1000);
    double xmin = 3.0, xmax = 7.0, dart;
    std::set<double> darts;
    for (std::size_t i=0; i<LEN; i++){
        dart = rootfind::pseudo_unif(xmin, xmax);
        TS_ASSERT_LESS_THAN_EQUALS(xmin, dart);
        TS_ASSERT_LESS_THAN_EQUALS(dart, xmax);
        TS_ASSERT_EQUALS(darts.count(dart), 0);
        darts.insert(dart);
    }
}

void test_find_brackets_move_lines(){
    printf("\nROOTFINDING TESTS");

    double xmin, xmax, xguess;

    // test situation
    //        /
    // ---[--O-x]----
    //      /
    xguess = over_root;
    rootfind::bisection_find_brackets_move(slopeup, xguess, xmin, xmax);
    TS_ASSERT_LESS_THAN( xmin, xmax );
    TS_ASSERT_LESS_THAN( xmin, root );
    TS_ASSERT_LESS_THAN( root, xmax );
    // specific to side of zero
    TS_ASSERT_LESS_THAN( xmin, xguess );
    TS_ASSERT_EQUALS(xguess, xmax);
    // specific to slope up
    TS_ASSERT_LESS_THAN( slopeup(xmin), slopeup(xmax) );
    TS_ASSERT_LESS_THAN( slopeup(xmin), 0.0 );
    TS_ASSERT_LESS_THAN( 0.0, slopeup(xmax));

    // test situation
    //        /
    // ---[x-O--]----
    //      /       
    xmin = xmax = 0.0;
    xguess = under_root;
    rootfind::bisection_find_brackets_move(slopeup, xguess, xmin, xmax);
    TS_ASSERT_LESS_THAN( xmin, xmax );
    TS_ASSERT_LESS_THAN( xmin, root );
    TS_ASSERT_LESS_THAN( root, xmax );
    // specific to side of zero
    TS_ASSERT_EQUALS( xmin, xguess );
    TS_ASSERT_LESS_THAN(xguess, xmax);
    // specific to slope up
    TS_ASSERT_LESS_THAN( slopeup(xmin), slopeup(xmax) );
    TS_ASSERT_LESS_THAN( slopeup(xmin), 0.0 );
    TS_ASSERT_LESS_THAN( 0.0, slopeup(xmax));

    /* test situation
    //      \        
    // ---[x-O--]----
    //        \    */
    xmin = xmax = 0.0;
    xguess = under_root;
    rootfind::bisection_find_brackets_move(slopedown, xguess, xmin, xmax);
    TS_ASSERT_LESS_THAN( xmin, xmax );
    TS_ASSERT_LESS_THAN( xmin, root );
    TS_ASSERT_LESS_THAN( root, xmax );
    // specific to side of zero
    TS_ASSERT_EQUALS( xmin, xguess );
    TS_ASSERT_LESS_THAN(xguess, xmax);
    // specific to slope up
    TS_ASSERT_LESS_THAN( slopedown(xmax), slopedown(xmin) );
    TS_ASSERT_LESS_THAN( slopedown(xmax), 0.0 );
    TS_ASSERT_LESS_THAN( 0.0, slopedown(xmin));

    /* test situation
    //      \        
    // ---[--O-x]----
    //        \     */      
    xmin = xmax = 0.0;
    xguess = over_root;
    rootfind::bisection_find_brackets_move(slopedown, xguess, xmin, xmax);
    TS_ASSERT_LESS_THAN( xmin, xmax );
    TS_ASSERT_LESS_THAN( xmin, root );
    TS_ASSERT_LESS_THAN( root, xmax );
    // specific to side of zero
    TS_ASSERT_LESS_THAN( xmin, xguess );
    TS_ASSERT_EQUALS(xguess, xmax);
    // specific to slope up
    TS_ASSERT_LESS_THAN( slopedown(xmax), slopedown(xmin) );
    TS_ASSERT_LESS_THAN( slopedown(xmax), 0.0 );
    TS_ASSERT_LESS_THAN( 0.0, slopedown(xmin));
}

void test_find_brackets_move_sine(){
    double xmin, xmax, xguess;

    /* test situation
    #  /\        /\
    # /  \      /  \
    #/-[x-o--]-/----\- 
    #      \  /      \
    #       \/      */ 
    xguess = under_root;
    rootfind::bisection_find_brackets_move(sinusoid, xguess, xmin, xmax);
    TS_ASSERT_LESS_THAN( xmin, xmax );
    TS_ASSERT_LESS_THAN( xmin, root );
    TS_ASSERT_LESS_THAN( root, xmax );
    // specific to side of zero
    TS_ASSERT_EQUALS( xmin, xguess );
    TS_ASSERT_LESS_THAN( xguess, xmax );
    // specific to slope down
    TS_ASSERT_LESS_THAN( sinusoid(xmax), sinusoid(xmin) );
    TS_ASSERT_LESS_THAN( sinusoid(xmax), 0.0 );
    TS_ASSERT_LESS_THAN( 0.0, sinusoid(xmin));

    /* test situation
    #  /\        /\
    # /  \      /  \
    #/-[--o-x]-/----\- 
    #      \  /      \
    #       \/      */ 
    xguess = over_root;
    rootfind::bisection_find_brackets_move(sinusoid, xguess, xmin, xmax);
    TS_ASSERT_LESS_THAN( xmin, xmax );
    TS_ASSERT_LESS_THAN( xmin, root );
    TS_ASSERT_LESS_THAN( root, xmax );
    // specific to side of zero
    TS_ASSERT_EQUALS( xmax, xguess );
    TS_ASSERT_LESS_THAN( xmin, xguess );
    // specific to slope down
    TS_ASSERT_LESS_THAN( sinusoid(xmax), sinusoid(xmin) );
    TS_ASSERT_LESS_THAN( sinusoid(xmax), 0.0 );
    TS_ASSERT_LESS_THAN( 0.0, sinusoid(xmin));

    /* test situation
    #  /\        /\
    # /  \      /  \
    #/----\-[x-o--]-\- 
    #      \  /      \
    #       \/      */ 
    xguess = 2.*under_root;
    rootfind::bisection_find_brackets_move(sinusoid, xguess, xmin, xmax);
    TS_ASSERT_LESS_THAN( xmin, xmax );
    TS_ASSERT_LESS_THAN( xmin, 2.*root );
    TS_ASSERT_LESS_THAN( 2.*root, xmax );
    // specific to side of zero
    TS_ASSERT_EQUALS( xmin, xguess );
    TS_ASSERT_LESS_THAN(xguess, xmax);
    // specific to slope up
    TS_ASSERT_LESS_THAN( sinusoid(xmin), sinusoid(xmax) );
    TS_ASSERT_LESS_THAN( sinusoid(xmin), 0.0 );
    TS_ASSERT_LESS_THAN( 0.0, sinusoid(xmax));

    /* test situation
    #  /\        /\
    # /  \      /  \
    #/----\-[--o-x]-\- 
    #      \  /      \
    #       \/      */ 
    xguess = 2.*over_root;
    rootfind::bisection_find_brackets_move(sinusoid, xguess, xmin, xmax);
    TS_ASSERT_LESS_THAN( xmin, xmax );
    TS_ASSERT_LESS_THAN( xmin, 2.*root );
    TS_ASSERT_LESS_THAN( 2.*root, xmax );
    // specific to side of zero
    TS_ASSERT_EQUALS( xmax, xguess );
    TS_ASSERT_LESS_THAN(xmin, xguess);
    // specific to slope up
    TS_ASSERT_LESS_THAN( sinusoid(xmin), sinusoid(xmax) );
    TS_ASSERT_LESS_THAN( sinusoid(xmin), 0.0 );
    TS_ASSERT_LESS_THAN( 0.0, sinusoid(xmax));

    /* test situation - concave down
    #  /-\         /-\
    # /   \       /   \
    #/-[x--o-]---/-----\- 
    #       \   /       \
    #        \_/       */ 
    xguess = 1.0;
    rootfind::bisection_find_brackets_move(sinusoid, xguess, xmin, xmax);
    TS_ASSERT_LESS_THAN( xmin, xmax );
    TS_ASSERT_LESS_THAN( xmin, root );
    TS_ASSERT_LESS_THAN( root, xmax );
    // specific to side of zero
    TS_ASSERT_EQUALS( xmin, xguess );
    TS_ASSERT_LESS_THAN(xguess, xmax);
    // specific to concave down (same as slope down)
    TS_ASSERT_LESS_THAN( sinusoid(xmax), sinusoid(xmin) );
    TS_ASSERT_LESS_THAN( sinusoid(xmax), 0.0 );
    TS_ASSERT_LESS_THAN( 0.0, sinusoid(xmin));

    /* test situation - concave up
    #  /-\         /_\
    # /   \       /   \
    #/-----\-[x--o-]---\- 
    #       \   /       \
    #        \_/        */ 
    xguess = 3.0;
    rootfind::bisection_find_brackets_move(sinusoid, xguess, xmin, xmax);
    TS_ASSERT_LESS_THAN( xmin, xmax );
    TS_ASSERT_LESS_THAN( xmin, 2.*root );
    TS_ASSERT_LESS_THAN( 2.*root, xmax );
    // specific to side of zero
    TS_ASSERT_EQUALS( xmin, xguess );
    TS_ASSERT_LESS_THAN(xguess, xmax);
    // specific to concave up (same as slpe up)
    TS_ASSERT_LESS_THAN( sinusoid(xmin), sinusoid(xmax) );
    TS_ASSERT_LESS_THAN( sinusoid(xmin), 0.0 );
    TS_ASSERT_LESS_THAN( 0.0, sinusoid(xmax));
}

void test_find_brackets_move_failures(){
    double xmin, xmax, xguess;

    /* test situation
    #----------------             
    #----x----------- 
    #              */
    xguess = under_root;
    rootfind::bisection_find_brackets_move(flat, xguess, xmin, xmax);
    TS_ASSERT(std::isnan(xmin));
    TS_ASSERT(std::isnan(xmax));

    /* test situation
    #  \   /
    #   \-/          
    #----x----------- 
    #
    #              */
    xguess = root;
    rootfind::bisection_find_brackets_move(parabola, xguess, xmin, xmax);
    TS_ASSERT(std::isnan(xmin));
    TS_ASSERT(std::isnan(xmax));
}

void test_find_brackets_newton(){
   double xmin, xmax, xguess;

    /* test situation
    #  /\        /\
    # /  \      /  \
    #/-[x-o--]-/----\- 
    #      \  /      \
    #       \/      */ 
    xguess = under_root;
    rootfind::bisection_find_brackets_newton(sinusoid, xguess, xmin, xmax);
    TS_ASSERT_LESS_THAN( xmin, xmax );
    TS_ASSERT_LESS_THAN( xmin, root );
    TS_ASSERT_LESS_THAN( root, xmax );
    // specific to side of zero
    TS_ASSERT_EQUALS( xmin, xguess );
    TS_ASSERT_LESS_THAN( xguess, xmax );
    // specific to slope down
    TS_ASSERT_LESS_THAN( sinusoid(xmax), sinusoid(xmin) );
    TS_ASSERT_LESS_THAN( sinusoid(xmax), 0.0 );
    TS_ASSERT_LESS_THAN( 0.0, sinusoid(xmin));

    /* test situation
    #  /\        /\
    # /  \      /  \
    #/-[--o-x]-/----\- 
    #      \  /      \
    #       \/      */ 
    xguess = over_root;
    rootfind::bisection_find_brackets_newton(sinusoid, xguess, xmin, xmax);
    TS_ASSERT_LESS_THAN( xmin, xmax );
    TS_ASSERT_LESS_THAN( xmin, root );
    TS_ASSERT_LESS_THAN( root, xmax );
    // specific to side of zero
    TS_ASSERT_EQUALS( xmax, xguess );
    TS_ASSERT_LESS_THAN( xmin, xguess );
    // specific to slope down
    TS_ASSERT_LESS_THAN( sinusoid(xmax), sinusoid(xmin) );
    TS_ASSERT_LESS_THAN( sinusoid(xmax), 0.0 );
    TS_ASSERT_LESS_THAN( 0.0, sinusoid(xmin));

    /* test situation
    #  /\        /\
    # /  \      /  \
    #/----\-[x-o--]-\- 
    #      \  /      \
    #       \/      */ 
    xguess = 2.*under_root;
    rootfind::bisection_find_brackets_newton(sinusoid, xguess, xmin, xmax);
    TS_ASSERT_LESS_THAN( xmin, xmax );
    TS_ASSERT_LESS_THAN( xmin, 2.*root );
    TS_ASSERT_LESS_THAN( 2.*root, xmax );
    // specific to side of zero
    TS_ASSERT_EQUALS( xmin, xguess );
    TS_ASSERT_LESS_THAN(xguess, xmax);
    // specific to slope up
    TS_ASSERT_LESS_THAN( sinusoid(xmin), sinusoid(xmax) );
    TS_ASSERT_LESS_THAN( sinusoid(xmin), 0.0 );
    TS_ASSERT_LESS_THAN( 0.0, sinusoid(xmax));

    /* test situation
    #  /\        /\
    # /  \      /  \
    #/----\-[--o-x]-\- 
    #      \  /      \
    #       \/      */ 
    xguess = 2.*over_root;
    rootfind::bisection_find_brackets_newton(sinusoid, xguess, xmin, xmax);
    TS_ASSERT_LESS_THAN( xmin, xmax );
    TS_ASSERT_LESS_THAN( xmin, 2.*root );
    TS_ASSERT_LESS_THAN( 2.*root, xmax );
    // specific to side of zero
    TS_ASSERT_EQUALS( xmax, xguess );
    TS_ASSERT_LESS_THAN(xmin, xguess);
    // specific to slope up
    TS_ASSERT_LESS_THAN( sinusoid(xmin), sinusoid(xmax) );
    TS_ASSERT_LESS_THAN( sinusoid(xmin), 0.0 );
    TS_ASSERT_LESS_THAN( 0.0, sinusoid(xmax));

    /* test situation - concave down
    #  /-\         /-\
    # /   \       /   \
    #/-[x--o-]---/-----\- 
    #       \   /       \
    #        \_/       */ 
    xguess = 1.0;
    rootfind::bisection_find_brackets_newton(sinusoid, xguess, xmin, xmax);
    TS_ASSERT_LESS_THAN( xmin, xmax );
    TS_ASSERT_LESS_THAN( xmin, root );
    TS_ASSERT_LESS_THAN( root, xmax );
    // specific to side of zero
    TS_ASSERT_EQUALS( xmin, xguess );
    TS_ASSERT_LESS_THAN(xguess, xmax);
    // specific to concave down (same as slope down)
    TS_ASSERT_LESS_THAN( sinusoid(xmax), sinusoid(xmin) );
    TS_ASSERT_LESS_THAN( sinusoid(xmax), 0.0 );
    TS_ASSERT_LESS_THAN( 0.0, sinusoid(xmin));

    /* test situation - concave up
    #  /-\         /_\
    # /   \       /   \
    #/-----\-[x--o-]---\- 
    #       \   /       \
    #        \_/        */ 
    xguess = 3.0;
    rootfind::bisection_find_brackets_newton(sinusoid, xguess, xmin, xmax);
    TS_ASSERT_LESS_THAN( xmin, xmax );
    TS_ASSERT_LESS_THAN( xmin, 2.*root );
    TS_ASSERT_LESS_THAN( 2.*root, xmax );
    // specific to side of zero
    TS_ASSERT_EQUALS( xmin, xguess );
    TS_ASSERT_LESS_THAN(xguess, xmax);
    // specific to concave up (same as slpe up)
    TS_ASSERT_LESS_THAN( sinusoid(xmin), sinusoid(xmax) );
    TS_ASSERT_LESS_THAN( sinusoid(xmin), 0.0 );
    TS_ASSERT_LESS_THAN( 0.0, sinusoid(xmax));
}

void test_bisection_search(){
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
    TS_ASSERT_DELTA( xguess, root, 1.e-10 );

    /* test situation
    #  /\        /\
    # /  \      /  \
    #/-[--o-x]-/----\- 
    #      \  /      \
    #       \/      */ 
    xguess = over_root;
    rootfind::bisection_find_brackets_newton(sinusoid, xguess, xmin, xmax);
    rootfind::bisection_search(sinusoid, xguess, xmin, xmax);
    TS_ASSERT_DELTA( xguess, root, 1.e-10 );

    /* test situation
    #  /\        /\
    # /  \      /  \
    #/----\-[x-o--]-\- 
    #      \  /      \
    #       \/      */ 
    xguess = 2.*under_root;
    rootfind::bisection_find_brackets_newton(sinusoid, xguess, xmin, xmax);
    rootfind::bisection_search(sinusoid, xguess, xmin, xmax);
    TS_ASSERT_DELTA( xguess, 2.0*root, 1.e-10 );

    /* test situation
    #  /\        /\
    # /  \      /  \
    #/----\-[--o-x]-\- 
    #      \  /      \
    #       \/      */ 
    xguess = 2.*over_root;
    rootfind::bisection_find_brackets_newton(sinusoid, xguess, xmin, xmax);
    rootfind::bisection_search(sinusoid, xguess, xmin, xmax);
    TS_ASSERT_DELTA( xguess, 2.0*root, 1.e-10 );

    /* test situation - concave down
    #  /-\         /-\
    # /   \       /   \
    #/-[x--o-]---/-----\- 
    #       \   /       \
    #        \_/       */ 
    xguess = 1.0;
    rootfind::bisection_find_brackets_newton(sinusoid, xguess, xmin, xmax);
    rootfind::bisection_search(sinusoid, xguess, xmin, xmax);
    TS_ASSERT_DELTA( xguess, root, 1.e-10 );

    /* test situation - concave up
    #  /-\         /_\
    # /   \       /   \
    #/-----\-[x--o-]---\- 
    #       \   /       \
    #        \_/        */ 
    xguess = 3.0;
    rootfind::bisection_find_brackets_newton(sinusoid, xguess, xmin, xmax);
    rootfind::bisection_search(sinusoid, xguess, xmin, xmax);
    TS_ASSERT_DELTA( xguess, 2.*root, 1.e-10 );
}

void test_newton_search_1d(){
    double xguess, tol = 1.e-10;
    
    // happy path -- the guess is the root
    xguess = root;
    rootfind::newton_search(sinusoid, xguess, 0.1, tol, 0);
    TS_ASSERT_EQUALS( xguess, root );

    // the root is nearby
    xguess = under_root;
    rootfind::newton_search(sinusoid, xguess, 0.1, tol, 0);
    TS_ASSERT_DELTA( xguess, root, tol );

    xguess = over_root;
    rootfind::newton_search(sinusoid, xguess, 0.1, tol, 0);
    TS_ASSERT_DELTA( xguess, root, tol );

    // if only allowed on step, will not find root within tolerance
    xguess = under_root;
    rootfind::newton_search(sinusoid, xguess, 0.1, tol, 1);
    TS_ASSERT_LESS_THAN( tol, fabs(xguess-root) );

    // finds root before ten iterations
    xguess = under_root;
    rootfind::newton_search(sinusoid, xguess, 0.1, tol, 10);
    TS_ASSERT_DELTA( xguess, root, tol );
}

void test_newton_search_target_1d(){
    double xguess, tol = 1.e-10;
    double xtarget = 1./3., ytarget=0.5;
    
    // happy path -- the guess is the target
    xguess = xtarget;
    rootfind::newton_search(sinusoid, ytarget, xguess, 0.1, tol, 0);
    TS_ASSERT_EQUALS( xguess, xtarget );
    TS_ASSERT_DELTA( sinusoid(xguess), ytarget, tol );

    // the target is nearby
    xguess = 0.1;
    rootfind::newton_search(sinusoid, ytarget, xguess, 0.1, tol, 0);
    TS_ASSERT_DELTA( xguess, xtarget, tol );

    // finds target before ten iterations
    xguess = 0.1;
    rootfind::newton_search(sinusoid, ytarget, xguess, 0.1, tol, 10);
    TS_ASSERT_DELTA( xguess, xtarget, tol );
}

void test_newton_search_1d_complex(){
    typedef std::complex<double> C;
    C croot {this->root, 0.5*this->root};
    std::function<C(C)> up_complex = [croot](C x)->C{return (x-croot);};
    C xguess;
    double tol = 1.e-10;
    
    // happy path -- the guess is the root
    xguess = root;
    rootfind::newton_search(up_complex, xguess, C{0.1,0.2}, tol, 0);
    TS_ASSERT_LESS_THAN( std::abs(xguess-croot), 1.e-15 );

    // the root is nearby
    xguess = under_root;
    rootfind::newton_search(up_complex, xguess, C{0.1,0.2}, tol, 0);
    TS_ASSERT_LESS_THAN( std::abs(xguess-croot), tol );

    xguess = over_root;
    rootfind::newton_search(up_complex, xguess, C{0.1,0.2}, tol, 0);
    TS_ASSERT_LESS_THAN( std::abs(xguess-croot), tol );

    // if only allowed one step, will not find root within tolerance
    xguess = under_root;
    rootfind::newton_search(up_complex, xguess, C{0.1,0.2}, tol, 1);
    TS_ASSERT_LESS_THAN( tol, std::abs(xguess-croot) );

    // finds root before ten iterations
    xguess = under_root;
    rootfind::newton_search(up_complex, xguess, C{0.1,0.2}, tol, 10);
    TS_ASSERT_LESS_THAN( std::abs(xguess-croot), tol );


    // test finding a target 
    C target {1.5*this->root, 3.0*this->root};
    C xroot  = croot + target;
    xguess = under_root;
    rootfind::newton_search(up_complex, target, xguess, C{0.2,0.1}, tol, 10);
    TS_ASSERT_LESS_THAN( std::abs(xguess-xroot), tol );
}

void test_newton_search_2d(){
    double xguess[2];
    double dx[2] = {0.1,0.2};
    double tol = 1.0e-10;

    // happy path -- the guess is the root
    xguess[0] = xroot[0];
    xguess[1] = xroot[1];
    rootfind::newton_search(upup, xguess, dx, tol, 0);
    TS_ASSERT_DELTA( xguess[0], xroot[0], 1.e-15 );
    TS_ASSERT_DELTA( xguess[1], xroot[1], 1.e-15 );
    xguess[0] = xroot[0];
    xguess[1] = xroot[1];
    rootfind::newton_search(updown, xguess, dx, tol, 0);
    TS_ASSERT_DELTA( xguess[0], xroot[0], 1.e-15 );
    TS_ASSERT_DELTA( xguess[1], xroot[1], 1.e-15 );
    xguess[0] = xroot[0];
    xguess[1] = xroot[1];
    rootfind::newton_search(downup, xguess, dx, tol, 0);
    TS_ASSERT_DELTA( xguess[0], xroot[0], 1.e-15 );
    TS_ASSERT_DELTA( xguess[1], xroot[1], 1.e-15 );
    xguess[0] = xroot[0];
    xguess[1] = xroot[1];
    rootfind::newton_search(downdown, xguess, dx, tol, 0);
    TS_ASSERT_DELTA( xguess[0], xroot[0], 1.e-15 );
    TS_ASSERT_DELTA( xguess[1], xroot[1], 1.e-15 );

    // the root can be found
    xguess[0] = 3.*xroot[0];
    xguess[1] = 5.*xroot[1];
    rootfind::newton_search(upup, xguess, dx, tol, 10);
    TS_ASSERT_DELTA( xguess[0], xroot[0], tol );
    TS_ASSERT_DELTA( xguess[1], xroot[1], tol );
    xguess[0] = 3.*xroot[0];
    xguess[1] = 5.*xroot[1];
    rootfind::newton_search(updown, xguess, dx, tol, 10);
    TS_ASSERT_DELTA( xguess[0], xroot[0], tol );
    TS_ASSERT_DELTA( xguess[1], xroot[1], tol );
    xguess[0] = 3.*xroot[0];
    xguess[1] = 5.*xroot[1];
    rootfind::newton_search(downup, xguess, dx, tol, 10);
    TS_ASSERT_DELTA( xguess[0], xroot[0], tol );
    TS_ASSERT_DELTA( xguess[1], xroot[1], tol );
    xguess[0] = 3.*xroot[0];
    xguess[1] = 5.*xroot[1];
    rootfind::newton_search(downdown, xguess, dx, tol, 10);
    TS_ASSERT_DELTA( xguess[0], xroot[0], tol );
    TS_ASSERT_DELTA( xguess[1], xroot[1], tol );
}

};

