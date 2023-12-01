#include "../src/MODES/ModeDriver.h"
#include "../src/ThrainMode.h"
#include <cxxtest/TestSuite.h>

// TODO: check STEP 3c -- can't find max bracket
// TODO: check failure if not using Mode
// TODO: check failure if not using ModeDriver
// TODO: check g-modes (k<0)

const double ZEROW2 = 0.16;

class ModeBaseTest : public CxxTest::TestSuite {
public:

    // test that result fails if given bad polytrope index
    void test_compare_JCD_bad_n(){
        TS_ASSERT(std::isnan(mode::compare_JCD(1.2, 1, 1, 0.01)));
    }

    // test that result fails if given bad L or K for lookup
    void test_compare_JCD_bad_lk(){
        // fail L too low
        TS_ASSERT(std::isnan(mode::compare_JCD(1.5, 0, 1, 0.01)));
        // fail L too high
        TS_ASSERT(std::isnan(mode::compare_JCD(1.5, 4, 1, 0.01)));
        // fail K too low
        TS_ASSERT(std::isnan(mode::compare_JCD(1.5, 2, 0, 0.01)));
        // fail K too high
        TS_ASSERT(std::isnan(mode::compare_JCD(1.5, 1, 36, 0.01)));
    }

    void test_compare_JCD(){
        for(std::size_t L=1; L<=3; L++){
            for(std::size_t K=1; K<36; K++){
                TS_ASSERT_EQUALS(0.0, mode::compare_JCD(1.5, L, K, mode::JCD1_5[L-1][K-1]/nug));
                TS_ASSERT_EQUALS(0.0, mode::compare_JCD(3.0, L, K, mode::JCD3_0[L-1][K-1]/nug));
                TS_ASSERT_EQUALS(0.0, mode::compare_JCD(4.0, L, K, mode::JCD4_0[L-1][K-1]/nug));
            }
        }
    } 

    void test_compare_Pekeris_k0(){
        // first trst special case for K=0
        int K = 0;
        for(int L=2; L<10; L++){
            double exp = sqrt(double(2*L*(L-1))/double(2*L+1));
            double res = mode::compare_Pekeris(exp, L, K, 0.0);
            TS_ASSERT_DELTA(0.0, res, 1.0e-12);
        }
        
        // now test using Gam1=0, so that dnl = -2.0
        double Gam1 = 0.0;
        for(int L=1; L<5; L++){
            for(int K=1; K<10; K++){
                double exp = sqrt(sqrt(4.0 + double(L*L+L)) - 2.0);
                double res = mode::compare_Pekeris(exp, L, K, Gam1);
                TS_ASSERT_DELTA(0.0, res, 1.0e-12);
            }
        }

        // now test a more realistic value of Gam1
        Gam1 = 1.66667;
        for(int L=2; L<5; L++){
            for(int K=1; K<10; K++){
                double dnl = Gam1*double(K)*(double(K+L)+0.5)-2.0;
                double exp = sqrt(sqrt(dnl*dnl + double(L*L+L)) + dnl);
                double res = mode::compare_Pekeris(exp, L, K, Gam1);
                TS_ASSERT_DELTA(0.0, res, 1.0e-12);
            }
        }
    }

    /* test finding min/max from within found lists */

    void test_min_from_set_exact(){
        // test will find exactly bracketing k's
        std::set<int> kbrackets {0,1,3,5};
        std::map<int, double> w2brackets {{0,4},{1,12},{3,14},{5,17}};
        double w2min=0.0;
        int ktarget = 2, kmin = -1;
        mode::get_min_from_set(kbrackets, w2brackets, ktarget, kmin, w2min);
        TS_ASSERT_EQUALS(kmin, 1);
        TS_ASSERT_EQUALS(w2min, 12);
    }

    void test_max_from_set_exact(){
        // test will find exactly bracketing k's
        std::set<int> kbrackets {0,1,3,5};
        std::map<int, double> w2brackets {{0,4},{1,12},{3,14},{5,17}};
        double w2max=0.0;
        int ktarget = 2, kmax = 11;
        mode::get_max_from_set(kbrackets, w2brackets, ktarget, kmax, w2max);
        TS_ASSERT_EQUALS(kmax, 3);
        TS_ASSERT_EQUALS(w2max, 14);
    }

    void test_min_from_set_rough(){
        // test will find roughly bracketing k's
        std::set<int> kbrackets {0,1,5,11};
        std::map<int, double> w2brackets {{0,4},{1,12},{5,14},{11,17}};
        double w2min=0.0;
        int ktarget = 3, kmin = 10;
        mode::get_min_from_set(kbrackets, w2brackets, ktarget, kmin, w2min);
        TS_ASSERT_EQUALS(kmin, 1);
        TS_ASSERT_EQUALS(w2min, 12);
    }

    void test_max_from_set_rough(){
        // test will find roughly bracketing k's
        std::set<int> kbrackets {0,1,5,11};
        std::map<int, double> w2brackets {{0,4},{1,12},{5,14},{11,17}};
        double w2max=0.0;
        int ktarget = 3, kmax = 11;
        mode::get_max_from_set(kbrackets, w2brackets, ktarget, kmax, w2max);
        TS_ASSERT_EQUALS(kmax, 5);
        TS_ASSERT_EQUALS(w2max, 14);
    }

    void test_max_from_set_nothing(){
        // test will find not find a max
        std::set<int> kbrackets {0,1,3,5};
        std::map<int, double> w2brackets {{0,4},{1,12},{3,14},{5,17}};
        double w2max=0.0;
        int ktarget = 7, kmax = 11;
        mode::get_max_from_set(kbrackets, w2brackets, ktarget, kmax, w2max);
        TS_ASSERT_EQUALS(kmax, 11)
        TS_ASSERT_EQUALS(w2max, 0.0);
    }

    void test_min_from_set_nothing(){
        // test will not find a min
        std::set<int> kbrackets {5,6,7};
        std::map<int, double> w2brackets {{5,4.0},{6,12.0},{7,14.0}};
        double w2min=0.0;
        int ktarget = 3, kmin = 10;
        mode::get_min_from_set(kbrackets, w2brackets, ktarget, kmin, w2min);
        TS_ASSERT_EQUALS(kmin, 10);
        TS_ASSERT_EQUALS(w2min, 0.0);
    }

    void test_min_from_set_largest(){
        // the target is larger than all members of set
        // will return the max of the set as the min
        std::set<int> kbrackets {5,6,7};
        std::map<int, double> w2brackets {{5,4.0},{6,12.0},{7,14.0}};
        double w2min=0.0;
        int ktarget = 10, kmin = 11;
        mode::get_min_from_set(kbrackets, w2brackets, ktarget, kmin, w2min);
        TS_ASSERT_EQUALS(kmin, 7);
        TS_ASSERT_EQUALS(w2min, 14.0);
    }

    void test_max_from_set_smallest(){
        // the target is smaller than all members of set
        // will return the min of the set as the max
        std::set<int> kbrackets {5,6,7};
        std::map<int, double> w2brackets {{5,4.0},{6,12.0},{7,14.0}};
        double w2max=0.0;
        int ktarget = 2, kmax = 1;
        mode::get_max_from_set(kbrackets, w2brackets, ktarget, kmax, w2max);
        TS_ASSERT_EQUALS(kmax, 5);
        TS_ASSERT_EQUALS(w2max, 4.0);
    }

};

class FakeModeDriver : public ModeDriver {
public:
    static const unsigned int num_var=0U;
    FakeModeDriver(Star* s, double x) : ModeDriver(num_var, s) {}
    int length() { return 0; }
    double Gamma1() { return 0; }
    double rad(int x) { return 0; }
    int CentralBC(double **y, double *yo, double s2, int l, int m=0) {return 0;}
    int SurfaceBC(double **y, double *ys, double s2, int l, int m=0) {return 0;}
    void getCoeff(double *CC, const int, const int, const double, const int) {}
    void setupBoundaries(){}
    double SSR(double x, int y, ModeBase* ) {return 0.0;}
    double tidal_overlap(ModeBase*) {return 0.0;}
    double innerproduct(ModeBase*,ModeBase*) {return 0.0;}
    void getBoundaryMatrix(int, double *, double*, double **, int*){};
	void varnames(std::string*){};	//names of variables to print out
};

class FakeMode : public ModeBase {
public:
    static const unsigned int num_var=0U;
    void printMode(const char *const c){}
    void writeMode(const char *const c){}
    double getRad(int x){return 0.0;}
    double getY(int i, int x){return 0.0;}
    double getYtilde(int i, int x){return 0.0;}
    void modeNumbers(int& x, int& y, int& z){}
    double SSR(){return 0.0;}
    double tidal_overlap(){return 0.0;}
    double getFreq() {return 0.0;}
    double getPeriod() {return 0.0;}

    // the ones we intend to mock
    int modeOrder(){return K;}
    double getOmega2(){return (K==0? ZEROW2 : double(K));}
    FakeMode() {}
    FakeMode(int k, int l, int m, ModeDriver* drv) {
        K = klist[iter];
        iter = (iter+1)%klist.size();
    }
    FakeMode(double w2, int l, int m, ModeDriver* drv) {
        K = int(floor(w2));
    }
    FakeMode(double w1, double w2, int l, int m, ModeDriver* drv) {
        K = int((w2+w1)/2);
    }

    static int iter;
    static std::vector<int> klist;
private:
    int K;
};
int FakeMode::iter(0);
std::vector<int> FakeMode::klist({});

class ModeFinderTest : public CxxTest::TestSuite {
public:

    // this function will make output files for better debugging
    static ModeFinderTest *createSuite (){
        system("mkdir -p ./output/../tests/modefinder/../tests");
        system("rm ./output/../tests/modefinder/../tests/modefinder.txt");
        system("touch ./output/../tests/modefinder/../tests/modefinder.txt");
        return new ModeFinderTest();
    }
    static void destroySuite(ModeFinderTest *suite) { 
        delete suite; 
    }

    Calculation::OutputData mode_finder_test_setup(
        const std::vector<int>& L,
        const std::vector<int>& K
    ){
        Calculation::OutputData data;
        data.calcname = "../tests/modefinder";
        data.freq0 = 1.0;
        data.star = nullptr;
        data.driver = new FakeModeDriver(data.star, 1.5);
        data.i_err = 0;
        for(int i=0; i<4; i++)
            data.error[i] = false;
        
        data.mode_num = K.size();
        data.mode_writ = data.mode_done = 0;
        data.l = L;
        data.k = K;
        data.w.reserve(data.mode_num);
        data.f.reserve(data.mode_num);
        data.mode.reserve(data.mode_num);
        data.period.reserve(data.mode_num);
        data.mode_SSR.reserve(data.mode_num);
        return data;
    }

    void do_test_mode_finder_fake_classes(
        const std::string test_name,
        const std::vector<int>& L, 
        const std::vector<int>& K,
        const std::vector<int>& Lexp,
        const std::vector<int>& Kexp,
        const std::vector<double>& w2exp
    ){
        Calculation::OutputData data = mode_finder_test_setup(L,K);
        FILE* modetest = fopen("./tests/tests/modefinder.txt", "a");
        fprintf(modetest, "\n#  TEST NAME: %s\n", test_name.c_str());
        fclose(modetest);
        printf("TEST %s\n", test_name.c_str());
        mode::mode_finder<FakeMode, FakeModeDriver>(data);

        TS_ASSERT_EQUALS(data.mode_num, Kexp.size());
        TS_ASSERT_EQUALS(data.k.size(), Kexp.size());
        TS_ASSERT_EQUALS(data.l.size(), Kexp.size());
        TS_ASSERT_EQUALS(data.k.size(), Lexp.size());
        for(std::size_t i=0; i<Kexp.size(); i++){
            TS_ASSERT_EQUALS(data.l[i], Lexp[i]);
            TS_ASSERT_EQUALS(data.k[i], Kexp[i]);
            TS_ASSERT_EQUALS(data.w[i], sqrt(double(w2exp[i])));
        }
    }

    // make sure that it does not try to look for a dipole f-mode
    void test_mode_finder_no_dipole_f(){
        std::string testname = "NO DIPOLE F";
        // STEP 1 removes the L=1,K=0 mode
        // STEP 2 finds only K=1,2 modes
        FakeMode::iter = 0;
        FakeMode::klist ={2,1};

        // ask for the dipole f-mode
        std::vector<int> L = {1,1};
        std::vector<int> K = {0,2};
        // expect to not find one
        std::vector<int> Lexp = {1,1};
        std::vector<int> Kexp = {1,2};
        std::vector<double> w2 = {1,2};

        do_test_mode_finder_fake_classes(testname, L,K,Lexp,Kexp, w2);
    }

    // make sure that it always looks for f-mode with L>=2
    void test_mode_akways_find_k0(){
        std::string testname = "ALWAYS FIND K0 - L=", name;
        // STEP 1 adds a K=0 mode for K>=2
        // STEP 2 finds K=0,1 modes
        FakeMode::iter = 0;
        FakeMode::klist ={0,1};

        std::vector<int> Kin {1};
        std::vector<int> Kout {0,1};
        std::vector<double> w2 {ZEROW2,1};

        // ask for several L values. expect the K=0 f-mode to be added in
        for(int L=2; L<10; L++){
            std::vector<int> Lin {L};
            std::vector<int> Lout {L,L};
            name = testname + std::to_string(L);
            do_test_mode_finder_fake_classes(name, Lin,Kin, Lout,Kout, w2);
        }
    }

    // simplest test, finds modes in STEP 2 in order
    void test_mode_finder_happy_path(){
        std::string testname = "HAPPY PATH";
        // STEP 1 does nothing
        // STEP 2 finds K=1,2,3,4,5 in order
        FakeMode::iter = 0;
        FakeMode::klist = {1,2,3,4,5};

        std::vector<int> L = {1,1,1,1,1};
        std::vector<int> K = {1,2,3,4,5};
        std::vector<double> w2 = {1,2,3,4,5};
        do_test_mode_finder_fake_classes(testname, L,K,L,K, w2);
    }

    // simple test, finds modes in STEP 2 out of order
    void test_mode_finder_happy_path_unordered(){
        std::string testname = "HAPPY PATH UNORDERED";
        // STEP 1 puts the unordered K in order
        // STEP 2 finds K=1,2,3,4,5 in order
        FakeMode::iter = 0;
        FakeMode::klist = {4,3,5,2,1};

        std::vector<int> L = {1,1,1,1,1};
        std::vector<int> K = {1,2,3,4,5};
        std::vector<double> w2 = {1,2,3,4,5};
        do_test_mode_finder_fake_classes(testname, L,K,L,K, w2);
    }

    // chck multiple Ls
    void test_mode_finder_many_L(){
        std::string testname = "MANY Ls";
        // STEP 1 creates maps for each L, adding 0s in L=2,3
        // STEP 2 finds K=0,1,2,3 for each L
        FakeMode::iter = 0;
        FakeMode::klist = {0,1,2,3};

        std::vector<int> L = {1,1,1, 2,2,2, 3,3,3};
        std::vector<int> K = {1,2,3, 1,2,3, 1,2,3};
        std::vector<int> Lexp = {1,1,1, 2,2,2,2, 3,3,3,3};
        std::vector<int> Kexp = {1,2,3, 0,1,2,3, 0,1,2,3};
        std::vector<double> w2 = {1,2,3, ZEROW2,1,2,3, ZEROW2,1,2,3};
        do_test_mode_finder_fake_classes(testname, L,K,Lexp,Kexp, w2);
    }

    void test_mode_finder_discover_in_bracket_search(){
        std::string testname = "FIND BRACKET IN SEACH";
        // STEP 1 create maps
        // STEP 2 all modes found except for 3
        // STEP 3 begin looking for brackets on 3, and find 3
        FakeMode::iter = 0;
        FakeMode::klist ={1,2,2, 3};

        std::vector<int> L = {1,1,1};
        std::vector<int> K = {1,2,3};
        std::vector<double> w2 = {1,2,3};
        do_test_mode_finder_fake_classes(testname, L,K,L,K, w2);
    }

    void test_mode_finder_exact_brackets(){
        std::string testname = "EXACT BRACKETS";
        // STEP 1 create maps
        // STEP 2 find modes to be brackets
        // STEP 3 setup exact brackets (K=1,3,5)
        // STEP 4 use exact brackets to find K=2,4
        FakeMode::iter = 0;
        // the mode will never produce K=2,4 except in bracketed constructor
        FakeMode::klist ={1,3,5,7};

        std::vector<int> L = {1,1,1};
        std::vector<int> K = {1,2,4};
        std::vector<double> w2 = {1,2,4};
        do_test_mode_finder_fake_classes(testname, L,K,L,K, w2);
    }

    void test_mode_finder_near_brackets(){
        std::string testname = "NEAR BRACKETS";
        // STEP 1 create maps
        // STEP 2 modes to benear brackets are found
        // STEP 3 setup near brackets from earlier
        // STEP 4 use near brackets to converge
        FakeMode::iter = 0;
        FakeMode::klist ={1,7,13,12};

        std::vector<int> L = {1,1,1};
        std::vector<int> K = {1,4,9};
        std::vector<double> w2 = {1,4,9};
        do_test_mode_finder_fake_classes(testname, L,K,L,K, w2);
    }

    void test_mode_finder_find_brackets(){
        std::string testname = "FINDING BRACKETS";
        // STEP 1 create maps
        // STEP 2 find K=1, brackets for K=4, not K=9
        // STEP 3 K=1 skip, K=4 brackets exist and find K=10, K=9 now has K=10 as bracket
        // STEP 4 K=1 skip, K=4,9 use brackets
        FakeMode::iter = 0;
        FakeMode::klist ={1,1,19, 11,};

        std::vector<int> L = {1,1,1};
        std::vector<int> K = {1,4,9};
        std::vector<double> w2 = {1,4,9};
        do_test_mode_finder_fake_classes(testname, L,K,L,K, w2);
    }

    void test_mode_finder_set_lower(){
        std::string testname = "TEST FIND LOWER";
        // STEP 2: find a mode to serve as upper bracket
        // STEP 3: setup upper, find default lower
        // STEP 4: use brackets to converge
        FakeMode::iter = 0;
        // all found modes are higher than those requested
        FakeMode::klist ={5,5,  5};

        std::vector<int> L = {1,1};
        std::vector<int> K = {1,3};
        std::vector<double> w2 = {1,3};
        do_test_mode_finder_fake_classes(testname, L,K,L,K, w2);
    }

    void test_mode_finder_set_upper(){
        std::string testname = "TEST FIND UPPER";
        // STEP 2: find a mode to serve as lower bracket
        // STEP 3: setup lower, quest for upper
        // STEP 4: use brackets to converge
        FakeMode::iter = 0;
        // all found modes are lower than those requested
        FakeMode::klist ={1,1,  1};

        std::vector<int> L = {1,1};
        std::vector<int> K = {1,3};
        std::vector<double> w2 = {1,3};
        do_test_mode_finder_fake_classes(testname, L,K,L,K, w2);
    }

    void test_mode_finder_unordered(){
        std::string testname = "TEST UNORDERED";
        FakeMode::iter = 0;
        // STEP 2 finds nothing in list
        // STEP 3 must search for brackets
        // STEP 4 will find everything (needs to be given one, due to round-up)
        FakeMode::klist ={9,9,8,8,7};

        std::vector<int> L = {1,1,1,1,1};
        std::vector<int> K = {1,2,3,4,5};
        std::vector<double> w2 = {1,2,3,4,5};
        do_test_mode_finder_fake_classes(testname, L,K,L,K, w2);
    }
};

class ControlledMode : public FakeMode {
// this one is like the previous, except we can always control which modes appear
public:
    static const unsigned int num_var=0U;

    // the ones we intend to mock
    int modeOrder(){return K;}
    double getOmega2(){return (K==0? ZEROW2 : double(K));}
    ControlledMode(int k, int l, int m, ModeDriver* drv) {
        K = klist[iter];
        iter = (iter+1)%klist.size();
    }
    ControlledMode(double w2, int l, int m, ModeDriver* drv) {
        K = klist[iter];
        iter = (iter+1)%klist.size();
    }
    ControlledMode(double w1, double w2, int l, int m, ModeDriver* drv) {
        K = klist[iter];
        iter = (iter+1)%klist.size();
    }

    static int iter;
    static std::vector<int> klist;
private:
    int K;

};
int ControlledMode::iter(0);
std::vector<int> ControlledMode::klist({});

class ControlledModeFinderTest : public ModeFinderTest {
public:

    // this will destroy one of the extraneous dirs created
    static ControlledModeFinderTest *createSuite (){
        return new ControlledModeFinderTest();
    }
    static void destroySuite(ControlledModeFinderTest *suite) { 
        system("rm -r ./tests/modefinder/");
        delete suite; 
    }

    void do_test_mode_finder_fake_classes(
        const std::string test_name,
        const std::vector<int>& L, 
        const std::vector<int>& K,
        const std::vector<int>& Lexp,
        const std::vector<int>& Kexp,
        const std::vector<double>& w2exp
    ){
        Calculation::OutputData data = mode_finder_test_setup(L,K);
        FILE* modetest = fopen("./tests/tests/modefinder.txt", "a");
        fprintf(modetest, "\n#  TEST NAME: %s\n", test_name.c_str());
        fclose(modetest);
        printf("TEST %s\n", test_name.c_str());
        mode::mode_finder<ControlledMode, FakeModeDriver>(data);

        TS_ASSERT_EQUALS(data.mode_num, Kexp.size());
        TS_ASSERT_EQUALS(data.k.size(), Kexp.size());
        TS_ASSERT_EQUALS(data.l.size(), Kexp.size());
        TS_ASSERT_EQUALS(data.k.size(), Lexp.size());
        for(std::size_t i=0; i<Kexp.size(); i++){
            TS_ASSERT_EQUALS(data.l[i], Lexp[i]);
            TS_ASSERT_EQUALS(data.k[i], Kexp[i]);
            TS_ASSERT_EQUALS(data.w[i], sqrt(double(w2exp[i])));
        }
    }

    void test_controlled_mode_finder_happy_path(){
        std::string testname = "CONTROLLED HAPPY PATH";
        // STEP 1 does nothing
        // STEP 2 finds K=1,2,3,4,5 in order
        ControlledMode::iter = 0;
        ControlledMode::klist = {1,2,3,4,5};

        std::vector<int> L = {1,1,1,1,1};
        std::vector<int> K = {1,2,3,4,5};
        std::vector<double> w2 = {1,2,3,4,5};
        do_test_mode_finder_fake_classes(testname, L,K,L,K, w2);
    }

    void test_controlled_mode_finder_random(){
        std::string testname = "CONTROLLED RANDOM";
        ControlledMode::iter = 0;
        ControlledMode::klist = {
            // STEP 2 will not find L=2, but will find brackets
            5,1, 
            // STEP 3 uses K=1,3 from above as brackets
            // STEP 4
            // STEP 4b will find K=5
            5,
            // STEP 4c will be skipped
            // STEP 4d will fail to move brackets
            // STEP 4e will be skipped
            // STEP 4f will try random position, find K=4
            4,
            // STEP 4b creates a new trial mode, will find K=4
            4,
            // STEP 4c - STEP 4f as above, find K=2 on second try
            4, 2
        };

        std::vector<int> L = {1};
        std::vector<int> K = {2};
        std::vector<int> Lexp = {1,1};
        std::vector<int> Kexp = {1,2};
        std::vector<double> w2 = {1,2};
        do_test_mode_finder_fake_classes(testname, L,K,Lexp,Kexp, w2);
    }

    void test_controlled_mode_finder_straggler(){
        std::string testname = "CONTROLLED STRAGGLER";
        ControlledMode::iter = 0;
        ControlledMode::klist = {
            // STEP 2 will not find L=2, but will find brackets
            5,1, 
            // STEP 3 uses K=1,3 from above as brackets
            // STEP 4b will find K=5
            5,
            // STEP 4c will be skipped
            // STEP 4d will fail to move brackets
            // STEP 4e will be skipped
            // STEP 4f will try random position, stuck on K=5
            5,5,5,5,5,
            // STEP 5a will be skipped
            // STEP 5b (i) will find brackets K=1,5
            // STEP 5b (ii) will try to bisect, then find K=2
            5, 2
        };


        std::vector<int> L = {1};
        std::vector<int> K = {2};
        std::vector<int> Lexp = {1,1};
        std::vector<int> Kexp = {1,2};
        std::vector<double> w2 = {1,2};
        do_test_mode_finder_fake_classes(testname, L,K,Lexp,Kexp, w2);
    }

    void test_controlled_mode_finder_unfindable(){
        std::string testname = "CONTROLLED UNFINDABLE";
        ControlledMode::iter = 0;
        ControlledMode::klist = { 1, 5 };

        std::vector<int> L = {1};
        std::vector<int> K = {2};
        std::vector<int> Lexp = {1,1};
        std::vector<int> Kexp = {1,2};
        std::vector<double> w2 = {1,0};
        do_test_mode_finder_fake_classes(testname, L,K,Lexp,Kexp, w2);
    }

};
