#include "../src/ThrainMain.h"
#include "../src/ThrainUnits.h"
#include <cxxtest/TestSuite.h>

class UnitsBaseTest : public CxxTest::TestSuite {
public:

    static UnitsBaseTest *createSuite (){
        printf("\nUNIT TESTS");
        return new UnitsBaseTest;
    }
    static void destroySuite(UnitsBaseTest *suite) { 
        delete suite; 
    }

    void setUp() {
        freopen("tests/artifacts/logio.txt", "a", stdout);
    }

    void tearDown() {
        freopen("/dev/tty", "w", stdout);
    }

    Calculation::OutputData setupFakeCalcData(units::Units unitType){
        Calculation::OutputData fakeData;
        fakeData.star = nullptr;
        fakeData.driver = nullptr;
        fakeData.i_err = 0;

        char fakeParams = units::ParamType::pmass|units::ParamType::pradius;
        double fakeMass = 1.0;
        // must be fluffy enough for zsurf to be real
        double fakeRadius = 10.0;
        fakeData.units = unitType;        
        fakeData.params = fakeParams;
        fakeData.mass = fakeMass;
        fakeData.radius = fakeRadius;
        fakeData.freq0 = 0.0;

        units::UnitSet fakeUnits {1.0, 2.0, 8.0, 7.0, 6.0};
        fakeData.unitset = fakeUnits;
        // ensure the fake untis were set with fake numbers
        TS_ASSERT(fakeData.unitset.G == 1.0);
        TS_ASSERT(fakeData.unitset.C == 2.0);
        TS_ASSERT(fakeData.unitset.base_length == 8.0);
        TS_ASSERT(fakeData.unitset.base_time == 7.0);
        TS_ASSERT(fakeData.unitset.base_mass == 6.0);

        units::format_units(fakeData);
        // assert nothing unchanging has changed
        TS_ASSERT(fakeData.mass == fakeMass);
        TS_ASSERT(fakeData.radius == fakeRadius);
        TS_ASSERT(fakeData.params == fakeParams);

        return fakeData;
    }

    void do_test_setUnits_noparams(units::Units unitType){

        // test with an undefined params
        Calculation::OutputData fakeData1;
        fakeData1.star = nullptr;
        fakeData1.driver = nullptr;
        fakeData1.i_err = 0;
        fakeData1.units = units::Unit::astro;
        TS_ASSERT_THROWS_NOTHING(units::format_units(fakeData1));
        // assert units safely zeroed out
        // TS_ASSERT(fakeData1.unitset.G == 0.0);
        // TS_ASSERT(fakeData1.unitset.C == 0.0);
        // TS_ASSERT(fakeData1.unitset.base_length == 0.0);
        // TS_ASSERT(fakeData1.unitset.base_time == 0.0);
        // TS_ASSERT(fakeData1.unitset.base_mass == 0.0);
        // assert all safely zeroed out
        TS_ASSERT(fakeData1.mass == 0.0);
        TS_ASSERT(fakeData1.radius == 0.0);
        TS_ASSERT(fakeData1.zsurf == 0.0);
        TS_ASSERT(fakeData1.logg == 0.0);

        // test with a zero-set params
        Calculation::OutputData fakeData2;
        fakeData2.star = nullptr;
        fakeData2.driver = nullptr;
        fakeData2.i_err = 0;
        fakeData2.params = 0;
        fakeData1.units = units::Units::astro;
        TS_ASSERT_THROWS_NOTHING(units::format_units(fakeData2));
        // assert all safely zeroed out
        TS_ASSERT(fakeData2.mass == 0.0);
        TS_ASSERT(fakeData2.radius == 0.0);
        TS_ASSERT(fakeData2.zsurf == 0.0);
        TS_ASSERT(fakeData2.logg == 0.0);
    }    

    void do_test_setUnits(units::Units unitType){
        Calculation::OutputData fakeData = setupFakeCalcData(unitType);
        // assert units set properly
        units::UnitSet theseUnits = units::unitSets.at(unitType); 
        TS_ASSERT(fakeData.unitset.G == theseUnits.G);
        TS_ASSERT(fakeData.unitset.C == theseUnits.C);
        TS_ASSERT(fakeData.unitset.base_length == theseUnits.base_length);
        TS_ASSERT(fakeData.unitset.base_time == theseUnits.base_time);
        TS_ASSERT(fakeData.unitset.base_mass == theseUnits.base_mass);
        // assert other input params updated
        double mcgs = fakeData.unitset.base_mass*fakeData.mass;
        double rcgs = fakeData.unitset.base_length*fakeData.radius;
        double logg = log10(G_CGS*mcgs*pow(rcgs,-2));
        double zsurf = 1./sqrt(1. - 2.*G_CGS*mcgs/(rcgs*pow(C_CGS,2))) - 1.0;
        TS_ASSERT_DELTA(fakeData.zsurf, zsurf, 1e-4);
        TS_ASSERT(fakeData.logg == logg);
        TS_ASSERT(fakeData.freq0 == sqrt(G_CGS*mcgs*pow(rcgs,-3)));
    }

    void test_setAllUnits_noparams(){
        do_test_setUnits_noparams(units::Units::astro);
        do_test_setUnits_noparams(units::Units::geo);
        do_test_setUnits_noparams(units::Units::SI);
        do_test_setUnits_noparams(units::Units::CGS);
    }

    void test_setAllUnits(){
        do_test_setUnits(units::Units::astro);
        do_test_setUnits(units::Units::geo);
        do_test_setUnits(units::Units::SI);
        do_test_setUnits(units::Units::CGS);
    }

    // test astro explicitly
    void test_setUnitsAstro() {
        Calculation::OutputData fakeData = setupFakeCalcData(units::Units::astro);
        // assert units set properly
        TS_ASSERT(fakeData.unitset.G == G_astro);
        TS_ASSERT(fakeData.unitset.C == C_astro);
        TS_ASSERT(fakeData.unitset.base_length == 1.0e5);
        TS_ASSERT(fakeData.unitset.base_time == 1.0);
        TS_ASSERT(fakeData.unitset.base_mass == MSOLAR);
    }

    // test geo explicitly
    void test_setUnitsGeo() {
        Calculation::OutputData fakeData = setupFakeCalcData(units::Units::geo);
        // assert units set properly
        TS_ASSERT(fakeData.unitset.G == 1.0);
        TS_ASSERT(fakeData.unitset.C == 1.0);
        TS_ASSERT(fakeData.unitset.base_length == 1.0/C_CGS);
        TS_ASSERT(fakeData.unitset.base_time == 1.0);
        TS_ASSERT(fakeData.unitset.base_mass == C_CGS/G_CGS);
    }

    // test all the param setter functions
    void test_paramSetters(){
        // use geometric units, because easier to predict
        Calculation::OutputData fakeData = setupFakeCalcData(units::Units::geo);

        double R = fakeData.radius;
        double M = fakeData.mass;

        // test making logg from R and M
        double res1 = units::getLoggFromRM(fakeData);
        double exp1 = log10( M/R/R *pow(C_CGS,3) );
        TS_ASSERT_EQUALS( res1, exp1 );

        // test making zsurf from R and M
        double res2 = units::getZsurfFromRM(fakeData);
        double exp2 = 1./sqrt(1. - 2.*M/R) - 1.;
        TS_ASSERT_EQUALS( res2, exp2 );

        // test making radius from zsurf and M
        fakeData.zsurf = units::getZsurfFromRM(fakeData);
        double radius3 = units::getRadiusFromZM(fakeData);
        TS_ASSERT_DELTA( radius3, fakeData.radius, 1.e-10 );

        // test making mass from zsurf and R
        fakeData.zsurf = units::getZsurfFromRM(fakeData);
        double mass4 = units::getMassFromRZ(fakeData);
        TS_ASSERT_DELTA( mass4, fakeData.mass, 1.e-10 );

        // test getting GM/R^2 from logg
        double res5 = units::GMR2FromLogg(fakeData);
        double exp5 = M/R/R;
        TS_ASSERT_DELTA( res5, exp5, 1.e-10 );

        // test making radius from logg and M
        fakeData.logg = units::getLoggFromRM(fakeData);
        double radius6 = units::getRadiusFromLoggM(fakeData);
        TS_ASSERT_DELTA( radius6, fakeData.radius, 1.e-10 );

        // test making mass from logg and R
        fakeData.logg = units::getLoggFromRM(fakeData);
        double mass7 = units::getMassFromRLogg(fakeData);
        TS_ASSERT_DELTA( mass7, fakeData.mass, 1.e-10 );  
    }
};