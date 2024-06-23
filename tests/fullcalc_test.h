#include "../src/ThrainMain.h"
#include "../src/ThrainMode.h"
#include "../src/ThrainUnits.h"
#include "../src/ThrainIO.h"
#include <cxxtest/TestSuite.h>

class FullCalculationTest : public CxxTest::TestSuite {
public:

    Calculation::InputData make_input_data_pmodes(std::string testname){
        Calculation::InputData in;
        in.calcname = testname;
        in.regime = regime::PN0;
        in.model = model::polytrope;
        in.units = units::Units::astro;
        in.modetype = modetype::nonradial;
        in.input_params = {0.0, 5e4};
        in.mode_num = 17;
        in.Ngrid = std::size_t(in.input_params[1]);
        in.l = {1,2,3};
        in.kl = {
            {1,{  1,2,3,4,5}},
            {2,{0,1,2,3,4,5}},
            {3,{0,1,2,3,4,5}}
        };
        in.mass = 1.0;
        in.radius = 1.2;
        in.params = units::ParamType::pmass|units::ParamType::pradius;
        in.adiabatic_index = 5./3.;
        return in;
    }

    Calculation::InputData make_input_data_gmodes(std::string testname){
        Calculation::InputData in;
        in.calcname = testname;
        in.regime = regime::PN0;
        in.model = model::polytrope;
        in.units = units::Units::astro;
        in.modetype = modetype::nonradial;
        in.input_params = {0.0, 5e4};
        in.mode_num = 18;
        in.Ngrid = std::size_t(in.input_params[1]);
        in.l = {1,2,3};
        in.kl = {
            {1,{1,-1,-2,-3,-4,-5}},
            {2,{0,-1,-2,-3,-4,-5}},
            {3,{0,-1,-2,-3,-4,-5}}
        };
        in.mass = 1.0;
        in.radius = 1.2;
        in.params = units::ParamType::pmass|units::ParamType::pradius;
        in.adiabatic_index = 5./3.;
        return in;
    }

    void read_entire_file(std::string filename, std::string& contents){
        FILE* infile = fopen(filename.c_str(), "r");
        if(!infile) {
            TS_FAIL("could not read in indicated file\n");
        }
        fseek(infile, 0, SEEK_END);
        std::size_t fsize = ftell(infile);
        fseek(infile, 0, SEEK_SET);

        char *read_buffer = new char[fsize + 1];
        fread(read_buffer, fsize, 1, infile);
        read_buffer[fsize] = 0;
        fclose(infile);
        contents = std::string(read_buffer);
        delete[] read_buffer;
    }

    // /* test uniform density (n=0) polytrope */
    // void test_full_calculation_uniform() {
    //     printf("\nTEST CALCULATION UNIFORM STAR");
    //     system("mkdir -p ./output/../tests/uniform/../tests");
    //     system("touch ./output/../tests/modefinder/../tests/uniform.txt");
    //     Calculation::InputData in = make_input_data_pmodes("../tests/uniform");
    //     Calculation::OutputData out;
    //     TS_ASSERT_EQUALS(0, io::setup_output(in, out));
    //     TS_ASSERT_EQUALS(0, create_star(out));
    //     TS_ASSERT_EQUALS(0, mode::create_modes(out));

    //     TS_ASSERT_LESS_THAN(out.star_SSR, 1.e-12);
    //     for(int i=0; i<out.mode_num; i++){
    //         printf("%d,%d\t", out.l[i], out.k[i]);
    //         if(out.k[i]!=0){
    //             printf("%le %le", out.err[0][i], out.err[1][i]);
    //             TS_ASSERT_LESS_THAN(out.err[0][i], 1.e-8);
    //             TS_ASSERT_LESS_THAN(out.err[1][i], 1.e-6);
    //         }
    //         printf("\n");
    //     }
    // }

//     /* test n=1 polytrope */
//     void test_full_calculation_polytrope_15() {
//         printf("TEST CALCULATION POLYTROPE n=1.5\n");
//         Calculation::InputData in = make_input_data_pmodes("../tests/poly1.5");
//         in.input_params[0] = 1.5; // change polytrope index
//         Calculation::OutputData out;
//         TS_ASSERT_EQUALS(0, io::setup_output(in, out));
//         TS_ASSERT_EQUALS(0, create_star(out));
//         TS_ASSERT_EQUALS(0, mode::create_modes(out));

//         TS_ASSERT_LESS_THAN(out.star_SSR, 1.e-12);
//         for(int i=0; i<out.mode_num; i++){
//             printf("%d,%d\t", out.l[i], out.k[i]);
//             printf("%le %le\n", out.err[0][i], out.err[1][i]);
//             TS_ASSERT_LESS_THAN(out.err[0][i], 1.e-8);
//             if(out.k[i]!=0){
//                 TS_ASSERT_LESS_THAN(out.err[1][i], 2.e-4);
//             }
//         }
//     }

//     /* test n=3 polytrope */
//     void test_full_calculation_polytrope_30_pmode() {
//         printf("TEST CALCULATION POLYTROPE n=3.0\n");
//         Calculation::InputData in = make_input_data_pmodes("../tests/poly3.0");
//         in.input_params[0] = 3.0; // change polytrope index
//         Calculation::OutputData out;
//         TS_ASSERT_EQUALS(0, io::setup_output(in, out));
//         TS_ASSERT_EQUALS(0, create_star(out));
//         TS_ASSERT_EQUALS(0, mode::create_modes(out));

//         std::map<int,ModeBase*> fmodes;

//         TS_ASSERT_LESS_THAN(out.star_SSR, 1.e-12);
//         for(int i=0; i<out.mode_num; i++){
//             printf("%d,%d\t", out.l[i], out.k[i]);
//             printf("%le %le\n", out.err[0][i], out.err[1][i]);
//             TS_ASSERT_LESS_THAN(out.err[0][i], 1.e-8);
//             if(out.k[i]!=0){
//                 TS_ASSERT_LESS_THAN(out.err[1][i], 2.e-4);
//             }

//             if(out.l[i]==1 && out.k[i]==1){
//                 fmodes[out.l[i]] = out.mode[i];
//             }
//             if(out.l[i]==2 && out.k[i] == 0){
//                 if(out.k[i]==0) fmodes[out.l[i]]=out.mode[i];
//             }
//             if(out.l[i]==3 && out.k[i] == 0){
//                 if(out.k[i]==0) fmodes[out.l[i]]=out.mode[i];
//             }
//         }
//         //test the c0 values
//         for(int i=0; i<out.mode_num; i++){
//             if( 
//                 (out.l[i]==1 && out.k[i]!=1) ||
//                 (out.l[i]==2 && out.k[i]!=0) ||
//                 (out.l[i]==3 && out.k[i]!=0)
//             ) {
//                 TS_ASSERT_LESS_THAN(out.driver->innerproduct(out.mode[i], fmodes[out.l[i]]), 1.e-10);
//             }
//         }
//     }
//     void test_full_calculation_polytrope_30_gmode() {
//         printf("TEST CALCULATION POLYTROPE n=3.0 GMODE\n");
//         Calculation::InputData in = make_input_data_gmodes("../tests/poly3.0g");
//         in.input_params[0] = 3.0; // change polytrope index
//         Calculation::OutputData out;
//         TS_ASSERT_EQUALS(0, io::setup_output(in, out));
//         TS_ASSERT_EQUALS(0, create_star(out));
//         TS_ASSERT_EQUALS(0, mode::create_modes(out));

//         std::map<int,ModeBase*> fmodes;

//         TS_ASSERT_LESS_THAN(out.star_SSR, 1.e-12);
//         for(int i=0; i<out.mode_num; i++){
//             printf("%d,%d\t", out.l[i], out.k[i]);
//             printf("%le %le\n", out.err[0][i], out.err[1][i]);
//             TS_ASSERT_LESS_THAN(out.err[0][i], 1.e-8);

//             if(out.l[i]==1 && out.k[i]==1){
//                 fmodes[out.l[i]] = out.mode[i];
//             }
//             if(out.l[i]==2 && out.k[i] == 0){
//                 if(out.k[i]==0) fmodes[out.l[i]]=out.mode[i];
//             }
//             if(out.l[i]==3 && out.k[i] == 0){
//                 if(out.k[i]==0) fmodes[out.l[i]]=out.mode[i];
//             }
//         }
//         //test the c0 values
//         for(int i=0; i<out.mode_num; i++){
//             if( 
//                 (out.l[i]==1 && out.k[i]!=1) ||
//                 (out.l[i]==2 && out.k[i]!=0) ||
//                 (out.l[i]==3 && out.k[i]!=0)
//             ) {
//                 TS_ASSERT_LESS_THAN(out.driver->innerproduct(out.mode[i], fmodes[out.l[i]]), 1.e-10);
//             }
//         }
//     }

//     /* test CHWD with uniform mu */
//     void test_full_calculation_simple_CHWD() {
//         printf("TEST CALCULATION CHWD SIMPLE\n");
//         Calculation::InputData in = make_input_data_pmodes("../tests/chwd_0");
//         in.model = model::CHWD;
//         double temp = in.input_params[1];
//         in.input_params = {1.581, 0, temp};
//         in.adiabatic_index = 0;
//         Calculation::OutputData out;
//         TS_ASSERT_EQUALS(0, io::setup_output(in, out));
//         TS_ASSERT_EQUALS(0, create_star(out));
//         TS_ASSERT_EQUALS(0, mode::create_modes(out));

//         std::map<int,ModeBase*> fmodes;

//         TS_ASSERT_LESS_THAN(out.star_SSR, 1.e-12);
//         for(int i=0; i<out.mode_num; i++){
//             printf("%d,%d\t", out.l[i], out.k[i]);
//             printf("%le\n", out.err[0][i]);
//             TS_ASSERT_LESS_THAN(out.err[0][i], 1.e-8);

//             if(out.l[i]==1 && out.k[i]==1){
//                 fmodes[out.l[i]] = out.mode[i];
//             }
//             if(out.l[i]==2 && out.k[i] == 0){
//                 if(out.k[i]==0) fmodes[out.l[i]]=out.mode[i];
//             }
//             if(out.l[i]==3 && out.k[i] == 0){
//                 if(out.k[i]==0) fmodes[out.l[i]]=out.mode[i];
//             }
//         }
//         //test the c0 values
//         for(int i=0; i<out.mode_num; i++){
//             if( 
//                 (out.l[i]==1 && out.k[i]!=1) ||
//                 (out.l[i]==2 && out.k[i]!=0) ||
//                 (out.l[i]==3 && out.k[i]!=0)
//             ) {
//                 TS_ASSERT_LESS_THAN(out.driver->innerproduct(out.mode[i], fmodes[out.l[i]]), 1.e-10);
//             }
//         }
//     }

//     /* test CHWD with uniform mu */
//     void test_full_calculation_sigmoidal_CHWD_pmode() {
//         printf("TEST CALCULATION CHWD SIGMOIDAL\n");
//         Calculation::InputData in = make_input_data_pmodes("../tests/chwd_1_p");
//         in.model = model::CHWD;
//         double temp = in.input_params[1];
//         in.input_params = {1.581, 1, temp};
//         Calculation::OutputData out;
//         TS_ASSERT_EQUALS(0, io::setup_output(in, out));
//         TS_ASSERT_EQUALS(0, create_star(out));
//         TS_ASSERT_EQUALS(0, mode::create_modes(out));

//         std::map<int,ModeBase*> fmodes;

//         TS_ASSERT_LESS_THAN(out.star_SSR, 1.e-02);
//         for(int i=0; i<out.mode_num; i++){
//             printf("%d,%d\t", out.l[i], out.k[i]);
//             printf("%le\n", out.err[0][i]);
//             TS_ASSERT_LESS_THAN(out.err[0][i], 1.e-8);

//             if(out.l[i]==1 && out.k[i]==1){
//                 fmodes[out.l[i]] = out.mode[i];
//             }
//             if(out.l[i]==2 && out.k[i] == 0){
//                 if(out.k[i]==0) fmodes[out.l[i]]=out.mode[i];
//             }
//             if(out.l[i]==3 && out.k[i] == 0){
//                 if(out.k[i]==0) fmodes[out.l[i]]=out.mode[i];
//             }
//         }
//         //test the c0 values
//         for(int i=0; i<out.mode_num; i++){
//             if( 
//                 (out.l[i]==1 && out.k[i]!=1) ||
//                 (out.l[i]==2 && out.k[i]!=0) ||
//                 (out.l[i]==3 && out.k[i]!=0)
//             ) {
//                 TS_ASSERT_LESS_THAN(out.driver->innerproduct(out.mode[i], fmodes[out.l[i]]), 1.e-10);
//             }
//         }
//     }

//     void test_full_calculation_sigmoidal_CHWD_gmode() {
//         printf("TEST CALCULATION CHWD SIGMOIDAL GMODE\n");
//         Calculation::InputData in = make_input_data_gmodes("../tests/chwd_1_g");
//         in.model = model::CHWD;
//         double temp = in.input_params[1];
//         in.input_params = {1.581, 1, temp};
//         Calculation::OutputData out;
//         TS_ASSERT_EQUALS(0, io::setup_output(in, out));
//         TS_ASSERT_EQUALS(0, create_star(out));
//         TS_ASSERT_EQUALS(0, mode::create_modes(out));

//         std::map<int,ModeBase*> fmodes;

//         TS_ASSERT_LESS_THAN(out.star_SSR, 1.e-02);
//         for(int i=0; i<out.mode_num; i++){
//             printf("%d,%d\t", out.l[i], out.k[i]);
//             printf("%le\n", out.err[0][i]);
//             TS_ASSERT_LESS_THAN(out.err[0][i], 1.e-8);

//             if(out.l[i]==1 && out.k[i]==1){
//                 fmodes[out.l[i]] = out.mode[i];
//             }
//             if(out.l[i]==2 && out.k[i] == 0){
//                 if(out.k[i]==0) fmodes[out.l[i]]=out.mode[i];
//             }
//             if(out.l[i]==3 && out.k[i] == 0){
//                 if(out.k[i]==0) fmodes[out.l[i]]=out.mode[i];
//             }
//         }
//         //test the c0 values
//         for(int i=0; i<out.mode_num; i++){
//             if( 
//                 (out.l[i]==1 && out.k[i]!=1) ||
//                 (out.l[i]==2 && out.k[i]!=0) ||
//                 (out.l[i]==3 && out.k[i]!=0)
//             ) {
//                 TS_ASSERT_LESS_THAN(out.driver->innerproduct(out.mode[i], fmodes[out.l[i]]), 1.e-10);
//             }
//         }
//     }

};
