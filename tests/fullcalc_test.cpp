#include "../src/ThrainMain.h"
#include "../src/ThrainMode.h"
#include "../src/ThrainUnits.h"
#include "../src/ThrainIO.h"
#include "../lib/filelib.h"
#include "doctest.h"

/* INTEGRATION TESTS
* run entire parts of the program and ensure correct answers
* according to known values, low errors, or previous calculation
*/

namespace {
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
      FAIL("could not read in indicated file\n");
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

  std::string const testfilename = "tests/output/fullcalc_test_log.txt";

  struct SetUpTearDown {
    SetUpTearDown() {
      ThrainLogger::setLogLevel(ThrainLogger::LogLevel::ERROR);
      ThrainConfig::reconfigure("tests/tests.config");
    }
    ~SetUpTearDown() {}
  };

} // namespace

TEST_SUITE("FullCalculation [slow]") {
  
/* test uniform density (n=0) polytrope */
TEST_CASE_FIXTURE(SetUpTearDown, "full_calculation_uniform [slow]") {
  printf("\nTEST CALCULATION UNIFORM STAR");
  ThrainLogger::setLogLevel(ThrainLogger::LogLevel::ERROR);
  ThrainConfig::reconfigure("tests/tests.config");

  Calculation::InputData in = make_input_data_pmodes("uniform");
  Calculation::OutputData out;
  CHECK_EQ(0, io::setup_output(in, out));
  // create the star, and print it out
  CHECK_EQ(0, create_star(out));
  io::write_stellar_output(out);
  // ensure the output was written
  CHECK(filelib::exists(ThrainConfig::summaryFileName(out.calcname)));

  // create the modes and print them out
  CHECK_EQ(0, mode::create_modes(out));
  io::write_mode_output(out);
  // ensure the modes were written
  CHECK(filelib::exists(ThrainConfig::calculationFileName(out.calcname, "modes", "mode_1.1.txt")));

  CHECK_LT(out.star_SSR, 1.e-12);
  for(int i=0; i<out.mode_num; i++){
    printf("%d,%d\t", out.l[i], out.k[i]);
    if(out.k[i]!=0){
      printf("%le %le", out.err[0][i], out.err[1][i]);
      CHECK_LT(out.err[0][i], 1.e-8);
      CHECK_LT(out.err[1][i], 1.e-6);
    }
    printf("\n");
  }
  // cleanup
  ThrainLogger::setLogLevel(ThrainLogger::LogLevel::INFO);
  filelib::remove(ThrainConfig::calculationDir(in.calcname));
}

/* test n=1 polytrope */
TEST_CASE_FIXTURE(SetUpTearDown, "full_calculation_polytrope_15 [slow]") {
  printf("TEST CALCULATION POLYTROPE n=1.5\n");
  Calculation::InputData in = make_input_data_pmodes("../tests/poly1.5");
  in.input_params[0] = 1.5; // change polytrope index
  Calculation::OutputData out;
  CHECK_EQ(0, io::setup_output(in, out));
  CHECK_EQ(0, create_star(out));
  CHECK_EQ(0, mode::create_modes(out));

  CHECK_LT(out.star_SSR, 1.e-12);
  for(int i=0; i<out.mode_num; i++){
    printf("%d,%d\t", out.l[i], out.k[i]);
    if(out.k[i]!=0){
      printf("%le %le", out.err[0][i], out.err[1][i]);
      CHECK_LT(out.err[0][i], 1.e-8);
      CHECK_LT(out.err[1][i], 2.e-4);
    }
    printf("\n");
  }
  ThrainLogger::setLogLevel(ThrainLogger::LogLevel::INFO);
}

/* test n=3 polytrope */
TEST_CASE_FIXTURE(SetUpTearDown, "full_calculation_polytrope_30_pmode [slow]") {
  printf("TEST CALCULATION POLYTROPE n=3.0\n");
  Calculation::InputData in = make_input_data_pmodes("../tests/poly3.0");
  in.input_params[0] = 3.0; // change polytrope index
  Calculation::OutputData out;
  CHECK_EQ(0, io::setup_output(in, out));
  CHECK_EQ(0, create_star(out));
  CHECK_EQ(0, mode::create_modes(out));

  std::map<int,ModeBase*> fmodes;

  CHECK_LT(out.star_SSR, 1.e-12);
  for(int i=0; i<out.mode_num; i++){
    printf("%d,%d\t", out.l[i], out.k[i]);
    printf("%le %le\n", out.err[0][i], out.err[1][i]);
    CHECK_LT(out.err[0][i], 1.e-8);
    if(out.k[i]!=0){
      CHECK_LT(out.err[1][i], 2.e-4);
    }

    if(out.l[i]==1 && out.k[i]==1){
      fmodes[out.l[i]] = out.mode[i];
    }
    if(out.l[i]==2 && out.k[i] == 0){
      if(out.k[i]==0) fmodes[out.l[i]]=out.mode[i];
    }
    if(out.l[i]==3 && out.k[i] == 0){
      if(out.k[i]==0) fmodes[out.l[i]]=out.mode[i];
    }
  }
  //test the c0 values
  for(int i=0; i<out.mode_num; i++){
    if( 
      (out.l[i]==1 && out.k[i]!=1) ||
      (out.l[i]==2 && out.k[i]!=0) ||
      (out.l[i]==3 && out.k[i]!=0)
    ) {
      CHECK_LT(out.driver->innerproduct(out.mode[i], fmodes[out.l[i]]), 1.e-10);
    }
  }
}

// TEST_CASE_FIXTURE(SetUpTearDown, "full_calculation_polytrope_30_gmode [slow]") {
//   printf("TEST CALCULATION POLYTROPE n=3.0 GMODE\n");
//   Calculation::InputData in = make_input_data_gmodes("../tests/poly3.0g");
//   in.input_params[0] = 3.0; // change polytrope index
//   Calculation::OutputData out;
//   CHECK_EQ(0, io::setup_output(in, out));
//   printf("before star calc\n"); fflush(stdout);
//   CHECK_EQ(0, create_star(out));
//   printf("after star calc\n"); fflush(stdout);
//   CHECK_EQ(0, mode::create_modes(out));
//   printf("after mode calc\n"); fflush(stdout);

//   std::map<int,ModeBase*> fmodes;

//   CHECK_LT(out.star_SSR, 1.e-12);
//   for(int i=0; i<out.mode_num; i++){
//     printf("%d,%d\t", out.l[i], out.k[i]);
//     printf("%le %le\n", out.err[0][i], out.err[1][i]);
//     fflush(stdout);
//     CHECK_LT(out.err[0][i], 1.e-8);

//     if(out.l[i]==1 && out.k[i]==1){
//       fmodes[out.l[i]] = out.mode[i];
//     }
//     if(out.l[i]==2 && out.k[i] == 0){
//       if(out.k[i]==0) fmodes[out.l[i]]=out.mode[i];
//     }
//     if(out.l[i]==3 && out.k[i] == 0){
//       if(out.k[i]==0) fmodes[out.l[i]]=out.mode[i];
//     }
//   }
//   //test the c0 values
//   for(int i=0; i<out.mode_num; i++){
//     if( 
//       (out.l[i]==1 && out.k[i]!=1) ||
//       (out.l[i]==2 && out.k[i]!=0) ||
//       (out.l[i]==3 && out.k[i]!=0)
//     ) {
//       CHECK_LT(out.driver->innerproduct(out.mode[i], fmodes[out.l[i]]), 1.e-10);
//     }
//   }
// }

/* test CHWD with uniform mu */
// TEST_CASE_FIXTURE(SetUpTearDown, "full_calculation_CHWD_simple [slow]") {
//   printf("TEST CALCULATION CHWD SIMPLE\n");
//   Calculation::InputData in = make_input_data_pmodes("../tests/chwd_0");
//   in.model = model::CHWD;
//   double temp = in.input_params[1];
//   in.input_params = {1.581, 0, temp};
//   in.adiabatic_index = 0;
//   Calculation::OutputData out;
//   CHECK_EQ(0, io::setup_output(in, out));
//   CHECK_EQ(0, create_star(out));
//   CHECK_EQ(0, mode::create_modes(out));

//   std::map<int,ModeBase*> fmodes;

//   CHECK_LT(out.star_SSR, 1.e-12);
//   for(int i=0; i<out.mode_num ; i++){
//     printf("%d,%d\t", out.l[i], out.k[i]);
//     printf("%le\n", out.err[0][i]);
//     fflush(stdout);
//     CHECK_LT(out.err[0][i], 1.e-8);

//     if(out.l[i]==1 && out.k[i]==1){
//       fmodes[out.l[i]] = out.mode[i];
//     }
//     if(out.l[i]==2 && out.k[i] == 0){
//       if(out.k[i]==0) fmodes[out.l[i]]=out.mode[i];
//     }
//     if(out.l[i]==3 && out.k[i] == 0){
//       if(out.k[i]==0) fmodes[out.l[i]]=out.mode[i];
//     }
//   }
//   //test the c0 values
//   for(int i=0; i<out.mode_num; i++){
//     if( 
//       (out.l[i]==1 && out.k[i]!=1) ||
//       (out.l[i]==2 && out.k[i]!=0) ||
//       (out.l[i]==3 && out.k[i]!=0)
//     ) {
//       CHECK_LT(out.driver->innerproduct(out.mode[i], fmodes[out.l[i]]), 1.e-10);
//     }
//   }
// }

/* test CHWD with sigmoidal mu -- pmode */
// TEST_CASE_FIXTURE(SetUpTearDown, "full_calculation_sigmoidal_CHWD_pmode [slow]") {
//   printf("TEST CALCULATION CHWD SIGMOIDAL\n");
//   Calculation::InputData in = make_input_data_pmodes("../tests/chwd_1_p");
//   in.model = model::CHWD;
//   double temp = in.input_params[1];
//   in.input_params = {1.581, 1, temp};
//   Calculation::OutputData out;
//   CHECK_EQ(0, io::setup_output(in, out));
//   CHECK_EQ(0, create_star(out));
//   CHECK_EQ(0, mode::create_modes(out));
//   std::map<int,ModeBase*> fmodes;

//   CHECK_LT(out.star_SSR, 1.e-02);
//   for(int i=0; i<out.mode_num; i++){
//     printf("%d,%d\t", out.l[i], out.k[i]);
//     printf("%le\n", out.err[0][i]);
//     CHECK_LT(out.err[0][i], 1.e-8);
//     if(out.l[i]==1 && out.k[i]==1){
//       fmodes[out.l[i]] = out.mode[i];
//     }
//     if(out.l[i]==2 && out.k[i] == 0){
//       if(out.k[i]==0) fmodes[out.l[i]]=out.mode[i];
//     }
//     if(out.l[i]==3 && out.k[i] == 0){
//       if(out.k[i]==0) fmodes[out.l[i]]=out.mode[i];
//     }
//   }
//   //test the c0 values
//   for(int i=0; i<out.mode_num; i++){
//     if( 
//       (out.l[i]==1 && out.k[i]!=1) ||
//       (out.l[i]==2 && out.k[i]!=0) ||
//       (out.l[i]==3 && out.k[i]!=0)
//     ) {
//       CHECK_LT(out.driver->innerproduct(out.mode[i], fmodes[out.l[i]]), 1.e-10);
//     }
//   }
// }

// /* test CHWD with sigmoidal mu -- gmode */
// TEST_CASE_FIXTURE(SetUpTearDown, "full_calculation_sigmoidal_CHWD_gmode [slow]") {
//   printf("TEST CALCULATION CHWD SIGMOIDAL GMODE\n");
//   Calculation::InputData in = make_input_data_gmodes("../tests/chwd_1_g");
//   in.model = model::CHWD;
//   double temp = in.input_params[1];
//   in.input_params = {1.581, 1, temp};
//   Calculation::OutputData out;
//   CHECK_EQ(0, io::setup_output(in, out));
//   CHECK_EQ(0, create_star(out));
//   CHECK_EQ(0, mode::create_modes(out));

//   std::map<int,ModeBase*> fmodes;

//   CHECK_LT(out.star_SSR, 1.e-02);
//   for(int i=0; i<out.mode_num; i++){
//     printf("%d,%d\t", out.l[i], out.k[i]);
//     printf("%le\n", out.err[0][i]);
//     CHECK_LT(out.err[0][i], 1.e-8);

//     if(out.l[i]==1 && out.k[i]==1){
//       fmodes[out.l[i]] = out.mode[i];
//     }
//     if(out.l[i]==2 && out.k[i] == 0){
//       if(out.k[i]==0) fmodes[out.l[i]]=out.mode[i];
//     }
//     if(out.l[i]==3 && out.k[i] == 0){
//       if(out.k[i]==0) fmodes[out.l[i]]=out.mode[i];
//     }
//   }
//   //test the c0 values
//   for(int i=0; i<out.mode_num; i++){
//     if( 
//       (out.l[i]==1 && out.k[i]!=1) ||
//       (out.l[i]==2 && out.k[i]!=0) ||
//       (out.l[i]==3 && out.k[i]!=0)
//     ) {
//       CHECK_LT(out.driver->innerproduct(out.mode[i], fmodes[out.l[i]]), 1.e-10);
//     }
//   }
// }
}
