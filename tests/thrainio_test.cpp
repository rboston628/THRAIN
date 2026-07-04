#include "../src/ThrainMain.h"
#include "../src/ThrainIO.h"
#include "doctest.h"

namespace {
  void make_test_input(std::string filename, std::string towrite){
    FILE* outfile = fopen(ThrainConfig::inputFileName(filename).c_str(), "w");
    REQUIRE(outfile);
    fprintf(outfile, "%s", towrite.c_str());
    fclose(outfile);
  }

  std::string remove_all_whitespace(std::string const& s){
    std::string result = s;
    result.erase(std::remove_if(result.begin(), result.end(), [](unsigned char x) {
        return std::isspace(x);
    }), result.end());
    return result;
  }

  void read_entire_file(std::string filename, std::string& contents){
    FILE* infile = fopen(ThrainConfig::echoedFileName(filename).c_str(), "r");
    REQUIRE(infile);
    fseek(infile, 0, SEEK_END);
    std::size_t fsize = ftell(infile);
    fseek(infile, 0, SEEK_SET);

    char *read_buffer = new char[fsize + 1];
    fread(read_buffer, fsize, 1, infile);
    read_buffer[fsize] = 0;
    fclose(infile);
    contents = std::string(read_buffer);
    // remove all whitespace from the string for comparison
    contents = remove_all_whitespace(contents);
    delete[] read_buffer;
  }
}

TEST_SUITE("ThrainIO [unit]") {

  /* test basic calculation setup */

  TEST_CASE("io: fail_bad_file") {
    printf("IO TEST - BAD FILE\n");
    Calculation::InputData data;
    std::string const badfilename = "nonexistent.txt";
    CHECK(!filelib::exists(ThrainConfig::inputFileName(badfilename)));
    CHECK_EQ(1, io::read_input(badfilename, data));
  }

  TEST_CASE("io: fail_bad_calcname") {
    printf("IO TEST - BAD CALCNAME\n");
    Calculation::InputData data;
    std::string const testfilename = "test_file.txt";
    std::string const filecontents = "Nombre: espanol\n";
    make_test_input(testfilename, filecontents);
    FILE* in = fopen(ThrainConfig::inputFileName(testfilename).c_str(), "r");
    CHECK_EQ(1, io::read_calcname(in, data));
    fclose(in);
    // also ensure read_input fails when these do
    CHECK_EQ(1, io::read_input(testfilename, data));
  }

   TEST_CASE("io: fail_no_calcname") {
    printf("IO TEST - NO CALCNAME\n");
    Calculation::InputData data;
    std::string const testfilename = "test_file.txt";
    std::string const filecontents = "Name:\n";
    make_test_input(testfilename, filecontents);
    FILE* in = fopen(ThrainConfig::inputFileName(testfilename).c_str(), "r");
    CHECK_EQ(1, io::read_calcname(in, data));
    fclose(in);
    // also ensure read_input fails when these do
    CHECK_EQ(1, io::read_input(testfilename, data));
  }

  TEST_CASE("io: read_calcname") {
    printf("IO TEST - READ CALCNAME\n");
    Calculation::InputData data;
    std::string const testfilename = "test_file.txt";
    std::string const filecontents = "Name: valid_test_name\n";
    make_test_input(testfilename, filecontents);
    FILE* in = fopen(ThrainConfig::inputFileName(testfilename).c_str(), "r");
    CHECK_EQ(0, io::read_calcname(in, data));
    CHECK_EQ(data.calcname, "valid_test_name");
    fclose(in);
  }

  /* tests of reading models */

  TEST_CASE("io: fail_bad_model") {
    printf("IO TEST - BAD MODEL\n");
    Calculation::InputData data;
    std::string const testfilename = "test_file.txt";
    std::string const filecontents = 
      "Name: valid_test_name\n"
      "Model: GR fakemodel\n";
    make_test_input(testfilename, filecontents);
    FILE* in = fopen(ThrainConfig::inputFileName(testfilename).c_str(), "r");
    CHECK_EQ(1, io::read_model(in, data));
    CHECK_EQ(data.regime, regime::PN0);
    fclose(in);
    // also ensure read_input fails when these do
    CHECK_EQ(1, io::read_input(testfilename, data));
  }

  TEST_CASE("io: fail_no_model") {
    printf("IO TEST - NO MODEL\n");
    Calculation::InputData data;
    std::string const testfilename = "test_file.txt";
    std::string const filecontents = 
      "Name: valid_test_name\n"
      "Model: newtonian\n";
    make_test_input(testfilename, filecontents);
    FILE* in = fopen(ThrainConfig::inputFileName(testfilename).c_str(), "r");
    CHECK_EQ(1, io::read_model(in, data));
    CHECK_EQ(data.regime, regime::PN0);
    fclose(in);
    // also ensure read_input fails when these do
    CHECK_EQ(1, io::read_input(testfilename, data));
  }

  void do_test_read_model(model::StellarModel m){
    Calculation::InputData data;
    std::string testfilename("test_file.txt");
    std::string filecontents = "Model: newtonian ";
    switch(m){
    case model::polytrope:
      filecontents += "polytrope\t\n";
      break;
    case model::CHWD:
      filecontents += "CHWD\t\n";
      break;
    case model::MESA:
      filecontents += "MESA\t\n";
      break;
    case model::SWD:
      filecontents += "SWD\t\n";
      break;
    }
    make_test_input(testfilename, filecontents);
    FILE* in = fopen(ThrainConfig::inputFileName(testfilename).c_str(), "r");
    CHECK_EQ(0, io::read_model(in, data));
    CHECK_EQ(data.regime, regime::PN0);
    CHECK_EQ(data.model, m);
    fclose(in);
  }

  void test_real_all_models(){
    printf("IO TEST - TEST ALL MODELS\n");
    do_test_read_model(model::polytrope);
    do_test_read_model(model::CHWD);
    do_test_read_model(model::MESA);
    do_test_read_model(model::SWD);
  }

  /* tests of reading the units */

  TEST_CASE("io: fail_no_units") {
    printf("IO TEST - NO UNITS\n");
    Calculation::InputData data;
    std::string const testfilename = "test_file.txt";
    std::string const filecontents = "Unidades: buenos\n";
    make_test_input(testfilename, filecontents);
    FILE* in = fopen(ThrainConfig::inputFileName(testfilename).c_str(), "r");
    CHECK_EQ(1, io::read_units(in, data));
    fclose(in);
  }

  TEST_CASE("io: fail_bad_units") {
    printf("IO TEST - BAD UNITS\n");
    Calculation::InputData data;
    std::string const testfilename = "test_file.txt";
    std::string const filecontents = "Units: fakeunits\n";
    make_test_input(testfilename, filecontents);
    FILE* in = fopen(ThrainConfig::inputFileName(testfilename).c_str(), "r");
    CHECK_EQ(1, io::read_units(in, data));
    fclose(in);
  }

  void do_test_units_types(units::Units u) {
    Calculation::InputData data;
    std::string const testfilename = "test_file.txt";
    std::string filecontents = "Units: ";
    switch(u){
    case units::Units::astro:
      filecontents += "astro\n";
      break;
    case units::Units::geo:
      filecontents += "geo\n";
      break;
    case units::Units::SI:
      filecontents += "SI\n";
      break;
    case units::Units::CGS:
      filecontents += "CGS\n";
      break;
    }
    make_test_input(testfilename,filecontents);
    FILE* in = fopen(ThrainConfig::inputFileName(testfilename).c_str(), "r");
    CHECK_EQ(0, io::read_units(in, data));
    CHECK_EQ(data.units, u);
    fclose(in);
  }

  TEST_CASE("io: read_all_units") {
    printf("IO TEST - ALL UNITS\n");
    do_test_units_types(units::Units::astro);
    do_test_units_types(units::Units::geo);
    do_test_units_types(units::Units::SI);
    do_test_units_types(units::Units::CGS);
  }

  /* tests of reading the modes */

  TEST_CASE("io: no_modetype") {
    printf("IO TEST - NO MODETYPE\n");
    Calculation::InputData data;
    std::string const testfilename = "test_file.txt";
    std::string const filecontents = "Frequencies: 5/3\n"; //will become nonradial mode
    make_test_input(testfilename, filecontents);
    FILE* in = fopen(ThrainConfig::inputFileName(testfilename).c_str(), "r");
    CHECK_EQ(1, io::read_frequencies(in, data));
    CHECK_EQ(data.modetype, modetype::nonradial);
    // CHECK_ISNAN( data.adiabatic_index );
    CHECK(std::isnan(data.adiabatic_index));
    fclose(in);
  }

  TEST_CASE("io: read_bad_modetype") {
    printf("IO TEST - BAD MODETYPE\n");
    Calculation::InputData data;
    std::string const testfilename = "test_file.txt";
    std::string const filecontents = "Frequencies: fakemode 5/3\n"; //will become nonradial mode
    make_test_input(testfilename, filecontents);
    FILE* in = fopen(ThrainConfig::inputFileName(testfilename).c_str(), "r");
    CHECK_EQ(0, io::read_frequencies(in, data));
    CHECK_EQ(data.modetype, modetype::nonradial);
    CHECK_EQ(data.adiabatic_index, 5./3.);
    CHECK_EQ(data.mode_num, 0);
    fclose(in);
  }

  void do_test_modetype(modetype::ModeType mtype) {
    Calculation::InputData data;
    std::string const testfilename = "test_file.txt"; 
    std::string filecontents = "Frequencies: ";
    switch(mtype){
    case modetype::radial:
      filecontents += "radial";
      break;
    case modetype::cowling:
      filecontents += "cowling";
      break;
    case modetype::nonradial:
      filecontents += "nonradial";
      break;
    } 
    filecontents += " 5/3\n";
    make_test_input(testfilename,filecontents);
    FILE* in = fopen(ThrainConfig::inputFileName(testfilename).c_str(), "r");
    CHECK_EQ(0, io::read_frequencies(in, data));
    CHECK_EQ(data.modetype, mtype);
    CHECK_EQ(data.adiabatic_index, 5./3.);
    fclose(in);
  }

  TEST_CASE("io: read_all_modetypes") {
    printf("IO TEST - ALL MODETYPES\n");
    do_test_modetype(modetype::radial);
    do_test_modetype(modetype::cowling);
    do_test_modetype(modetype::nonradial);
  }

  void do_test_adiabatic_index(std::string index, double res) {
    Calculation::InputData data;
    std::string const testfilename = "test_file.txt";
    std::string const filecontents = "Frequencies: nonradial " + index;
    make_test_input(testfilename, filecontents);
    FILE* in = fopen(ThrainConfig::inputFileName(testfilename).c_str(), "r");
    CHECK_EQ(0, io::read_frequencies(in, data));
    CHECK_EQ(data.adiabatic_index, res);
    fclose(in);
  }

  void test_read_adiabatic_index(){
    printf("IO TEST - ALL ADIABATIC INDICES\n");
    do_test_adiabatic_index("5/3", 5./3.);
    do_test_adiabatic_index("4/3", 4./3.);
    do_test_adiabatic_index("5/2", 2.5);
    do_test_adiabatic_index("3.0", 3.0);
    do_test_adiabatic_index("2.0", 2.0);
    do_test_adiabatic_index("1.5", 1.5);
    do_test_adiabatic_index("5.0", 5.0);
    do_test_adiabatic_index("0", 0.0);
    do_test_adiabatic_index("1", 1.0);
    do_test_adiabatic_index("2", 2.0);
  }

  TEST_CASE("io: read_mode_list") {
    printf("IO TEST - READ MODE LIST\n");
    Calculation::InputData data;
    std::string const testfilename = "test_file.txt";
    std::string filecontents = "Frequencies: nonradial 5/3\n";
    int Lvalues[3] = {1,2,3};
    int Kvalues[5] = {1,2,3,4,5};
    for(int L : Lvalues){
      for(int K : Kvalues){
      filecontents += strmakef("%d,%d\n", L,K);
      }
    }
    make_test_input(testfilename, filecontents);
    FILE* in = fopen(ThrainConfig::inputFileName(testfilename).c_str(), "r");
    CHECK_EQ(0, io::read_frequencies(in, data));
    CHECK_EQ(data.mode_num, 17); //L=2,3 will have K=0 added
    CHECK_EQ(3, data.l.size());
    CHECK_EQ(3, data.kl.size());
  
    for(int i=0; i<3; i++){
      for(int j=0; j<5; j++){
      CHECK_EQ(data.l.count(Lvalues[i]), 1);
      CHECK_LE(5, data.kl.at(Lvalues[i]).size());
      CHECK_LE(data.kl.at(Lvalues[i]).size(),6);
      CHECK_EQ(data.kl.at(Lvalues[i]).count(Kvalues[j]), 1);
      }
    }
    fclose(in);
  }

  TEST_CASE("io: read_mode_list_unordered") {
    printf("IO TEST - READ UNORDERED MODE LIST\n");
    // tests that duplicate L,K will be removed
    // and that they will be put in order
    Calculation::InputData data;
    std::string const testfilename = "test_file.txt";
    std::string filecontents = "Frequencies: nonradial 5/3\n";
    int Lvalues[] = {3,1,2   /* duplicates to be ignored */, 3,1,2};
    int Kvalues[] = {4,2,3,5,1 /* duplicates to be ignored */, 5,5,1};
    for(int L : Lvalues){
      for(int K : Kvalues){
      filecontents += strmakef("%d,%d\n", L,K);
      }
    }
    make_test_input(testfilename, filecontents);
    FILE* in = fopen(ThrainConfig::inputFileName(testfilename).c_str(), "r");
    CHECK_EQ(0, io::read_frequencies(in, data));
    CHECK_EQ(data.mode_num, 17); //L=2,3 will have K=0 added
    CHECK_EQ(3, data.l.size());
    CHECK_EQ(3, data.kl.size());
    
    int Lout[17] = {1,1,1,1,1, 2,2,2,2,2,2, 3,3,3,3,3,3};
    int Kout[17] = {1,2,3,4,5, 0,1,2,3,4,5, 0,1,2,3,4,5};
    std::vector<int> Ldat, Kdat;
    for(auto it=data.l.begin(); it!=data.l.end(); it++){
      auto kforl = data.kl[*it];
      for(auto jt=kforl.begin(); jt!=kforl.end(); jt++){
      Ldat.push_back(*it);
      Kdat.push_back(*jt);
      }
    }
    for(int i=0; i<17; i++){
      CHECK_EQ(Lout[i], Ldat[i]);
      CHECK_EQ(Kout[i], Kdat[i]);
    }
    fclose(in);
  }

  TEST_CASE("io: read_mode_list_no_dipole_f") {
    printf("IO TEST - REMOVE L=1 F-MODE\n");
    // tests that duplicate L,K will be removed
    // and that they will be put in order
    Calculation::InputData data;
    std::string const testfilename = "test_file.txt";
    std::string filecontents = "Frequencies: nonradial 5/3\n";
    int Lvalues[] = {1,2};
    int Kvalues[] = {0,1,2,3}; //the k=0 will be dropped when l=1
    for(int L : Lvalues){
      for(int K : Kvalues){
      filecontents += strmakef("%d,%d\n", L,K);
      }
    }
      make_test_input(testfilename, filecontents);
    FILE* in = fopen(ThrainConfig::inputFileName(testfilename).c_str(), "r");
    CHECK_EQ(0, io::read_frequencies(in, data));
    
    // there is no l=1,k=0 mode present
    int Lout[7] = {1,1,1, 2,2,2,2}; 
    int Kout[7] = {1,2,3, 0,1,2,3};
    std::vector<int> Ldat, Kdat;
    for(auto it=data.l.begin(); it!=data.l.end(); it++){
      auto kforl = data.kl[*it];
      for(auto jt=kforl.begin(); jt!=kforl.end(); jt++){
      Ldat.push_back(*it);
      Kdat.push_back(*jt);
      }
    }
    for(int i=0; i<7; i++){
      CHECK_EQ(Lout[i], Ldat[i]);
      CHECK_EQ(Kout[i], Kdat[i]);
    }
    fclose(in);
  }

  TEST_CASE("io: read_bad_mode_list") {
    printf("IO TEST - READ BAD MODE LIST\n");
    Calculation::InputData data;
    std::string const testfilename = "test_file.txt";
    std::string const filecontents = 
      "Frequencies: nonradial 5/3\n"
      "1,1\n"
      "1,2\n"
      "a,b\n"
      "1,3\n";
    make_test_input(testfilename, filecontents);
    FILE* in = fopen(ThrainConfig::inputFileName(testfilename).c_str(), "r");
    CHECK_EQ(0, io::read_frequencies(in, data));
    CHECK_EQ(data.mode_num, 3);
    CHECK_EQ(data.l.size(), 1);
    CHECK_EQ(data.kl[1].size(), 3);
    std::set<int> kforl = data.kl[1];
    std::vector<int> k(kforl.begin(), kforl.end());
    for(int i=0; i<3; i++){
      CHECK_EQ(k[i], i+1); // 1, 2, 3
    }
    fclose(in);
  }

  /* test reading polytrope input data */

  TEST_CASE("io: fail_bad_params") {
    printf("IO TEST - BAD PARAMS\n");
    Calculation::InputData data;
    std::string const testfilename = "test_file.txt";
    std::string const filecontents = 
      "# test a comment\n\n"
      "Name: valid_test_name\n"
      "Model: newtonian polytrope 1.5 1\n"
      "Params: fakeparam1 1.0 fakeparm2 10.0\n";
    make_test_input(testfilename, filecontents);
    CHECK_EQ(1, io::read_input(testfilename, data));
    CHECK_EQ(data.calcname, "valid_test_name");
    CHECK_EQ(data.regime, regime::PN0); 
    CHECK_EQ(data.model, model::polytrope);
    CHECK_EQ(data.input_params.size(), 2);
    CHECK_EQ(data.input_params[0], 1.5);
    CHECK_EQ(data.input_params[1], 1);
    CHECK_EQ(data.Ngrid, 1);
  }

  TEST_CASE("io: open_file_polytrope") {
    printf("IO TEST - OPEN FILE\n");
    Calculation::InputData data;
    std::string testfilename("test_file.txt");
    std::string filecontents = 
      "# test a comment\n\n"
      "Name: valid_test_name\n"
      "Model: newtonian polytrope 1.5 1\n"
      "Params: mass 1.0 radius 10.0\n"
      "Units: geo\n"
      "Frequencies: cowling 4/3\n\n";
    make_test_input(testfilename,filecontents);
    CHECK_EQ(0, io::read_input(testfilename, data));
    CHECK_EQ(data.calcname, "valid_test_name");
    CHECK_EQ(data.regime, regime::PN0); 
    CHECK_EQ(data.model, model::polytrope);
    CHECK_EQ(data.input_params.size(), 2);
    CHECK_EQ(data.input_params[0], 1.5);
    CHECK_EQ(data.input_params[1], 1);
    CHECK_EQ(data.Ngrid, 1);
    CHECK_EQ(data.radius, 10.0);
    CHECK_EQ(data.mass, 1.0);
    CHECK_EQ(data.units, units::Units::geo);
    CHECK_EQ(data.modetype, modetype::cowling);
    CHECK_EQ(data.adiabatic_index, 4./3.);
    CHECK_EQ(data.mode_num, 0);
  }

  TEST_CASE("io: open_file_polytrope_RZ") {
    Calculation::InputData data;
    std::string testfilename("test_file.txt");
    std::string filecontents = 
      "# test a comment\n\n"
      "Name: valid_test_name\n"
      "Model: newtonian polytrope 1.5 1\n"
      "Params: zsurf 0.01 radius 10.0\n"
      "Units: geo\n"
      "Frequencies: cowling 4/3\n\n";
    make_test_input(testfilename,filecontents);
    CHECK_EQ(0, io::read_input(testfilename, data));
    CHECK_EQ(data.radius, 10.0);
    CHECK_EQ(data.zsurf, 0.01);
  }

  TEST_CASE("io: open_file_polytrope_MZ") {
    Calculation::InputData data;
    std::string testfilename("test_file.txt");
    std::string filecontents = 
      "# test a comment\n\n"
      "Name: valid_test_name\n"
      "Model: newtonian polytrope 1.5 1\n"
      "Params: zsurf 0.01 mass 1.0\n"
      "Units: geo\n"
      "Frequencies: cowling 4/3\n\n";
    make_test_input(testfilename,filecontents);
    CHECK_EQ(0, io::read_input(testfilename, data));
    CHECK_EQ(data.mass, 1.0);
    CHECK_EQ(data.zsurf, 0.01);
  }

  TEST_CASE("io: open_file_polytrope_RLogg") {
    Calculation::InputData data;
    std::string testfilename("test_file.txt");
    std::string filecontents = 
      "# test a comment\n\n"
      "Name: valid_test_name\n"
      "Model: newtonian polytrope 1.5 1\n"
      "Params: logg 10.0 radius 10.0\n"
      "Units: geo\n"
      "Frequencies: cowling 4/3\n\n";
    make_test_input(testfilename,filecontents);
    CHECK_EQ(0, io::read_input(testfilename, data));
    CHECK_EQ(data.radius, 10.0);
    CHECK_EQ(data.logg, 10.0);
  }

  TEST_CASE("io: open_file_polytrope_MLogg") {
    Calculation::InputData data;
    std::string testfilename("test_file.txt");
    std::string filecontents = 
      "# test a comment\n\n"
      "Name: valid_test_name\n"
      "Model: newtonian polytrope 1.5 1\n"
      "Params: logg 10.0 mass 1.0\n"
      "Units: geo\n"
      "Frequencies: cowling 4/3\n\n";
    make_test_input(testfilename,filecontents);
    CHECK_EQ(0, io::read_input(testfilename, data));
    CHECK_EQ(data.mass, 1.0);
    CHECK_EQ(data.logg, 10.0);
  }

  TEST_CASE("io: open_file_polytrope_fail_RR") {
    Calculation::InputData data;
    std::string testfilename("test_file.txt");
    std::string filecontents = 
      "# test a comment\n\n"
      "Name: valid_test_name\n"
      "Model: newtonian polytrope 1.5 1\n"
      "Params: radius 10.0 radius 10.0\n"
      "Units: geo\n"
      "Frequencies: cowling 4/3\n\n";
    make_test_input(testfilename,filecontents);
    CHECK_EQ(1, io::read_input(testfilename, data));
  }

  TEST_CASE("io: open_file_polytrope_fail_teff") {
    Calculation::InputData data;
    std::string testfilename("test_file.txt");
    std::string filecontents = 
      "# test a comment\n\n"
      "Name: valid_test_name\n"
      "Model: newtonian polytrope 1.5 1\n"
      "Params: radius 10.0 teff 10000.0\n"
      "Units: geo\n"
      "Frequencies: cowling 4/3\n\n";
    make_test_input(testfilename,filecontents);
    CHECK_EQ(1, io::read_input(testfilename, data));
  }

  /* test reading CHWD input  */
  /* test reading MESA input  */
  /* test reading SWD input  */

  /* test echo */
  
  TEST_CASE("io: echo_input") {
    // create a test file, echo it, read in the echo, check match
    Calculation::InputData data;
    std::string testfilename("test_file.txt");
    std::string filecontents = 
      "Name:\ttest_file\n"
      "Model:\tnewtonian polytrope 1.50 1\n"
      "Params:\tmass 1\tradius 10\n"
      "Units:\tgeo\n\n"
      "Frequencies:\tcowling\t1.500\n";
    make_test_input(testfilename,filecontents);
    CHECK_EQ(0, io::read_input(testfilename, data));
    std::string echoedcontents;
    std::string readfilename("test_file");
    CHECK_EQ(0, io::echo_input(data));
    CAPTURE(ThrainConfig::echoedFileName(readfilename));
    read_entire_file(readfilename, echoedcontents);
    // remove whitespace from the original file contents for comparison
    filecontents = remove_all_whitespace(filecontents);
    CAPTURE(filecontents);
    CAPTURE(echoedcontents);
    CHECK_EQ(filecontents, echoedcontents);
    // clean up the file tree
    filelib::remove(ThrainConfig::inputFileName(testfilename));
    filelib::remove(ThrainConfig::echoedFileName(readfilename));
  }

  /* test setup output, which created output object */

  TEST_CASE("io: setup_output") {
    Calculation::InputData indata;
    std::string testfilename("test_file.txt");
    std::string filecontents =
      "# test a comment\n\n"
      "Name: TEST_SETUP_OUTPUT\n"
      "Model: newtonian polytrope 1.5 1\n"
      "Params: mass 1.0 radius 10.0\n"
      "Units: geo\n"
      "Frequencies: nonradial 5/3\n\n";
    int Lvalues[] = {2,1,2}; // out of order and includes duplicates
    int Kvalues[] = {3,1,2,3}; // out of order and includes duplicates
    char freqLine[10] = "";
    for(int L : Lvalues){
      for(int K : Kvalues){
      snprintf(freqLine,10, "%d,%d\n", L,K);
      filecontents += std::string(freqLine);
      }
    }
    make_test_input(testfilename,filecontents);
    CHECK_EQ(0, io::read_input(testfilename, indata));
    CHECK_EQ(indata.mode_num, 7);

    Calculation::OutputData outdata;
    outdata.star = nullptr;
    outdata.driver = nullptr;
    CHECK_EQ(0, io::setup_output(indata, outdata));
    CHECK_EQ(outdata.calcname, "TEST_SETUP_OUTPUT");
    CHECK_EQ(outdata.regime, regime::PN0); 
    CHECK_EQ(outdata.model, model::polytrope);
    CHECK_EQ(outdata.modetype, modetype::nonradial);
    CHECK_EQ(outdata.units, units::Units::geo);
    CHECK_EQ(outdata.input_params.size(), 2);
    CHECK_EQ(outdata.input_params[0], 1.5);
    CHECK_EQ(outdata.input_params[1], 1);
    CHECK_EQ(outdata.Ngrid, 1);
    CHECK_EQ(outdata.radius, 10.0);
    CHECK_EQ(outdata.mass, 1.0);
    CHECK_EQ(outdata.zsurf, units::getZsurfFromRM(outdata));
    CHECK_EQ(outdata.logg, units::getLoggFromRM(outdata));
    CHECK_EQ(outdata.mode_num, 7);
    CHECK_EQ(outdata.adiabatic_index, 5./3.);
    int Lout[] = {1,1,1, 2,2,2,2}; // in order, no duplicates
    int Kout[] = {1,2,3, 0,1,2,3}; // in order, no duplicates
    for(int i=0; i<7; i++){
      CHECK_EQ(Lout[i], outdata.l[i]);
      CHECK_EQ(Kout[i], outdata.k[i]);
    }
    CHECK_EQ(outdata.w.capacity(), 7);
    CHECK_EQ(outdata.f.capacity(), 7);
    CHECK_EQ(outdata.period.capacity(), 7);
    CHECK_EQ(outdata.mode_SSR.capacity(), 7);
    CHECK_EQ(outdata.i_err, 2);
    CHECK_EQ(outdata.error[0], true);
    CHECK_EQ(outdata.error[1], false);
    CHECK_EQ(outdata.error[2], false);
    CHECK_EQ(outdata.error[3], true);
  }
}
