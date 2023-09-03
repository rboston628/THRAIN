#include "../src/ThrainMain.h"
#include "../src/ThrainIO.h"
#include <cxxtest/TestSuite.h>

class IOBaseTest : public CxxTest::TestSuite {
public:

    void make_test_input(std::string filename, std::string towrite){
        FILE* outfile = fopen(filename.c_str(), "w");
        if(!outfile){
            TS_FAIL("could not open test input\n");
        }
        fprintf(outfile, "%s", towrite.c_str());
        fclose(outfile);
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

    /* test in basic calculation setup */

    void test_fail_bad_file() {
        Calculation::InputData data;
        char badfilename[] = "tests/nonexistent.txt";
        TS_ASSERT_EQUALS(1, io::read_input(badfilename, data));
    }

    void test_fail_bad_calcname() {
        Calculation::InputData data;
        std::string testfilename("tests/test_file.txt");
        std::string filecontents = "Nombre: espanol\n";
        make_test_input(testfilename,filecontents);
        FILE* in = fopen(testfilename.c_str(), "r");
        TS_ASSERT_EQUALS(1, io::read_calcname(in, data));
        fclose(in);
    }

     void test_fail_no_calcname() {
        Calculation::InputData data;
        std::string testfilename("tests/test_file.txt");
        std::string filecontents = "Name:\n";
        make_test_input(testfilename,filecontents);
        FILE* in = fopen(testfilename.c_str(), "r");
        TS_ASSERT_EQUALS(1, io::read_calcname(in, data));
        printf("CALCNAME = %s\n", data.calcname.c_str());
        fclose(in);
    }

    void test_read_calcname() {
        Calculation::InputData data;
        std::string testfilename("tests/test_file.txt");
        std::string filecontents = "Name: valid_test_name\n";
        make_test_input(testfilename,filecontents);
        FILE* in = fopen(testfilename.c_str(), "r");
        TS_ASSERT_EQUALS(0, io::read_calcname(in, data));
        TS_ASSERT_EQUALS(data.calcname, "valid_test_name");
        fclose(in);
    }

    /* tests of reading models */

    void test_fail_bad_model() {
        Calculation::InputData data;
        std::string testfilename("tests/test_file.txt");
        std::string filecontents = "Model: GR fakemodel\n";
        make_test_input(testfilename,filecontents);
        FILE* in = fopen(testfilename.c_str(), "r");
        TS_ASSERT_EQUALS(1, io::read_model(in, data));
        TS_ASSERT_EQUALS(data.regime, regime::PN0);
        fclose(in);
    }

    void test_fail_no_model() {
        Calculation::InputData data;
        std::string testfilename("tests/test_file.txt");
        std::string filecontents = "Model: newtonian\n";
        make_test_input(testfilename,filecontents);
        FILE* in = fopen(testfilename.c_str(), "r");
        TS_ASSERT_EQUALS(1, io::read_model(in, data));
        TS_ASSERT_EQUALS(data.regime, regime::PN0);
        fclose(in);
    }

    void do_test_read_model(model::StellarModel m){
        Calculation::InputData data;
        std::string testfilename("tests/test_file.txt");
        std::string filecontents = "Model: newtonian ";
        switch(m){
        case model::polytrope:
            filecontents += "polytrope\n";
            break;
        case model::CHWD:
            filecontents += "CHWD\n";
            break;
        case model::MESA:
            filecontents += "MESA\n";
            break;
        case model::SWD:
            filecontents += "SWD\n";
            break;
        }
        make_test_input(testfilename,filecontents);
        FILE* in = fopen(testfilename.c_str(), "r");
        TS_ASSERT_EQUALS(0, io::read_model(in, data));
        TS_ASSERT_EQUALS(data.regime, regime::PN0);
        TS_ASSERT_EQUALS(data.model, m);
        fclose(in);
    }

    void test_real_all_models(){
        do_test_read_model(model::polytrope);
        do_test_read_model(model::CHWD);
        do_test_read_model(model::MESA);
        do_test_read_model(model::SWD);
    }

    /* tests of reading the units */

    void test_fail_no_units() {
        Calculation::InputData data;
        std::string testfilename("tests/test_file.txt");
        std::string filecontents = "Unidades: buenos\n";
        make_test_input(testfilename,filecontents);
        FILE* in = fopen(testfilename.c_str(), "r");
        TS_ASSERT_EQUALS(1, io::read_units(in, data));
        fclose(in);
    }

    void test_fail_bad_units() {
        Calculation::InputData data;
        std::string testfilename("tests/test_file.txt");
        std::string filecontents = "Units: fakeunits\n";
        make_test_input(testfilename,filecontents);
        FILE* in = fopen(testfilename.c_str(), "r");
        TS_ASSERT_EQUALS(1, io::read_units(in, data));
        fclose(in);
    }

    void do_test_units_types(units::Units u) {
        Calculation::InputData data;
        std::string testfilename("tests/test_file.txt");
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
        FILE* in = fopen(testfilename.c_str(), "r");
        TS_ASSERT_EQUALS(0, io::read_units(in, data));
        TS_ASSERT_EQUALS(data.units, u);
        fclose(in);
    }

    void test_read_all_units() {
        do_test_units_types(units::Units::astro);
        do_test_units_types(units::Units::geo);
        do_test_units_types(units::Units::SI);
        do_test_units_types(units::Units::CGS);
    }

    /* tests of reading the modes */

    void test_no_modetype() {
        Calculation::InputData data;
        std::string testfilename("tests/test_file.txt");
        std::string filecontents = "Frequencies: 5/3\n"; //will become nonradial mode
        make_test_input(testfilename,filecontents);
        FILE* in = fopen(testfilename.c_str(), "r");
        TS_ASSERT_EQUALS(1, io::read_frequencies(in, data));
        TS_ASSERT_EQUALS(data.modetype, modetype::nonradial);
        TS_ASSERT_IS_NAN(data.adiabatic_index);
        fclose(in);
    }

    void test_read_bad_modetype() {
        Calculation::InputData data;
        std::string testfilename("tests/test_file.txt");
        std::string filecontents = "Frequencies: fakemode 5/3\n"; //will become nonradial mode
        make_test_input(testfilename,filecontents);
        FILE* in = fopen(testfilename.c_str(), "r");
        TS_ASSERT_EQUALS(0, io::read_frequencies(in, data));
        TS_ASSERT_EQUALS(data.modetype, modetype::nonradial);
        TS_ASSERT_EQUALS(data.adiabatic_index, 5./3.);
        TS_ASSERT_EQUALS(data.mode_num, 0);
        fclose(in);
    }

    void do_test_modetype(modetype::ModeType mtype) {
        Calculation::InputData data;
        std::string testfilename("tests/test_file.txt");
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
        FILE* in = fopen(testfilename.c_str(), "r");
        TS_ASSERT_EQUALS(0, io::read_frequencies(in, data));
        TS_ASSERT_EQUALS(data.modetype, mtype);
        TS_ASSERT_EQUALS(data.adiabatic_index, 5./3.);
        fclose(in);
    }

    void test_read_all_modetypes() {
        do_test_modetype(modetype::radial);
        do_test_modetype(modetype::cowling);
        do_test_modetype(modetype::nonradial);
    }

    void do_test_adiabatic_index(std::string index, double res) {
        Calculation::InputData data;
        std::string testfilename("tests/test_file.txt");
        std::string filecontents = "Frequencies: nonradial ";
        filecontents += index;
        make_test_input(testfilename,filecontents);
        FILE* in = fopen(testfilename.c_str(), "r");
        TS_ASSERT_EQUALS(0, io::read_frequencies(in, data));
        TS_ASSERT_EQUALS(data.adiabatic_index, res);
        fclose(in);
    }

    void test_read_adiabatic_index(){
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

    void test_read_mode_list() {
        Calculation::InputData data;
        std::string testfilename("tests/test_file.txt");
        std::string filecontents = "Frequencies: nonradial 5/3\n";
        int Lvalues[3] = {1,2,3};
        int Kvalues[5] = {1,2,3,4,5};
        for(int L : Lvalues){
            for(int K : Kvalues){
                filecontents += strmakef("%d,%d\n", L,K);
            }
        }
        make_test_input(testfilename,filecontents);
        FILE* in = fopen(testfilename.c_str(), "r");
        TS_ASSERT_EQUALS(0, io::read_frequencies(in, data));
        TS_ASSERT_EQUALS(data.mode_num, 17); //L=2,3 will have K=0 added
        TS_ASSERT_EQUALS(3, data.l.size());
        TS_ASSERT_EQUALS(3, data.kl.size());
        
        for(int i=0; i<3; i++){
            for(int j=0; j<5; j++){
                TS_ASSERT_EQUALS(data.l.count(Lvalues[i]), 1);
                TS_ASSERT_LESS_THAN_EQUALS(5, data.kl.at(Lvalues[i]).size());
                TS_ASSERT_LESS_THAN_EQUALS(data.kl.at(Lvalues[i]).size(),6);
                TS_ASSERT_EQUALS(data.kl.at(Lvalues[i]).count(Kvalues[j]), 1);
            }
        }
        fclose(in);
    }

    void test_read_mode_list_unordered() {
        // tests that duplicate L,K will be removed
        // and that they will be put in order
        Calculation::InputData data;
        std::string testfilename("tests/test_file.txt");
        std::string filecontents = "Frequencies: nonradial 5/3\n";
        int Lvalues[] = {3,1,2     /* duplicates to be ignored */, 3,1,2};
        int Kvalues[] = {4,2,3,5,1 /* duplicates to be ignored */, 5,5,1};
        for(int L : Lvalues){
            for(int K : Kvalues){
                filecontents += strmakef("%d,%d\n", L,K);
            }
        }
        make_test_input(testfilename,filecontents);
        FILE* in = fopen(testfilename.c_str(), "r");
        TS_ASSERT_EQUALS(0, io::read_frequencies(in, data));
        TS_ASSERT_EQUALS(data.mode_num, 17); //L=2,3 will have K=0 added
        TS_ASSERT_EQUALS(3, data.l.size());
        TS_ASSERT_EQUALS(3, data.kl.size());
        
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
            TS_ASSERT_EQUALS(Lout[i], Ldat[i]);
            TS_ASSERT_EQUALS(Kout[i], Kdat[i]);
        }
        fclose(in);
    }

    void test_read_mode_list_no_dipole_f() {
        // tests that duplicate L,K will be removed
        // and that they will be put in order
        Calculation::InputData data;
        std::string testfilename("tests/test_file.txt");
        std::string filecontents = "Frequencies: nonradial 5/3\n";
        int Lvalues[] = {1,2};
        int Kvalues[] = {0,1,2,3}; //the k=0 will be dropped when l=1
        for(int L : Lvalues){
            for(int K : Kvalues){
                filecontents += strmakef("%d,%d\n", L,K);
            }
        }
        make_test_input(testfilename,filecontents);
        FILE* in = fopen(testfilename.c_str(), "r");
        TS_ASSERT_EQUALS(0, io::read_frequencies(in, data));
        
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
            TS_ASSERT_EQUALS(Lout[i], Ldat[i]);
            TS_ASSERT_EQUALS(Kout[i], Kdat[i]);
        }
        fclose(in);
    }

    void test_read_bad_mode_list() {
        Calculation::InputData data;
        std::string testfilename("tests/test_file.txt");
        std::string filecontents = 
            "Frequencies: nonradial 5/3\n"
            "1,1\n"
            "1,2\n"
            "a,b\n"
            "1,3\n";
        make_test_input(testfilename,filecontents);
        FILE* in = fopen(testfilename.c_str(), "r");
        TS_ASSERT_EQUALS(0, io::read_frequencies(in, data));
        TS_ASSERT_EQUALS(data.mode_num, 3);
        TS_ASSERT_EQUALS(data.l.size(), 1);
        TS_ASSERT_EQUALS(data.kl[1].size(), 3);
        std::set<int> kforl = data.kl[1];
        std::vector<int> k(kforl.begin(), kforl.end());
        for(int i=0; i<3; i++){
            TS_ASSERT_EQUALS(k[i], i+1); // 1, 2, 3
        }
        fclose(in);
    }

    /* test reading polytrope input data */

    void test_fail_bad_params() {
        Calculation::InputData data;
        std::string testfilename("tests/test_file.txt");
        std::string filecontents = 
            "# test a comment\n\n"
            "Name: valid_test_name\n"
            "Model: newtonian polytrope 1.5 1\n"
            "Params: fakeparam1 1.0 fakeparm2 10.0\n";
        make_test_input(testfilename,filecontents);
        TS_ASSERT_EQUALS(1, io::read_input(testfilename.c_str(), data));
        TS_ASSERT_EQUALS(data.calcname, "valid_test_name");
        TS_ASSERT_EQUALS(data.regime, regime::PN0); 
        TS_ASSERT_EQUALS(data.model, model::polytrope);
        TS_ASSERT_EQUALS(data.input_params.size(), 2);
        TS_ASSERT_EQUALS(data.input_params[0], 1.5);
        TS_ASSERT_EQUALS(data.input_params[1], 1);
        TS_ASSERT_EQUALS(data.Ngrid, 1); 
    }

    void test_open_file_polytrope() {
        Calculation::InputData data;
        std::string testfilename("tests/test_file.txt");
        std::string filecontents = 
            "# test a comment\n\n"
            "Name: valid_test_name\n"
            "Model: newtonian polytrope 1.5 1\n"
            "Params: mass 1.0 radius 10.0\n"
            "Units: geo\n"
            "Frequencies: cowling 4/3\n\n";
        make_test_input(testfilename,filecontents);
        TS_ASSERT_EQUALS(0, io::read_input(testfilename.c_str(), data));
        TS_ASSERT_EQUALS(data.calcname, "valid_test_name");
        TS_ASSERT_EQUALS(data.regime, regime::PN0); 
        TS_ASSERT_EQUALS(data.model, model::polytrope);
        TS_ASSERT_EQUALS(data.input_params.size(), 2);
        TS_ASSERT_EQUALS(data.input_params[0], 1.5);
        TS_ASSERT_EQUALS(data.input_params[1], 1);
        TS_ASSERT_EQUALS(data.Ngrid, 1);
        TS_ASSERT_EQUALS(data.radius, 10.0);
        TS_ASSERT_EQUALS(data.mass, 1.0);
        TS_ASSERT_EQUALS(data.units, units::Units::geo);
        TS_ASSERT_EQUALS(data.modetype, modetype::cowling);
        TS_ASSERT_EQUALS(data.adiabatic_index, 4./3.);
        TS_ASSERT_EQUALS(data.mode_num, 0);
    }

    void test_open_file_polytrope_RZ() {
        Calculation::InputData data;
        std::string testfilename("tests/test_file.txt");
        std::string filecontents = 
            "# test a comment\n\n"
            "Name: valid_test_name\n"
            "Model: newtonian polytrope 1.5 1\n"
            "Params: zsurf 0.01 radius 10.0\n"
            "Units: geo\n"
            "Frequencies: cowling 4/3\n\n";
        make_test_input(testfilename,filecontents);
        TS_ASSERT_EQUALS(0, io::read_input(testfilename.c_str(), data));
        TS_ASSERT_EQUALS(data.radius, 10.0);
        TS_ASSERT_EQUALS(data.zsurf, 0.01);
    }

    void test_open_file_polytrope_MZ() {
        Calculation::InputData data;
        std::string testfilename("tests/test_file.txt");
        std::string filecontents = 
            "# test a comment\n\n"
            "Name: valid_test_name\n"
            "Model: newtonian polytrope 1.5 1\n"
            "Params: zsurf 0.01 mass 1.0\n"
            "Units: geo\n"
            "Frequencies: cowling 4/3\n\n";
        make_test_input(testfilename,filecontents);
        TS_ASSERT_EQUALS(0, io::read_input(testfilename.c_str(), data));
        TS_ASSERT_EQUALS(data.mass, 1.0);
        TS_ASSERT_EQUALS(data.zsurf, 0.01);
    }

    void test_open_file_polytrope_RLogg() {
        Calculation::InputData data;
        std::string testfilename("tests/test_file.txt");
        std::string filecontents = 
            "# test a comment\n\n"
            "Name: valid_test_name\n"
            "Model: newtonian polytrope 1.5 1\n"
            "Params: logg 10.0 radius 10.0\n"
            "Units: geo\n"
            "Frequencies: cowling 4/3\n\n";
        make_test_input(testfilename,filecontents);
        TS_ASSERT_EQUALS(0, io::read_input(testfilename.c_str(), data));
        TS_ASSERT_EQUALS(data.radius, 10.0);
        TS_ASSERT_EQUALS(data.logg, 10.0);
    }

    void test_open_file_polytrope_MLogg() {
        Calculation::InputData data;
        std::string testfilename("tests/test_file.txt");
        std::string filecontents = 
            "# test a comment\n\n"
            "Name: valid_test_name\n"
            "Model: newtonian polytrope 1.5 1\n"
            "Params: logg 10.0 mass 1.0\n"
            "Units: geo\n"
            "Frequencies: cowling 4/3\n\n";
        make_test_input(testfilename,filecontents);
        TS_ASSERT_EQUALS(0, io::read_input(testfilename.c_str(), data));
        TS_ASSERT_EQUALS(data.mass, 1.0);
        TS_ASSERT_EQUALS(data.logg, 10.0);
    }

    void test_open_file_polytrope_fail_RR() {
        Calculation::InputData data;
        std::string testfilename("tests/test_file.txt");
        std::string filecontents = 
            "# test a comment\n\n"
            "Name: valid_test_name\n"
            "Model: newtonian polytrope 1.5 1\n"
            "Params: radius 10.0 radius 10.0\n"
            "Units: geo\n"
            "Frequencies: cowling 4/3\n\n";
        make_test_input(testfilename,filecontents);
        TS_ASSERT_EQUALS(1, io::read_input(testfilename.c_str(), data));
    }

    void test_open_file_polytrope_fail_teff() {
        Calculation::InputData data;
        std::string testfilename("tests/test_file.txt");
        std::string filecontents = 
            "# test a comment\n\n"
            "Name: valid_test_name\n"
            "Model: newtonian polytrope 1.5 1\n"
            "Params: radius 10.0 teff 10000.0\n"
            "Units: geo\n"
            "Frequencies: cowling 4/3\n\n";
        make_test_input(testfilename,filecontents);
        TS_ASSERT_EQUALS(1, io::read_input(testfilename.c_str(), data));
    }

    /* test reading CHWD input  */
    /* test reading MESA input  */
    /* test reading SWD input  */

    /* test echo */
    
    void test_echo_input() {
        // create a test file, echo it, read in the echo, check match
        Calculation::InputData data;
        std::string testfilename("tests/test_file.txt");
        std::string filecontents = 
            "Name:\t../tests/test_file\n"
            "Model:\tnewtonian polytrope 1.50 1\n"
            "Params:\tmass 1\tradius 10\n"
            "Units:\tgeo\n\n"
            "Frequencies:\tcowling\t1.500\n";
        make_test_input(testfilename,filecontents);
        TS_ASSERT_EQUALS(0, io::read_input(testfilename.c_str(), data));
        system("mkdir -p tests/tests");
        std::string echoedcontents;
        std::string readfilename("tests/tests/test_file_in.txt");
        TS_ASSERT_EQUALS(0, io::echo_input(data));
        read_entire_file(readfilename, echoedcontents);
        TS_ASSERT_EQUALS(filecontents, echoedcontents);
        // clean up the file tree
        system("rm -r tests/tests");
        system("rm -r tests/test_file");
    }

    /* test setup output, which created output object */

    void test_setup_output() {
        Calculation::InputData indata;
        std::string testfilename("tests/test_file.txt");
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
        TS_ASSERT_EQUALS(0, io::read_input(testfilename.c_str(), indata));
        TS_ASSERT_EQUALS(indata.mode_num, 7);

        Calculation::OutputData outdata;
        outdata.star = nullptr;
        outdata.driver = nullptr;
        TS_ASSERT_EQUALS(0, io::setup_output(indata, outdata));
        TS_ASSERT_EQUALS(outdata.calcname, "TEST_SETUP_OUTPUT");
        TS_ASSERT_EQUALS(outdata.regime, regime::PN0); 
        TS_ASSERT_EQUALS(outdata.model, model::polytrope);
        TS_ASSERT_EQUALS(outdata.modetype, modetype::nonradial);
        TS_ASSERT_EQUALS(outdata.units, units::Units::geo);
        TS_ASSERT_EQUALS(outdata.input_params.size(), 2);
        TS_ASSERT_EQUALS(outdata.input_params[0], 1.5);
        TS_ASSERT_EQUALS(outdata.input_params[1], 1);
        TS_ASSERT_EQUALS(outdata.Ngrid, 1);
        TS_ASSERT_EQUALS(outdata.radius, 10.0);
        TS_ASSERT_EQUALS(outdata.mass, 1.0);
        TS_ASSERT_EQUALS(outdata.zsurf, units::getZsurfFromRM(outdata));
        TS_ASSERT_EQUALS(outdata.logg, units::getLoggFromRM(outdata));
        TS_ASSERT_EQUALS(outdata.mode_num, 7);
        TS_ASSERT_EQUALS(outdata.adiabatic_index, 5./3.);
        int Lout[] = {1,1,1, 2,2,2,2}; // in order, no duplicates
        int Kout[] = {1,2,3, 0,1,2,3}; // in order, no duplicates
        for(int i=0; i<7; i++){
            TS_ASSERT_EQUALS(Lout[i], outdata.l[i]);
            TS_ASSERT_EQUALS(Kout[i], outdata.k[i]);
        }
        TS_ASSERT_EQUALS(outdata.w.capacity(), 7);
        TS_ASSERT_EQUALS(outdata.f.capacity(), 7);
        TS_ASSERT_EQUALS(outdata.period.capacity(), 7);
        TS_ASSERT_EQUALS(outdata.mode_SSR.capacity(), 7);
        TS_ASSERT_EQUALS(outdata.i_err, 2);
        TS_ASSERT_EQUALS(outdata.error[0], true);
        TS_ASSERT_EQUALS(outdata.error[1], false);
        TS_ASSERT_EQUALS(outdata.error[2], false);
        TS_ASSERT_EQUALS(outdata.error[3], true);
    }

};
