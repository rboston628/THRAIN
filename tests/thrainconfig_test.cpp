#include "../src/ThrainMain.h"
#include "../src/ThrainMode.h"
#include "../src/ThrainUnits.h"
#include "../src/ThrainIO.h"
#include "../lib/filelib.h"
#include "./test_stars/DummyStar.h"
#include "doctest.h"

#include <iostream>
#include <fstream>

namespace {
  struct SetupTeardown {
    SetupTeardown() {
      original_config = ThrainConfig::getConfigOptions();
    }
    ~SetupTeardown() {
      ThrainConfig::reconfigure(original_config);
    }
    std::unordered_map<std::string, std::string> original_config;
  };
}

TEST_SUITE("ThrainConfig") {

  TEST_CASE_FIXTURE(SetupTeardown, "restore defaults") {
    // put it int oa state that is NOT the default
    std::unordered_map<std::string, std::string> newconfig{
      {"input_directory", "A"},
      {"output_directory", "B"},
      {"default_calc_name", "C"}
    };
    ThrainConfig::reconfigure(newconfig);
    // now restore the defaults
    ThrainConfig::restoreDefaults();
    auto config = ThrainConfig::getConfigOptions();
    // we cannot check against the defaults directly, but we can check they have changed
    for (const auto& keyval : newconfig) {
      CHECK_NE(config[keyval.first], keyval.second);
    }
  }

  TEST_CASE_FIXTURE(SetupTeardown, "get config options") {
    std::unordered_map<std::string, std::string> newconfig{
      {"input_directory", "A"},
      {"output_directory", "B"},
      {"default_calc_name", "C"}
    };
    ThrainConfig::reconfigure(newconfig);
    auto config = ThrainConfig::getConfigOptions();
    for (const auto& keyval : newconfig) {
      CHECK_EQ(config[keyval.first], keyval.second);
    }
  }

  TEST_CASE_FIXTURE(SetupTeardown, "reconfigure from file") {
    ThrainConfig::restoreDefaults();
    CHECK_NE(ThrainConfig::inputDir(), "testin");
    CHECK_NE(ThrainConfig::outputDir(), "testout");
    CHECK_NE(ThrainConfig::defaultCalcName(), "testname");
    ThrainConfig::reconfigure("tests/inputs/config_test.config");
    CHECK_EQ(ThrainConfig::inputDir(), "testin");
    CHECK_EQ(ThrainConfig::outputDir(), "testout");
    CHECK_EQ(ThrainConfig::defaultCalcName(), "testname");
  }

  TEST_CASE_FIXTURE(SetupTeardown, "reconfigure from map") {
    CHECK_NE(ThrainConfig::inputDir(), "map_input/");
    CHECK_NE(ThrainConfig::outputDir(), "map_output/");
    CHECK_NE(ThrainConfig::defaultCalcName(), "map_calc");
    std::unordered_map<std::string, std::string> config_map{
      {"input_directory", "map_input/"},
      {"output_directory", "map_output/"},
      {"default_calc_name", "map_calc"}
    };
    ThrainConfig::reconfigure(config_map);
    CHECK_EQ(ThrainConfig::inputDir(), "map_input/");
    CHECK_EQ(ThrainConfig::outputDir(), "map_output/");
    CHECK_EQ(ThrainConfig::defaultCalcName(), "map_calc");
  }

  TEST_CASE_FIXTURE(SetupTeardown, "resolveCalcName") {
    std::string name = "testcalc";
    std::string resolved = ThrainConfig::resolveCalcName(nullptr, name);
    CHECK_EQ(resolved, ThrainConfig::outputDir() + ThrainConfig::defaultCalcName() + "/" + name + "/");
    std::string custom = "customdir/";
    resolved = ThrainConfig::resolveCalcName(custom.c_str(), name);
    CHECK_EQ(resolved.back(), '/');
    CHECK_EQ(resolved, custom);
  }
}


/// This suite is meant to test the integration of ThrainConfig with the rest of the code, 
/// in particular ThrainIO.  Ensure ThrainIO and others are saving to locations specified by the config
TEST_SUITE("ThrainConfigIntegration") {
  // THRAIN IO READ_INPUT
  TEST_CASE_FIXTURE(SetupTeardown, "dir through thrainio read input") {
    // copy over the example input file to new directory
    ThrainConfig::reconfigure("./tests/tests.config");
    std::string source_calcdata = ThrainConfig::inputFileName("example_input.txt");
    std::string dest_dir = ThrainConfig::outputDir() + "config/";
    std::string dest_name = "copy_example_input.txt";
    {
      std::ifstream source(source_calcdata, std::ios::binary);
      filelib::makedir(dest_dir);
      std::ofstream destination(dest_dir + dest_name, std::ios::binary);
      REQUIRE(source);
      REQUIRE(destination);
      destination << source.rdbuf();
    }
    // setup the config to look there
    ThrainConfig::reconfigure({{"input_directory", dest_dir}});
    // now read it and verify
    Calculation::InputData calcdata;
    io::read_input(dest_name.c_str(), calcdata);
    CHECK_EQ(calcdata.calcname, "test_input_file");
    CHECK_EQ(calcdata.regime, regime::PN0);
    CHECK_EQ(calcdata.model, model::polytrope);
    CHECK_EQ(calcdata.units, units::Units::astro);
    CHECK_EQ(calcdata.mass, 0.6);
    CHECK_EQ(calcdata.radius, 8282);
    // cleanup
    filelib::remove(dest_dir);
  }
  // THRAIN IO ECHO_INPUT
  TEST_CASE_FIXTURE(SetupTeardown, "dir through thrainio echo") {
    // retrieve an input file to load a calcdata object
    Calculation::InputData calcdata;
    ThrainConfig::reconfigure("./tests/tests.config");
    io::read_input("example_input.txt", calcdata);
    calcdata.calcname = "testcalc";
    // predict where it will save
    std::string expected_path = ThrainConfig::echoedFileName(calcdata.calcname);
    // verify that echo creates the file in the expcted place
    if(filelib::exists(expected_path)) filelib::remove(expected_path);
    io::echo_input(calcdata);
    CHECK(filelib::exists(expected_path));
    // cleanup
    filelib::remove(ThrainConfig::calculationDir(calcdata.calcname));
  }
  // THRAIN IO WRITE_STELLAR_OUTPUT
  TEST_CASE_FIXTURE(SetupTeardown, "dir through thrainio write stellar") {
    // retrieve an input file to load a calcdata object
    Calculation::InputData calcdataIn;
    ThrainConfig::reconfigure("./tests/tests.config");
    io::read_input("example_input.txt", calcdataIn);
    calcdataIn.calcname = "testcalc";
    // from it create an output
    Calculation::OutputData calcdataOut;
    io::setup_output(calcdataIn, calcdataOut);
    calcdataOut.star = new DummyStar(100);
    calcdataOut.driver = nullptr;
    // predict where it will save
    std::string expected_path = ThrainConfig::summaryFileName(calcdataOut.calcname);
    if(filelib::exists(expected_path)) filelib::remove(expected_path);
    // save it and check it went to right place
    io::write_stellar_output(calcdataOut);
    CHECK(filelib::exists(expected_path));
    // cleanup
    filelib::remove(ThrainConfig::calculationDir(calcdataIn.calcname));
  }
  // THRAIN IO WRITE_MODE_OUTPUT
  TEST_CASE_FIXTURE(SetupTeardown, "dir through thrainio write mode") {
    // retrieve an input file to load a calcdata object
    Calculation::InputData calcdataIn;
    ThrainConfig::reconfigure("./tests/tests.config");
    io::read_input("example_input.txt", calcdataIn);
    calcdataIn.calcname = "testcalc";
    // from it create an output
    Calculation::OutputData calcdataOut;
    io::setup_output(calcdataIn, calcdataOut);
    calcdataOut.star = new DummyStar(100);
    calcdataOut.driver = nullptr;
    calcdataOut.mode_writ = 0;
    calcdataOut.mode_done = 0;
    ThrainConfig::reconfigure("./tests/tests.config");
    std::string expected_path = ThrainConfig::summaryFileName(calcdataOut.calcname);
    if(filelib::exists(expected_path)) filelib::remove(expected_path);
    // will return error code because the file doesn't exist
    CHECK_EQ(io::write_mode_output(calcdataOut), 1);
    // now make the file and then try
    filelib::touch(expected_path);
    REQUIRE(filelib::exists(expected_path));
    io::write_mode_output(calcdataOut);
    // have to check that the file has "Stellar Pulsation Results" printed in it
    {
      std::ifstream infile(expected_path);
      REQUIRE(infile);
      std::string line;
      bool found = false;
      while (std::getline(infile, line)) {
        if (line.find("Stellar Pulsation Results") != std::string::npos) {
          found = true;
          break;
        }
      }
      REQUIRE(found);
    }
    // cleanup
    filelib::remove(ThrainConfig::calculationDir(calcdataOut.calcname));
  }
  // THRAIN IO WRITE_TIDAL_OVERLAP
  TEST_CASE_FIXTURE(SetupTeardown, "dir through thrainio write tidal") {
    // retrieve an input file to load a calcdata object
    Calculation::InputData calcdataIn;
    ThrainConfig::reconfigure("./tests/tests.config");
    io::read_input("example_input.txt", calcdataIn);
    calcdataIn.calcname = "testcalc";
    // from it create an output
    Calculation::OutputData calcdataOut;
    io::setup_output(calcdataIn, calcdataOut);
    calcdataOut.star = new DummyStar(100);
    calcdataOut.driver = nullptr;
    calcdataOut.mode_writ = 0;
    calcdataOut.mode_done = 0;
    // predict where it will save
    std::string expected_path = ThrainConfig::calculationFileName(calcdataOut.calcname, "tidal_overlap.txt");
    if(filelib::exists(expected_path)) filelib::remove(expected_path);
    io::write_tidal_overlap(calcdataOut);
    CHECK(filelib::exists(expected_path));
    // cleanup
    filelib::remove(ThrainConfig::calculationDir(calcdataOut.calcname));
  }
  // POLYTROPE WRITES TO CORRECT LOCATION
  TEST_CASE_FIXTURE(SetupTeardown, "dir through polytrope write") {
    // retrieve an input file to load a calcdata object
    Calculation::InputData calcdataIn;
    ThrainConfig::reconfigure("./tests/tests.config");
    io::read_input("example_input.txt", calcdataIn);
    calcdataIn.calcname = "testcalc";
    // from it create an output calcdata object
    Calculation::OutputData calcdataOut;
    io::setup_output(calcdataIn, calcdataOut);
    // add a small polytrope star to the output
    calcdataOut.star = new Polytrope(1.5, 100);
    // make sure the driver is set to nullptr and mode_writ and mode_done are set to 0
    calcdataOut.driver = nullptr;
    calcdataOut.mode_writ = 0;
    calcdataOut.mode_done = 0;
    // setup the config so that the star prints at tests/output/test_config_poly
    ThrainConfig::reconfigure({"output_directory", "./tests/output/"});
    calcdataOut.calcname = "test_config_poly";
    // call to print the star
    io::write_stellar_output(calcdataOut);
    // verify the star outputs exist: 
    // - star/star.txt (and .png)
    // - star/BV.txt (and .png)
    // - star/wave_coefficients/ (should be at least four files here)
    // cleanup
  }



}

