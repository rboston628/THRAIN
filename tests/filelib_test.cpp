#include "doctest.h"
#include <cstdio>
#include <cstdlib>

#include "../lib/filelib.h"

namespace {
  std::string const test_output_dir ("tests/output/");
}

TEST_SUITE("FileLib") {

  TEST_CASE("exists on non-existing file") {
    std::string const testfile (test_output_dir + "hoopydoopydeedaw.txt");
    REQUIRE( !filelib::exists(testfile) );
  }

  TEST_CASE("exists on existing file") {
    std::string const testfile ("./src/ThrainMain.cpp");
    REQUIRE( filelib::exists(testfile) );
  }

  TEST_CASE("exists on non-existing directory") {
    std::string const testdir (test_output_dir + "hoopydoopydeedaw");
    REQUIRE( !filelib::exists(testdir) );
  }

  TEST_CASE("exists on existing directory") {
    std::string const testdir ("./src");
    REQUIRE( filelib::exists(testdir) );
  }

  TEST_CASE("makedir") {
    if (!filelib::exists(test_output_dir)) {
      filelib::makedir(test_output_dir);
      REQUIRE( filelib::exists(test_output_dir) );
    }
    std::string const testdir (test_output_dir + "filelib_test_dir");
    filelib::remove(testdir);
    REQUIRE( !filelib::exists(testdir) );
    filelib::makedir(testdir);
    REQUIRE( filelib::exists(testdir) );
    // cleanup
    filelib::remove(testdir);
  }

  TEST_CASE("makedir recursive") {
    std::string const parentdir (test_output_dir + "level1/");
    std::string const childdir (parentdir + "level2/");
    filelib::remove(parentdir);
    REQUIRE( !filelib::exists(parentdir) );
    REQUIRE( !filelib::exists(childdir) );
    filelib::makedir(childdir);
    REQUIRE( filelib::exists(parentdir) );
    REQUIRE( filelib::exists(childdir) );
    // cleanup
    filelib::remove(parentdir);
  }

  TEST_CASE("remove") {
    std::string const testdir (test_output_dir + "filelib_test_dir");
    filelib::makedir(testdir);
    REQUIRE( filelib::exists(testdir) );
    filelib::remove(testdir);
    REQUIRE( !filelib::exists(testdir) );
  }

  TEST_CASE("touch") {
    std::string const testfile (test_output_dir + "filelib_test_file.txt");
    filelib::remove(testfile);
    REQUIRE( !filelib::exists(testfile) );
    filelib::touch(testfile);
    REQUIRE( filelib::exists(testfile) );
    // cleanup
    filelib::remove(testfile);
  }

  TEST_CASE("makedir on existing directory") {
    std::string const testdir (test_output_dir + "filelib_test_dir");
    filelib::makedir(testdir);
    REQUIRE( filelib::exists(testdir) );
    filelib::makedir(testdir);
    REQUIRE( filelib::exists(testdir) );
    // cleanup
    filelib::remove(testdir);
  }

  TEST_CASE("remove on non-existing directory") {
    std::string const testdir (test_output_dir + "filelib_test_dir");
    filelib::remove(testdir);
    REQUIRE( !filelib::exists(testdir) );
    filelib::remove(testdir);
    REQUIRE( !filelib::exists(testdir) );
  }

  TEST_CASE("touch on existing file") {
    std::string const testfile (test_output_dir + "filelib_test_file.txt");
    filelib::touch(testfile);
    REQUIRE( filelib::exists(testfile) );
    filelib::touch(testfile);
    REQUIRE( filelib::exists(testfile) );
    // cleanup
    filelib::remove(testfile);
  }
}