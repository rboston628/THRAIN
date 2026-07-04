#include "doctest.h"
#include <cstdio>
#include <cstdlib>

#include "../lib/filelib.h"

namespace {
  std::string const test_output_dir ("tests/output/");
}

TEST_SUITE("FileLib [unit]") {

  TEST_CASE("filelib: exists on non-existing file") {
    std::string const testfile (test_output_dir + "hoopydoopydeedaw.txt");
    REQUIRE( !filelib::exists(testfile) );
  }

  TEST_CASE("filelib: exists on existing file") {
    std::string const testfile ("./src/ThrainMain.cpp");
    REQUIRE( filelib::exists(testfile) );
  }

  TEST_CASE("filelib: exists on non-existing directory") {
    std::string const testdir (test_output_dir + "hoopydoopydeedaw");
    REQUIRE( !filelib::exists(testdir) );
  }

  TEST_CASE("filelib: exists on existing directory") {
    std::string const testdir ("./src");
    REQUIRE( filelib::exists(testdir) );
  }

  TEST_CASE("filelib: makedir") {
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

  TEST_CASE("filelib: makedir recursive") {
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

  TEST_CASE("filelib: remove") {
    std::string const testdir (test_output_dir + "filelib_test_dir");
    filelib::makedir(testdir);
    REQUIRE( filelib::exists(testdir) );
    filelib::remove(testdir);
    REQUIRE( !filelib::exists(testdir) );
  }

  TEST_CASE("filelib: touch") {
    std::string const testfile (test_output_dir + "filelib_test_file.txt");
    filelib::remove(testfile);
    REQUIRE( !filelib::exists(testfile) );
    filelib::touch(testfile);
    REQUIRE( filelib::exists(testfile) );
    // cleanup
    filelib::remove(testfile);
  }

  TEST_CASE("filelib: makedir on existing directory") {
    std::string const testdir (test_output_dir + "filelib_test_dir");
    filelib::makedir(testdir);
    REQUIRE( filelib::exists(testdir) );
    filelib::makedir(testdir);
    REQUIRE( filelib::exists(testdir) );
    // cleanup
    filelib::remove(testdir);
  }

  TEST_CASE("filelib: remove on non-existing directory") {
    std::string const testdir (test_output_dir + "filelib_test_dir");
    filelib::remove(testdir);
    REQUIRE( !filelib::exists(testdir) );
    filelib::remove(testdir);
    REQUIRE( !filelib::exists(testdir) );
  }

  TEST_CASE("filelib: touch on existing file") {
    std::string const testfile (test_output_dir + "filelib_test_file.txt");
    filelib::touch(testfile);
    REQUIRE( filelib::exists(testfile) );
    filelib::touch(testfile);
    REQUIRE( filelib::exists(testfile) );
    // cleanup
    filelib::remove(testfile);
  }

  TEST_CASE("filelib: no infinite recursion") {
    std::string const testdir ("filelib_test_dir");
    // this will segfA=ault due to recursion depth if not passing
    filelib::makedir(testdir);
    REQUIRE( filelib::exists(testdir) );
    // cleanup
    filelib::remove(testdir);
  }
}