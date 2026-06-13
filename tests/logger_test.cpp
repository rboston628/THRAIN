// Tests of the logger

#include "../src/constants.h"
#include "../lib/string.h"
#include "doctest.h"
#include <filesystem>

// define this for the tests below to work
ThrainLogger::LogLevel operator+(ThrainLogger::LogLevel l, unsigned char i) {
  unsigned char val = static_cast<unsigned char>(l);
  return static_cast<ThrainLogger::LogLevel>(val + i);
}

namespace {
  std::string const testfilename ("tests/tests/log_test.txt");

  void log_all_levels(ThrainLogger::LogLevel const level) {
    // set the log level
    ThrainLogger::setLogLevel(level);
    // log at all levels
    ThrainLogger::debug("%d\n", 1);
    ThrainLogger::info("%d\n", 2);
    ThrainLogger::warning("%d\n", 3);
    ThrainLogger::error("%d\n", 4);
  }

  void log_func_all_levels(ThrainLogger::LogLevel const level) {
    // set the log level
    ThrainLogger::setLogLevel(level);
    // log at all levels
    for (auto l = ThrainLogger::LogLevel::DEBUG; l <= ThrainLogger::LogLevel::ERROR; l=l+1) {
      ThrainLogger::log(l, "%d\n", static_cast<int>(l) + 1);
    }
    ThrainLogger::setLogLevel(level + 1);
    // this will not be printed
    ThrainLogger::log(level, "not printed");
    ThrainLogger::unsetOutputFile();
  }

  void loginline_func_all_levels(ThrainLogger::LogLevel const level) {
    // set the log level
    ThrainLogger::setLogLevel(level);
    // log at all levels
    for (auto l = ThrainLogger::LogLevel::DEBUG; l <= ThrainLogger::LogLevel::ERROR; l=l+1) {
      ThrainLogger::logInline(l, "%d", static_cast<int>(l) + 1);
    }
    ThrainLogger::setLogLevel(level + 1);
    // this will not be printed
    ThrainLogger::logInline(level, "not printed");
    ThrainLogger::unsetOutputFile();
  }

  void check_logged_level(ThrainLogger::LogLevel const level) {
    FILE *in = fopen(testfilename.c_str(), "r");
    REQUIRE(in != nullptr);
    char line[128] = {'\0'};
    int res;
    for (int i=static_cast<int>(level); i<=3; i++) {
      int n = fscanf(in, "%*4d-%*2d-%*2d %*2d:%*2d:%*2d [%*[A-Z]]: %d\n", &res);
      REQUIRE(n == 1);
      CHECK_EQ(res, i + 1);
    }
    char *s = fgets(line, 120, in);
    if (s == nullptr) {
      REQUIRE_EQ(line[0], '\0');
    }
    CHECK_EQ(line[0], '\0');
    fclose(in);
  }

  void check_file(char const *expected) {
    FILE *in = fopen(testfilename.c_str(), "r");
    REQUIRE(in != nullptr);
    char line[128] = {'\0'};
    char *s = fgets(line, 128, in);
    fclose(in);

    if (expected[0] == '\0') {
      CHECK(s == nullptr || line[0] == '\0');
    } else {
      REQUIRE(s != nullptr);
      CHECK_EQ(std::string(line), std::string(expected));
    }
  }
} // namespace

struct SetUpTearDown {
  SetUpTearDown() { 
    std::filesystem::path(testfilename).parent_path().create_directories();
    ThrainLogger::setOutputFile(testfilename, "w");
  }
  ~SetUpTearDown() { remove(testfilename.c_str()); }
};

TEST_SUITE("Logger") {

TEST_CASE_FIXTURE(SetUpTearDown, "set_output") {
  char expected[] = "test";
  ThrainLogger::logInline(ThrainLogger::LogLevel::INFO, expected);
  ThrainLogger::unsetOutputFile();
  ThrainLogger::logInline(ThrainLogger::LogLevel::INFO, "logged to stdout");
  
  // check output
  check_file(expected);
}

// debug

TEST_CASE_FIXTURE(SetUpTearDown, "debug") {
  log_all_levels(ThrainLogger::LogLevel::DEBUG);
  // set the log level to info
  ThrainLogger::setLogLevel(ThrainLogger::LogLevel::INFO);
  // this will not be printed
  ThrainLogger::debug("not printed");
  ThrainLogger::unsetOutputFile();

  check_logged_level(ThrainLogger::LogLevel::DEBUG);
}

TEST_CASE_FIXTURE(SetUpTearDown, "log_debug") {
  log_func_all_levels(ThrainLogger::LogLevel::DEBUG);
  check_logged_level(ThrainLogger::LogLevel::DEBUG);
}

TEST_CASE_FIXTURE(SetUpTearDown, "loginline_debug") {
  loginline_func_all_levels(ThrainLogger::LogLevel::DEBUG);
  check_file("1234");
}

// info

TEST_CASE_FIXTURE(SetUpTearDown, "info") {
  log_all_levels(ThrainLogger::LogLevel::INFO);
  // set the log level to warning
  ThrainLogger::setLogLevel(ThrainLogger::LogLevel::WARNING);
  // this will not be printed
  ThrainLogger::info("not printed");
  ThrainLogger::unsetOutputFile();
  check_logged_level(ThrainLogger::LogLevel::INFO);
}

TEST_CASE_FIXTURE(SetUpTearDown, "log_info") {
  log_func_all_levels(ThrainLogger::LogLevel::INFO);
  check_logged_level(ThrainLogger::LogLevel::INFO);
}

TEST_CASE_FIXTURE(SetUpTearDown, "loginline_info") {
  loginline_func_all_levels(ThrainLogger::LogLevel::INFO);
  check_file("234");
}

// warning

TEST_CASE_FIXTURE(SetUpTearDown, "warning") {
  log_all_levels(ThrainLogger::LogLevel::WARNING);
  // set the log level to error
  ThrainLogger::setLogLevel(ThrainLogger::LogLevel::ERROR);
  // this will not be printed
  ThrainLogger::warning("not printed");
  ThrainLogger::unsetOutputFile();

  check_logged_level(ThrainLogger::LogLevel::WARNING);
}

TEST_CASE_FIXTURE(SetUpTearDown, "log_warning") {
  log_func_all_levels(ThrainLogger::LogLevel::WARNING);
  check_logged_level(ThrainLogger::LogLevel::WARNING);
}

TEST_CASE_FIXTURE(SetUpTearDown, "loginline_warning") {
  loginline_func_all_levels(ThrainLogger::LogLevel::WARNING);
  check_file("34");
}

// error

TEST_CASE_FIXTURE(SetUpTearDown, "error") {
  log_all_levels(ThrainLogger::LogLevel::ERROR);
  // set the log level to error
  ThrainLogger::setLogLevel(ThrainLogger::LogLevel::MUTE);
  // this will not be printed
  ThrainLogger::error("not printed");
  ThrainLogger::unsetOutputFile();

  check_logged_level(ThrainLogger::LogLevel::ERROR);
}

TEST_CASE_FIXTURE(SetUpTearDown, "log_error") {
  log_func_all_levels(ThrainLogger::LogLevel::ERROR);
  check_logged_level(ThrainLogger::LogLevel::ERROR);
}

TEST_CASE_FIXTURE(SetUpTearDown, "loginline_error") {
  loginline_func_all_levels(ThrainLogger::LogLevel::ERROR);
  check_file("4");
}

// mute

TEST_CASE_FIXTURE(SetUpTearDown, "mute") {
  log_all_levels(ThrainLogger::LogLevel::MUTE);
  ThrainLogger::unsetOutputFile();
  check_file("");
}

TEST_CASE_FIXTURE(SetUpTearDown, "log_mute") {
  log_func_all_levels(ThrainLogger::LogLevel::MUTE);
  ThrainLogger::unsetOutputFile();
  check_file("");
}

TEST_CASE_FIXTURE(SetUpTearDown, "loginline_mute") {
  loginline_func_all_levels(ThrainLogger::LogLevel::MUTE);
  ThrainLogger::unsetOutputFile();
  check_file("");
}

} // end TEST_SUITE
