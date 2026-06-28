// Tests of the logger

#include "../src/constants.h"
#include "../lib/string.h"
#include "../src/ThrainConfig.h"
#include "../lib/filelib.h"
#include <cxxtest/TestSuite.h>

// define this for the tests below to work
ThrainLogger::LogLevel operator+(ThrainLogger::LogLevel l, unsigned char i) {
  unsigned char val = static_cast<unsigned char>(l);
  return static_cast<ThrainLogger::LogLevel>(val + i);
}

namespace {
  std::string const testfilename (ThrainConfig::outputDir() + "log_test.txt");
}

class LoggerTest : public CxxTest::TestSuite {
public:

static LoggerTest *createSuite (){
  filelib::makedir(ThrainConfig::outputDir());
  printf("\n## LOGGER TESTS ##\n");
  return new LoggerTest();
}

static void destroySuite(LoggerTest *suite) { 
  printf("\n################");
  filelib::remove(testfilename.c_str());
  delete suite; 
}

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
  // redirect outout to file for checking
  ThrainLogger::setOutputFile(testfilename, "w");
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
  // redirect outout to file for checking
  ThrainLogger::setOutputFile(testfilename, "w");
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
  char line[128] = {'\0'};
  int Y, M, D, h, m, s, res;
  for (int i=static_cast<int>(level); i<=3; i++) {
    fscanf(in, "%*4d-%*2d-%*2d %*2d:%*2d:%*2d [%*[A-Z]]: %d\n", &res);
    TS_ASSERT_EQUALS(res, i + 1);
  }
  fgets(line, 120, in);
  TS_ASSERT_EQUALS(line[0], '\0');
}

void check_file(char const *expected) {
  FILE *in = fopen(testfilename.c_str(), "r");
  char line[128] = {'\0'};
  fgets(line, 128, in);
  fclose(in);
  TS_ASSERT_EQUALS(line, expected);
}

void test_set_output() {
  ThrainLogger::setOutputFile(testfilename, "w");
  char expected[] = "test";
  ThrainLogger::logInline(ThrainLogger::LogLevel::INFO, expected);
  ThrainLogger::unsetOutputFile();
  ThrainLogger::logInline(ThrainLogger::LogLevel::INFO, "logged to stdout");
  
  // check output
  check_file(expected);
}

// debug

void test_debug() {
  ThrainLogger::setOutputFile(testfilename, "w");
  
  log_all_levels(ThrainLogger::LogLevel::DEBUG);
  // set the log level to info
  ThrainLogger::setLogLevel(ThrainLogger::LogLevel::INFO);
  // this will not be printed
  ThrainLogger::debug("not printed");
  ThrainLogger::unsetOutputFile();

  check_logged_level(ThrainLogger::LogLevel::DEBUG);
}

void test_log_debug() {
  log_func_all_levels(ThrainLogger::LogLevel::DEBUG);
  check_logged_level(ThrainLogger::LogLevel::DEBUG);
}

void test_loginline_debug() {
  loginline_func_all_levels(ThrainLogger::LogLevel::DEBUG);
  check_file("1234");
}

// info

void test_info() {
  ThrainLogger::setOutputFile(testfilename, "w");
  log_all_levels(ThrainLogger::LogLevel::INFO);
  // set the log level to warning
  ThrainLogger::setLogLevel(ThrainLogger::LogLevel::WARNING);
  // this will not be printed
  ThrainLogger::info("not printed");
  ThrainLogger::unsetOutputFile();
  check_logged_level(ThrainLogger::LogLevel::INFO);
}

void test_log_info() {
  log_func_all_levels(ThrainLogger::LogLevel::INFO);
  check_logged_level(ThrainLogger::LogLevel::INFO);
}

void test_loginline_info() {
  loginline_func_all_levels(ThrainLogger::LogLevel::INFO);
  check_file("234");
}

// warning

void test_warning() {
  ThrainLogger::setOutputFile(testfilename, "w");
  
  log_all_levels(ThrainLogger::LogLevel::WARNING);
  // set the log level to error
  ThrainLogger::setLogLevel(ThrainLogger::LogLevel::ERROR);
  // this will not be printed
  ThrainLogger::warning("not printed");
  ThrainLogger::unsetOutputFile();

  check_logged_level(ThrainLogger::LogLevel::WARNING);
}

void test_log_warning() {
  log_func_all_levels(ThrainLogger::LogLevel::WARNING);
  check_logged_level(ThrainLogger::LogLevel::WARNING);
}

void test_loginline_warning() {
  loginline_func_all_levels(ThrainLogger::LogLevel::WARNING);
  check_file("34");
}

// error

void test_error() {
  ThrainLogger::setOutputFile(testfilename, "w");
  
  log_all_levels(ThrainLogger::LogLevel::ERROR);
  // set the log level to error
  ThrainLogger::setLogLevel(ThrainLogger::LogLevel::MUTE);
  // this will not be printed
  ThrainLogger::error("not printed");
  ThrainLogger::unsetOutputFile();

  check_logged_level(ThrainLogger::LogLevel::ERROR);
}

void test_log_error() {
  log_func_all_levels(ThrainLogger::LogLevel::ERROR);
  check_logged_level(ThrainLogger::LogLevel::ERROR);
}

void test_loginline_error() {
  loginline_func_all_levels(ThrainLogger::LogLevel::ERROR);
  check_file("4");
}

// mute

void test_mute() {
  ThrainLogger::setOutputFile(testfilename, "w");
  log_all_levels(ThrainLogger::LogLevel::MUTE);
  ThrainLogger::unsetOutputFile();
  check_file("");
}

void test_log_mute() {
  ThrainLogger::setOutputFile(testfilename, "w");
  log_func_all_levels(ThrainLogger::LogLevel::MUTE);
  ThrainLogger::unsetOutputFile();
  check_file("");
}

void test_loginline_mute() {
  ThrainLogger::setOutputFile(testfilename, "w");
  loginline_func_all_levels(ThrainLogger::LogLevel::MUTE);
  ThrainLogger::unsetOutputFile();
  check_file("");
}

};
