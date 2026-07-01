#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"
#include "../src/ThrainConfig.h"
#include "../lib/logger.h"

int main(int argc, char **argv) {
  // use the test configuration file. this ensures read/write from files within test dir
  ThrainConfig::reconfigure("./tests/tests.config");

  // setup this directory, to be used for output
  filelib::makedir(ThrainConfig::outputDir());
  // setup the artifact directory, used by some tests
  filelib::makedir(ThrainConfig::outputDir() + "artifacts/");

  // set the log level.  defaults to mute.
  int out = 1;
  ThrainLogger::setLogLevel(ThrainLogger::LogLevel::MUTE);
  for (int in = 1; in < argc; ++in) {
    if (strstr(argv[in], "--verbose")) {
        ThrainLogger::setLogLevel(ThrainLogger::LogLevel::DEBUG);
    } else {
        argv[out++] = argv[in];
    }
  }

  // forward the remaining arguments to doctest
  argc = out;
  argv[argc] = nullptr;
  return doctest::Context(argc, argv).run();
}