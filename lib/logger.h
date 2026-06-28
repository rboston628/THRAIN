#ifndef THRAINLOGGER
#define THRAINLOGGER

#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <mutex>
#include <ctime>
#include <cstdarg>
#include "string.h"
#include "filelib.h"

class ThrainLogger {

public:
  enum class LogLevel : unsigned char { DEBUG = 0, INFO, WARNING, ERROR, MUTE = 0xFF };

private:
  bool logToFile;
  FILE *fp;
  std::mutex logMutex;
  LogLevel logLevel;

  // singleton
  static ThrainLogger& instance() {
    static ThrainLogger logger;
    return logger;
  }

  // private constructor, defaults to INFO level to stdout
  ThrainLogger() : logToFile(false), fp(stdout), logLevel(LogLevel::INFO), logMutex() {}

public:
  // delete copy and move operations
  ThrainLogger(ThrainLogger const&) = delete;
  ThrainLogger& operator=(ThrainLogger const&) = delete;
  ThrainLogger(ThrainLogger&&) = delete;
  ThrainLogger& operator=(ThrainLogger&&) = delete;

  // deconstructor
  ~ThrainLogger() {
    closeOpenLogFile();
  }

  static void setOutputFile(std::string const& filename, char const *access = "a") {
    // check if the file exists
    FILE *otherfp = fopen(filename.c_str(), access);
    if (!otherfp) {
      // if it failed, try making the directory
      std::string dir = filename.substr(0, filename.find_last_of("/\\"));
      filelib::makedir(dir);
      otherfp = fopen(filename.c_str(), access);
      if (!otherfp) {
        // if we still can't open the file, then we failed
        ThrainLogger::error("Failed to open log file %s; logger output file not changed\n", filename.c_str());
      }
      fprintf(stderr, "Failed to open log file %s; logger output file not changed\n", filename.c_str());
    } else {
      // close any open log file
      instance().closeOpenLogFile();
      // set the file
      std::lock_guard<std::mutex> lock(instance().logMutex);
      instance().fp = otherfp;
      instance().logToFile = true;
      // NOTE the file pointer is stored in singleton and deleted on destruction or reset
      otherfp = nullptr;
    }
  } // cppcheck-suppress resourceLeak

  static void unsetOutputFile() {
    instance().closeOpenLogFile();
    std::lock_guard<std::mutex> lock(instance().logMutex);
    instance().fp = stdout;
    instance().logToFile = false;
  }

  static void setLogLevel(LogLevel level) {
    instance().logLevel = level;
  }

  static void debug(char const *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    instance().innerlog(LogLevel::DEBUG, fmt, args);
    va_end(args);
  }

  static void info(char const *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    instance().innerlog(LogLevel::INFO, fmt, args);
    va_end(args);
  }

  static void warning(char const *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    instance().innerlog(LogLevel::WARNING, fmt, args);
    va_end(args);
  }

  static void error(char const *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    instance().innerlog(LogLevel::ERROR, fmt, args);
    va_end(args);
  }

  static void log(LogLevel level, char const *fmt, ...) {
    if(level != LogLevel::MUTE) {
      va_list args;
      va_start(args, fmt);
      instance().innerlog(level, fmt, args);
      va_end(args);
    }
  }

  static void logInline(LogLevel level, char const *fmt, ...) {
    if (level != LogLevel::MUTE && level >= instance().logLevel) {
      std::lock_guard<std::mutex> lock(instance().logMutex);
      va_list args;
      va_start(args, fmt);
      vfprintf(instance().fp, fmt, args);
      va_end(args);
      fflush(instance().fp);
    }
  }

private:

  void closeOpenLogFile() {
    if (logToFile && fp) {
      std::lock_guard<std::mutex> lock(logMutex);
      fclose(fp);
      logToFile = false;
    }
  }

  void innerlog(LogLevel level, char const *fmt, va_list args) {
    if (level >= logLevel ) {
      std::lock_guard<std::mutex> lock(logMutex);
      std::string log_fmt = strmakef("%s [%s]: %s", currentDateTime().c_str(), levelToString(level).c_str(), fmt);
      vfprintf(fp, log_fmt.c_str(), args);
      fflush(fp);
    }
  }

  std::string levelToString(LogLevel level) {
    switch (level) {
      case LogLevel::DEBUG:   return "DEBUG";
      case LogLevel::INFO:    return "INFO";
      case LogLevel::WARNING: return "WARNING";
      case LogLevel::ERROR:   return "ERROR";
      default:                return "UNKNOWN";
    }
  }

  std::string currentDateTime() {
    std::time_t now = std::time(nullptr);
    std::tm tm_info;
#if defined(_WIN32)
    localtime_s(&tm_info, &now);
#else
    localtime_r(&now, &tm_info);
#endif
    char buf[20] = {0};
    strftime(buf, 20, "%Y-%m-%d %H:%M:%S", &tm_info);
    return std::string(buf);
  }
}; // class ThrainLogger

#endif