#ifndef THRAINCONFIG_H
#define THRAINCONFIG_H

#include <string>
#include <stdexcept>
#include <cstdio>
#include <unordered_set>
#include <unordered_map>
#include "logger.h"

class ThrainConfig {
private:
  // singleton
  static ThrainConfig& instance() {
    static ThrainConfig config;
    return config;
  }

  // private constructor to make singleton
  ThrainConfig() : config_file("defaults"), 
    config_options({{"input_directory", "./"}, {"output_directory", "./output/"}, {"default_calc_name", "models"}}) {}

public:
  // complete fule of five for singleton pattern
  ThrainConfig(ThrainConfig const&) = delete;
  ThrainConfig& operator=(ThrainConfig const&) = delete;
  ThrainConfig(ThrainConfig&&) = delete;
  ThrainConfig& operator=(ThrainConfig&&) = delete;

  // setters

  /// @brief Reconfigures the application from a configuration file
  /// @param config_file path to the configuration file
  /// The configuration file should be a text file with lines of the form "key value"
  static void reconfigure(char const *const config_file) {
    std::unorderded_map<std::string, std::string> config_map;
    FILE *config_fp = fopen(config_file, "r");
    if (!config_fp) {
      ThrainLogger::error("Could not open config file %s\n", config_file);
      return;
    } else {
      char key[256], value[256];
      while (fscanf(config_fp, "%255s %255s", key, value) == 2) {
        config_map[key] = value;
      }
      fclose(config_fp);
    }
    instance().setConfig(config_map);
    instance().config_file = config_file;
  }

  static void reconfigure(std::unordered_map<std::string, std::string> const& config_map) {
    instance().setConfig(config_map);
    instance().config_file = "config from map";
  }

  // getters

  /// Returns path to file used to set the current configuration
  static std::string configFile() { return instance().config_file; }
  /// Returns the input directory as set in the configuration
  static std::string inputDir() { return instance().config_options.at("input_directory"); }
  /// Returns the output directory as set in the configuration
  static std::string outputDir() { return instance().config_options.at("output_directory"); }
  static std::string defaultCalcName() { return instance().config_options.at("default_calc_name"); }
  /// Determine correct location for writing output; directory terminates with a "/"
  static std::string resolveCalcName(char const *const c, std::string const& name) {
    std::string ret;
    if (c) return ret = c;
    else ret = outputDir() + defaultCalcName() + "/" + name + "/";
    if (ret.back() != '/') ret += '/';
    return ret;
  }

private:
  // loads the configuration from a file
  void setConfig(std::unordered_map<std::string, std::string> const& config_map) {
    for (const auto& [key, value] : config_map) {
        if (valid_keys.contains(key)) {
          instance().config_options[key] = value;
        } else {
          ThrainLogger::warning("Unknown config key %s in config map\n", key.c_str());
      }
    }
  }

  // reference to how the config was set
  std::string config_file;
  // configuration options
  static const std::unordered_set<std::string> valid_keys = {"input_directory", "output_directory", "default_calc_name"};
  std::unordered_map<std::string, std::string> config_options;
};
#endif