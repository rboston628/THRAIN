#ifndef THRAINCONFIG_CPP
#define THRAINCONFIG_CPP

#include "ThrainConfig.h"

#include <stdexcept>
#include <cstdio>
#include "../lib/logger.h"

/// @brief Reconfigures the application from a configuration file
/// @param config_file path to the configuration file
/// The configuration file should be a text file with lines of the form "key value"
void ThrainConfig::reconfigure(char const *const config_file) {
  std::unordered_map<std::string, std::string> config_map;
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

void ThrainConfig::reconfigure(std::unordered_map<std::string, std::string> const& config_map) {
  instance().setConfig(config_map);
  instance().config_file = "config from map";
}

/// Determine correct location for writing output; directory terminates with a "/"
std::string ThrainConfig::resolveCalcName(char const *const c, std::string const& name) {
  std::string ret;
  if (c) return ret = c;
  else ret = outputDir() + defaultCalcName() + "/" + name + "/";
  if (ret.back() != '/') ret += '/';
  return ret;
}

// loads the configuration from a file
void ThrainConfig::setConfig(std::unordered_map<std::string, std::string> const& config_map) {
  for (const auto& keyval : config_map) {
    if (valid_keys.count(keyval.first)) {
      instance().config_options[keyval.first] = keyval.second;
    } else {
      ThrainLogger::warning("Unknown config key %s in config map\n", keyval.first.c_str());
    }
  }
}

std::unordered_set<std::string> const ThrainConfig::valid_keys{"input_directory", "output_directory", "default_calc_name"};

#endif