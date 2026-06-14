#ifndef THRAINCONFIG_H
#define THRAINCONFIG_H

#include <string>
#include <unordered_set>
#include <unordered_map>

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
  static void reconfigure(char const *const config_file);

  static void reconfigure(std::unordered_map<std::string, std::string> const& config_map);

  // getters

  /// Returns path to file used to set the current configuration
  static std::string configFile() { return instance().config_file; }
  /// Returns the input directory as set in the configuration
  static std::string inputDir() { return instance().config_options.at("input_directory"); }
  /// Returns the output directory as set in the configuration
  static std::string outputDir() { return instance().config_options.at("output_directory"); }
  static std::string defaultCalcName() { return instance().config_options.at("default_calc_name"); }
  /// Determine correct location for writing output; directory terminates with a "/"
  static std::string resolveCalcName(char const *const c, std::string const& name);

private:
  // loads the configuration from a file
  void setConfig(std::unordered_map<std::string, std::string> const& config_map);

  // reference to how the config was set
  std::string config_file;
  // configuration options
  static const std::unordered_set<std::string> valid_keys;
  std::unordered_map<std::string, std::string> config_options;
};

#endif