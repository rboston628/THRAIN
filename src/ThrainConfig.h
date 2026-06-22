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
    config_options(default_config) {}

  // complete fule of five for singleton pattern
  ThrainConfig(ThrainConfig const&) = delete;
  ThrainConfig& operator=(ThrainConfig const&) = delete;
  ThrainConfig(ThrainConfig&&) = delete;
  ThrainConfig& operator=(ThrainConfig&&) = delete;

public:
  // setters

  /// @brief Reconfigures the application from a configuration file
  /// @param config_file path to the configuration file
  /// The configuration file should be a text file with lines of the form "key value"
  static void reconfigure(char const *const config_file);

  static void reconfigure(std::unordered_map<std::string, std::string> const& config_map);

  // return the config to the default parameters
  static void restoreDefaults() {
    instance().config_options = default_config;
    instance().config_file = "defaults";
  }

  // getters

  /// Returns path to file used to set the current configuration
  static std::string configFile() { return instance().config_file; }
  /// Returns the input directory as set in the configuration
  static std::string inputDir() { return trailingSlash(instance().config_options.at("input_directory")); }
  /// Returns the output directory as set in the configuration
  static std::string outputDir() { return trailingSlash(instance().config_options.at("output_directory")); }
  static std::string defaultCalcName() { return instance().config_options.at("default_calc_name"); }
  
  static std::string inputFileName(std::string const& name) {
    return inputDir() + name;
  }
  static std::string calculationDir(std::string const& calcname) {
    return outputDir() + trailingSlash(calcname);
  }
  static std::string echoedFileName(std::string const& calcname) {
    return calculationDir(calcname) + calcname + "_in.txt";
  }
  static std::string summaryFileName(std::string const& calcname) {
    return calculationDir(calcname) + calcname + ".txt";
  }
  static std::string calculationSubdir(std::string const& calcname, std::string const& subdir) {
    return calculationDir(calcname) + trailingSlash(subdir);
  }
  static std::string calculationFileName(std::string const& calcname, std::string const& filename) {
    return calculationDir(calcname) + filename;
  }
  static std::string calculationFileName(std::string const& calcname, std::string const& subdir, std::string const& filename) {
    return calculationDir(calcname) + trailingSlash(subdir) + filename;
  } 
  
  /// Determine correct location for writing output; directory terminates with a "/"
  static std::string resolveCalcName(std::string const& c, std::string const& name);

  static std::unordered_map<std::string, std::string> getConfigOptions() { return instance().config_options; }

private:
  static std::string trailingSlash(std::string const& s) {
    return s.back() == '/' ? s : s + '/';
  }
  // loads the configuration from a file
  void setConfig(std::unordered_map<std::string, std::string> const& config_map);

  // reference to how the config was set
  std::string config_file;
  // configuration options
  static const std::unordered_set<std::string> valid_keys;
  static const std::unordered_map<std::string, std::string> default_config;
  std::unordered_map<std::string, std::string> config_options;
};

#endif