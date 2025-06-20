#pragma once

#include <string>
#include <map>
#include <vector>

class ConfigManager {
public:
    bool loadConfig(const std::string& configFile);
    void loadIsotopes(const std::string& isotopeFile);
    void loadIonList(const std::string& ionFile, std::vector<std::vector<std::string>>& ionlist);

    void printConfig() const;

    const std::map<std::string, std::string>& getOptions() const;
    std::string  getIsotope(const std::string& key) const;
    std::string  get(const std::string& key) const;
    bool     getBool(const std::string& key) const;
    double getDouble(const std::string& key) const;
    int       getInt(const std::string& key) const;

private:
    std::map<std::string, std::string> options;
    std::map<std::string, std::string> isotopes;
};
