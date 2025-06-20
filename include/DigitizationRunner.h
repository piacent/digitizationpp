#pragma once

#include "ConfigManager.h"
#include "Globals.h"
#include "cygnolib.h"
#include <string>
#include <vector>
#include <map>
#include <memory>
#include "TH2F.h"

// Forward declarations
class TFile;

class DigitizationRunner {
public:
    DigitizationRunner(const std::string& configFile, const std::string& inputDir, const std::string& outputDir);
    void run();

private:
    void initializeGlobals();
    void prepareCameraSettings();
    void checkDimensionConsistency();
    void processRootFiles();
    bool isValidInputFile(const std::string& filename) const;

    void setSeed(int seed = 10);
    void loadIonList4SRIM(ConfigManager& config, std::vector<std::vector<std::string>>& SRIM_events, const std::string& filename, const std::string& infolder);
    void initSourceDir();
    void SaveValues(std::shared_ptr<TFile>& outfile);
    void fillVigmap(TH2F& VignMap);
    bool is_NR(std::vector<int> pdgID_hits, int pdg);
    void AddBckg(std::vector<std::vector<int>>& background);


    ConfigManager config;
    std::string configFile;
    std::string inputDir;
    std::string outputDir;
    int runCount=0;

    int NMAX_EVENTS = 10;
    std::string SOURCE_DIR = "";

};
