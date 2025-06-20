#pragma once

#include "ConfigManager.h"
#include "Globals.h"
#include <string>
#include <vector>
#include <map>
#include <memory>

// Forward declarations
class TFile;

class DigitizationRunner {
public:
    DigitizationRunner(ConfigManager& config, const std::string& inputDir, const std::string& outputDir);
    void run();

private:
    void initializeGlobals();
    void prepareCameraSettings();
    void checkDimensionConsistency();
    void processRootFiles();
    bool isValidInputFile(const std::string& filename) const;
    void setSeed(int seed = 10);

    ConfigManager& config;
    std::string inputDir;
    std::string outputDir;
    int runCount;
};
