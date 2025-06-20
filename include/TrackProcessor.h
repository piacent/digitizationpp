#pragma once

#include <vector>
#include <string>
#include <map>

class TrackProcessor {
public:
    TrackProcessor(const std::map<std::string, std::string>& config);

    void computeWithSaturation(const std::vector<double>& x_hits,
                               const std::vector<double>& y_hits,
                               const std::vector<double>& z_hits,
                               const std::vector<double>& energy_hits,
                               float energy,
                               bool isNR,
                               std::vector<std::vector<double>>& image);

    void computeWithoutSaturation(const std::vector<double>& x_hits,
                                  const std::vector<double>& y_hits,
                                  const std::vector<double>& z_hits,
                                  const std::vector<double>& energy_hits,
                                  std::vector<std::vector<double>>& image);

private:
    const std::map<std::string, std::string>& config;
    void smearHits(const std::vector<double>& x,
                   const std::vector<double>& y,
                   const std::vector<double>& z,
                   const std::vector<double>& nel,
                   std::vector<float>& Sx,
                   std::vector<float>& Sy,
                   std::vector<float>& Sz);

    std::vector<double> computeIonizationElectrons(const std::vector<double>& energyDep,
                                                   const std::vector<double>& z_hit);
};
