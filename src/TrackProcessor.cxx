#include "TrackProcessor.h"
#include "Globals.h"
#include <cmath>
#include <numeric>
#include <random>
#include <stdexcept>
#include <iostream>

TrackProcessor::TrackProcessor(const std::map<std::string, std::string>& configMgr)
    : config(configMgr) {}

void TrackProcessor::computeWithSaturation(const std::vector<double>& x_hits,
                                           const std::vector<double>& y_hits,
                                           const std::vector<double>& z_hits,
                                           const std::vector<double>& energy_hits,
                                           float energy,
                                           bool isNR,
                                           std::vector<std::vector<double>>& image)
{
    // Simplified saturation application — this would be extended to match original logic
    std::vector<float> Sx, Sy, Sz;

    auto nel = computeIonizationElectrons(energy_hits, z_hits);
    smearHits(x_hits, y_hits, z_hits, nel, Sx, Sy, Sz);

    for (size_t i = 0; i < Sx.size(); ++i) {
        int xi = static_cast<int>(Sx[i] * x_pix / std::stod(config.at("x_dim")));
        int yi = static_cast<int>(Sy[i] * y_pix / std::stod(config.at("y_dim")));
        if (xi >= 0 && xi < x_pix && yi >= 0 && yi < y_pix)
            image[xi][yi] += 1.0;
    }
}

void TrackProcessor::computeWithoutSaturation(const std::vector<double>& x_hits,
                                              const std::vector<double>& y_hits,
                                              const std::vector<double>& z_hits,
                                              const std::vector<double>& energy_hits,
                                              std::vector<std::vector<double>>& image)
{
    computeWithSaturation(x_hits, y_hits, z_hits, energy_hits, 0.0f, false, image);
}

void TrackProcessor::smearHits(const std::vector<double>& x,
                               const std::vector<double>& y,
                               const std::vector<double>& z,
                               const std::vector<double>& nel,
                               std::vector<float>& Sx,
                               std::vector<float>& Sy,
                               std::vector<float>& Sz)
{
    std::default_random_engine gen;
    std::normal_distribution<float> dist(0.0, 1.0);

    double sigma_xy = std::stod(config.at("diffusion_sigma"));
    for (size_t i = 0; i < x.size(); ++i) {
        Sx.push_back(x[i] + sigma_xy * dist(gen));
        Sy.push_back(y[i] + sigma_xy * dist(gen));
        Sz.push_back(z[i]); // z remains unchanged
    }
}

std::vector<double> TrackProcessor::computeIonizationElectrons(const std::vector<double>& energyDep,
                                                               const std::vector<double>& z_hit)
{
    double ion_pot = std::stod(config.at("ion_pot"));
    std::vector<double> nel;
    for (size_t i = 0; i < energyDep.size(); ++i) {
        double drift = z_hit[i] + std::stod(config.at("z_extra"));
        double absorption = std::stod(config.at("absorption_l"));
        double electrons = (energyDep[i] / ion_pot) * std::exp(-drift / absorption);
        nel.push_back(std::max(0.0, electrons));
    }
    return nel;
}
