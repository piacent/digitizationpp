#pragma once

#include "ConfigManager.h"
#include "TH2F.h"
#include <vector>
#include <string>
#include <map>

class TrackProcessor {
public:
    TrackProcessor(ConfigManager& config);

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


    void TrackVignetting(std::vector<std::vector<double>>& image, int xpix, int ypix, const TH2F& VignMap);

private:
    ConfigManager& config;

    void cloud_smearing3D(  const std::vector<double>& x_hits_tr,
                            const std::vector<double>& y_hits_tr,
                            const std::vector<double>& z_hits_tr,
                            const std::vector<double>& energy_hits,
                            std::vector<float>& S3D_x,
                            std::vector<float>& S3D_y,
                            std::vector<float>& S3D_z);

    std::vector<float> smear(const std::vector<double>& axis_hit,
                             const std::vector<double>& axis_sigma,
                             const std::vector<double>& nel);

    void smear_parallel(const std::vector<double>& x_axis_hit,
                        const std::vector<double>& y_axis_hit,
                        const std::vector<double>& z_axis_hit,
                        const std::vector<double>& x_axis_sigma,
                        const std::vector<double>& y_axis_sigma,
                        const std::vector<double>& z_axis_sigma,
                        const std::vector<double>& nel,
                        std::vector<float>& X,
                        std::vector<float>& Y,
                        std::vector<float>& Z);

    void ph_smearing2D( const std::vector<double>& x_hits_tr,
                        const std::vector<double>& y_hits_tr,
                        const std::vector<double>& z_hits_tr,
                        const std::vector<double>& energy_hits,
                        std::vector<float>& S2D_x,
                        std::vector<float>& S2D_y);

    std::vector<double> compute_sigma(const double diff_const, const double diff_coeff, const std::vector<double>& dz);

    std::vector<double> NelGEM1(const std::vector<double>& N_ioniz_el);
    std::vector<double> NelGEM2(const std::vector<double>& energyDep,const std::vector<double>& z_hit);

    double Nph_saturation(int nel, double A, double beta);



};
