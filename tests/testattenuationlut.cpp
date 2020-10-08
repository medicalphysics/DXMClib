

#include "dxmc/attenuationlut.h"

#include <cassert>
#include <iostream>
#include <vector>

using namespace dxmc;

constexpr double ERRF = 1E-4;

bool isEqual(double a, double b)
{
    return std::abs(a - b) < ERRF;
}

bool testAttenuationLutXcom()
{
    std::vector<Material> mats;
    mats.emplace_back("H62.95C12.91N1.17O22.78Na0.03P0.04S0.06Cl0.03K0.03");
    mats.emplace_back(13);
    mats[0].setStandardDensity(1.03);

    AttenuationLut lut;
    lut.setEnergyResolution(1.0);
    

    std::vector<double> energies {
        1.0000, 10.0000, 15.0000, 20.0000, 30.0000, 40.0000, 50.0000, 54.6000, 60.0000, 80.0000, 100.0000, 150.0000, 200
    };

    lut.generate(mats, 1, 200);
    std::cout << "Test XCOM Soft tissue\n";
    std::cout << "Energy, Rayl, Compt, Photo, Total\n";
    for (const auto e : energies) {
        std::cout << e << ", ";
        std::cout << lut.rayleightAttenuation(0, e) << ", ";
        std::cout << lut.comptonAttenuation(0, e) << ", ";
        std::cout << lut.photoelectricAttenuation(0, e) << ", ";
        std::cout << lut.totalAttenuation(0, e) << "\n";
    }
    std::cout << "Test XCOM Aluminium\n";
    std::cout << "Energy, Rayl, Compt, Photo, Total\n";
    for (const auto e : energies) {
        std::cout << e << ", ";
        std::cout << lut.rayleightAttenuation(1, e) << ", ";
        std::cout << lut.comptonAttenuation(1, e) << ", ";
        std::cout << lut.photoelectricAttenuation(1, e) << ", ";
        std::cout << lut.totalAttenuation(1, e) << "\n";
    }



    return false;
}

bool testAttenuationLutTG195()
{
    constexpr double energy = 54.6;
    AttenuationLut lut;
    lut.setEnergyResolution(1.0);

    Material air("C0.02N78.44O21.08Ar0.47");
    air.setStandardDensity(0.001205);
    Material soft("H62.95C12.91N1.17O22.78Na0.03P0.04S0.06Cl0.03K0.03");
    soft.setStandardDensity(1.03);

    std::vector<Material> mats = { air, soft };

    lut.generate(mats, 1, 150);
    const auto att_diff = lut.photoComptRayAttenuation(1, energy);
    const auto att_tot = lut.totalAttenuation(1, energy);

    const double xcom_total = 2.124E-1;
    const std::array<double, 3> xcom_diff { 1.915E-02, 1.781E-01, 1.518E-02 };

    assert(std::abs(att_tot - lut.totalAttenuation(1, energy)) < 1E-7);
    assert(std::abs(att_diff[0] - lut.photoelectricAttenuation(1, energy)) < 1E-7);
    assert(std::abs(att_diff[1] - lut.comptonAttenuation(1, energy)) < 1E-7);
    assert(std::abs(att_diff[2] - lut.rayleightAttenuation(1, energy)) < 1E-7);

    std::cout << "Test attenuation LUT\n";
    std::cout << "attenuation coefficients for 54.6 keV photon in soft tissue:\n";
    std::cout << "Total, XCOM: " << xcom_total << " DXMC: " << att_tot << " difference [%]:" << (xcom_total - att_tot) / xcom_total * 100 << "\n";
    std::cout << "Photoelectric, XCOM: " << xcom_diff[0] << " DXMC: " << att_diff[0] << " difference [%]:" << (xcom_diff[0] - att_diff[0]) / xcom_diff[0] * 100 << "\n";
    std::cout << "Compton, XCOM: " << xcom_diff[1] << " DXMC: " << att_diff[1] << " difference [%]:" << (xcom_diff[1] - att_diff[1]) / xcom_diff[1] * 100 << "\n";
    std::cout << "Rayleight, XCOM: " << xcom_diff[2] << " DXMC: " << att_diff[2] << " difference [%]:" << (xcom_diff[2] - att_diff[2]) / xcom_diff[2] * 100 << "\n";

    const double xray_total = soft.getTotalAttenuation(energy);
    const std::array<double, 3> xray_diff { soft.getPhotoelectricAttenuation(energy), soft.getComptonAttenuation(energy), soft.getRayleightAttenuation(energy) };

    std::cout << "Total, xraylib: " << xray_total << " DXMC: " << att_tot << " difference [%]:" << (xray_total - att_tot) / xray_total * 100 << "\n";
    std::cout << "Photoelectric, xraylib: " << xray_diff[0] << " DXMC: " << att_diff[0] << " difference [%]:" << (xray_diff[0] - att_diff[0]) / xray_diff[0] * 100 << "\n";
    std::cout << "Compton, xraylib: " << xray_diff[1] << " DXMC: " << att_diff[1] << " difference [%]:" << (xray_diff[1] - att_diff[1]) / xray_diff[1] * 100 << "\n";
    std::cout << "Rayleight, xraylib: " << xray_diff[2] << " DXMC: " << att_diff[2] << " difference [%]:" << (xray_diff[2] - att_diff[2]) / xray_diff[2] * 100 << "\n";

    std::cout << "\n";

    bool success = isEqual(att_tot, xcom_total);
    for (std::size_t i = 0; i < 3; ++i)
        success = success && isEqual(att_diff[i], xcom_diff[i]);

    success = success && isEqual(att_tot, xray_total);
    for (std::size_t i = 0; i < 3; ++i)
        success = success && isEqual(att_diff[i], xray_diff[i]);
    assert(success);
    return success;
}

int main(int argc, char* argv[])
{
    testAttenuationLutXcom();
    auto success = testAttenuationLutTG195();
    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
