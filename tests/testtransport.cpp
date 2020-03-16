


#include "dxmc/transport.h"
#include <iostream>

bool testCompton(double keV)
{
    constexpr std::size_t nHist = 100;
    constexpr std::size_t samples = 1e6;

    std::uint64_t seed[2];
    randomSeed(seed);
    
    std::vector<std::size_t> hist(nHist, 0);

    Particle p;
    p.energy = keV;
    double cosang;
    for (std::size_t i = 0; i < samples; ++i)
    {
        transport::comptonScatter(p, seed, cosang);
        const auto ind = static_cast<std::size_t>((cosang + 1.0) * 50.0);
        hist[ind] += 1;
    }

    std::vector<std::size_t> histG(nHist, 0);
    for (std::size_t i = 0; i < samples; ++i)
    {
        transport::comptonScatterGeant(p, seed, cosang);
        const auto ind = static_cast<std::size_t>((cosang + 1.0) * 50.0);
        histG[ind] += 1;
    }
    std::vector<std::size_t> histE(nHist, 0);
    for (std::size_t i = 0; i < samples; ++i)
    {
        transport::comptonScatterEGS(p, seed, cosang);
        const auto ind = static_cast<std::size_t>((cosang + 1.0) * 50.0);
        histE[ind] += 1;
    }


    AttenuationLut lut;
    std::vector<Material> mats;
    mats.emplace_back(Material("Water, Liquid"));
    lut.generate(mats);

    std::vector<double> histF(nHist, 0);
    for (std::size_t i = 0; i < nHist; ++i)
    {
        double theta = std::acos(i / 50.0 - 1.0);
        double att =  DCSb_Compt_CP(mats[0].name().c_str(), keV, theta, nullptr);
        
        histRay[i] = att;
    }


    std::cout << "cos angle, hist, hist Geant, hist EGS, rayleight\n";
    for (int i = 0; i < nHist; ++i)
    {
        std::cout << -1 + i / 50.0 << ", " << hist[i] << ", " << histG[i] << ", " << histE[i] << ", " << histRay[i] << "\n";
    }
    std::cout << std::endl;


    return true;


}



int main (int argc, char *argv[])
{

    testCompton(60.0);
    return 0;
}

