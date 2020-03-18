


#include "dxmc/transport.h"
#include "xraylib.h"
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
        const auto ind = static_cast<std::size_t>((cosang + 1.0) * nHist*0.5);
        hist[ind] += 1;
    }

    std::vector<std::size_t> histG(nHist, 0);
    for (std::size_t i = 0; i < samples; ++i)
    {
        transport::comptonScatterGeant(p, seed, cosang);
        const auto ind = static_cast<std::size_t>((cosang + 1.0) * nHist * 0.5);
        histG[ind] += 1;
    }
    std::vector<std::size_t> histE(nHist, 0);
    for (std::size_t i = 0; i < samples; ++i)
    {
        transport::comptonScatterEGS(p, seed, cosang);
        const auto ind = static_cast<std::size_t>((cosang + 1.0) * nHist * 0.5);
        histE[ind] += 1;
    }
    
    Material mat("Water, Liquid");

    std::vector<double> histF(nHist, 0);
    for (std::size_t i = 0; i < nHist; ++i)
    {
        cosang = (2.0 * i) / nHist - 1.0;
        double theta = std::acos(cosang);
        const auto ind = static_cast<std::size_t>(i);
        double att = DCS_Compt_CP((mat.name()).c_str(), keV, theta, nullptr);
        histF[ind] = att;
    }




    std::cout << "cos angle, hist, hist Geant, hist EGS, analytical\n";
    for (int i = 0; i < nHist; ++i)
    {
        std::cout << 2.0 * i / nHist - 1.0 << ", " << hist[i] << ", " << histG[i] << ", " << histE[i] << ", " << histF[i] << "\n";
    }
    std::cout << std::endl;


    return true;


}

bool testRayleight(double keV, const Material& mat)
{
    constexpr std::size_t nHist = 100;
    constexpr std::size_t samples = 1e5;

    std::uint64_t seed[2];
    randomSeed(seed);

    std::vector<std::size_t> hist(nHist, 0);


    AttenuationLut lut;
    std::vector<Material> mats;
    mats.push_back(mat);
    lut.generate(mats);



    Particle p;
    p.energy = keV;
    double cosang;
    for (std::size_t i = 0; i < samples; ++i)
    {
        transport::rayleightScatter(p, 0, lut, seed, cosang);
        auto ind = static_cast<std::size_t>((cosang + 1.0) * (nHist) * 0.5);
        ind = ind >= nHist ? nHist - 1 : ind;
        hist[ind] += 1;
    }

    std::vector<double> ana(nHist, 0.0);
    for (std::size_t i = 0; i < nHist; ++i)
    {
        cosang = (2.0 * i) / nHist - 1.0;
        double theta = std::acos(cosang);
        const auto ind = static_cast<std::size_t>(i);
        
        double att = DCS_Rayl_CP(mat.name().c_str(), keV, theta, nullptr);
        ana[ind] = att;
    }

    auto hist_sum = std::reduce(hist.cbegin(), hist.cend(), 0.0);

    auto ana_sum = std::reduce(ana.cbegin(), ana.cend(), 0.0);

    std::cout << "cos angle, Rayleight, Analytical\n";
    for (int i = 0; i < nHist; ++i)
    {
        std::cout << 2.0 * i / nHist - 1.0 << ", " << hist[i]/hist_sum<<", "<< ana[i]/ana_sum <<"\n";
    }
    std::cout << std::endl;


    return true;





}

int main (int argc, char *argv[])
{

    //testCompton(30.0);
    //Material mat("Water, Liquid");
    Material mat(1);
    testRayleight(10, mat);
    return 0;
}

