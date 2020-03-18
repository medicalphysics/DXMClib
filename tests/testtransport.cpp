


#include "dxmc/transport.h"
#include "xraylib.h"
#include <iostream>

bool testCompton(double keV, const Material& mat)
{
    constexpr std::size_t nHist = 100;
    constexpr std::size_t samples = 1e5;

    std::uint64_t seed[2];
    randomSeed(seed);
    
 
    AttenuationLut lut;
    std::vector<Material> mats;
    mats.push_back(mat);
    lut.generate(mats, 0.,keV);
    lut.generate(mats);

    Particle p;
    p.energy = keV;
    double cosang;
    std::vector<std::size_t> histL(nHist, 0);
    for (std::size_t i = 0; i < samples; ++i)
    {
        transport::comptonScatterLivermore(p, 0, lut, seed, cosang);
        auto ind = static_cast<std::size_t>((cosang + 1.0) * (nHist) * 0.5);
        histL[ind] += 1;
    }
    const double histL_sum = std::accumulate(histL.cbegin(), histL.cend(), 0.0);

    std::vector<std::size_t> hist(nHist, 0);
    for (std::size_t i = 0; i < samples; ++i)
    {
        transport::comptonScatter(p, seed, cosang);
        auto ind = static_cast<std::size_t>((cosang + 1.0) * (nHist) * 0.5);
        hist[ind] += 1;
    }
    const double hist_sum = std::accumulate(hist.cbegin(), hist.cend(), 0.0);


    std::vector<double> histF(nHist, 0);
    for (std::size_t i = 0; i < nHist; ++i)
    {
        cosang = (2.0 * i) / nHist - 1.0;
        double theta = std::acos(cosang);
        const auto ind = static_cast<std::size_t>(i);
        double att = DCS_Compt_CP((mat.name()).c_str(), keV, theta, nullptr);
        histF[ind] = att;
    }
    const double histF_sum = std::accumulate(histF.cbegin(), histF.cend(), 0.0);



    std::cout << "cos angle, hist Livermore, hist, analytical\n";
    for (int i = 0; i < nHist; ++i)
    {
        std::cout << 2.0 * i / nHist - 1.0 << ", " <<histL[i]/histL_sum << ", "<< hist[i]/hist_sum << ", " << histF[i]/histF_sum << "\n";
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

    AttenuationLut lut;
    std::vector<Material> mats;
    mats.push_back(mat);
    lut.generate(mats);



    Particle p;
    p.energy = keV;
    double cosang;

    std::vector<std::size_t> histL(nHist, 0);
    for (std::size_t i = 0; i < samples; ++i)
    {
        transport::rayleightScatterLivermore(p, 0, lut, seed, cosang);
        auto ind = static_cast<std::size_t>((cosang + 1.0) * (nHist) * 0.5);
        ind = ind >= nHist ? nHist - 1 : ind;
        histL[ind] += 1;
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

    auto histL_sum = std::reduce(histL.cbegin(), histL.cend(), 0.0);
    auto ana_sum = std::reduce(ana.cbegin(), ana.cend(), 0.0);

    std::cout << "cos angle,  RayleightLivermore, Analytical\n";
    for (int i = 0; i < nHist; ++i)
    {
        std::cout << 2.0 * i / nHist - 1.0 <<  ", " << histL[i] / histL_sum << ", " << ana[i] / ana_sum << "\n";
    }
    std::cout << std::endl;


    return true;

}

int main (int argc, char *argv[])
{

    //testCompton(30.0);
    //Material mat("Water, Liquid");
    Material mat(6);
    //Material mat(8);
    testRayleight(50, mat);
    //testCompton(10, mat);
    return 0;
}

