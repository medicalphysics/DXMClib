

#include "dxmc/constants.h"
#include "dxmc/dxmcrandom.h"
#include "dxmc/interpolation.h"
#include "dxmc/material.h"
#include "dxmc/attenuationinterpolator.h"

#include "xraylib.h"

#include <assert.h>
#include <iostream>
bool testSpline()
{

    std::vector<float> x(30);
    std::vector<float> y(x.size());
    for (int i = 0; i < x.size(); ++i) {
        x[i] = i * (360.0 / x.size()) * dxmc::DEG_TO_RAD<float>();
        y[i] = std::sin(x[i]);
    }
    dxmc::CubicSplineInterpolator S(x, y);
    float diff = 0;
    for (int i = 0; i < x.size() - 1; ++i) {
        std::cout << x[i] << ", ";
        std::cout << y[i] << ", ";
        std::cout << (x[i] + x[i + 1]) / 2 << ", ";
        std::cout << S((x[i] + x[i + 1]) / 2) << ", ";
        std::cout << std::sin((x[i] + x[i + 1]) / 2) << "\n";
        const auto d = S((x[i] + x[i + 1]) / 2) - std::sin((x[i] + x[i + 1]) / 2);
        diff += d * d;
    }

    const auto rms = std::sqrt(diff) / x.size();
    return rms < 0.001;
}

template <typename T>
bool testSplineAttenuation()
{
    //testing spline interpolation for photoelectric effect

    dxmc::Material mat(82);
    const T minEnergy = T { 1 };
    const T maxEnergy = 150;

    auto par = [&](const double e) { return static_cast<T>(mat.getComptonAttenuation(e)); }; // mat.getPhotoelectricAttenuation(e)); 
       

    const auto bindingEnergies = mat.getBindingEnergies(minEnergy);

    const std::size_t resolution = 150;

    std::vector<T> energies = bindingEnergies;
    for (std::size_t i = 0; i < resolution; ++i) {
        const T e = minEnergy + (i * (maxEnergy - minEnergy)) / (resolution - 1);
        energies.push_back(e);
    }
    std::sort(energies.begin(), energies.end());
    std::vector<T> att(energies.size());
    std::transform(energies.cbegin(), energies.cend(), att.begin(), [&](const T e) -> T { return par(e); });

    //spline

    const std::size_t spline_resolution = 30-bindingEnergies.size()*2;
    const T splineEdgeOffset = T { 0.001 };
    std::vector<T> spline_energies;
    const T logmin = std::log10(minEnergy);
    const T logmax = std::log10(maxEnergy);
    for (std::size_t i = 0; i < spline_resolution; ++i) {
        const T e = logmin + (i * (logmax - logmin)) / (spline_resolution - 1);
        spline_energies.push_back(e);
    }



    for (const T e : bindingEnergies) {

        const T val1 = std::log10(e - splineEdgeOffset);
        const T val2 = std::log10(e);
        
        spline_energies.push_back(val1);
        spline_energies.push_back(val2);

    }

    std::sort(spline_energies.begin(), spline_energies.end());
    std::vector<T> spline_att(spline_energies.size());
    for (std::size_t i = 0; i < spline_energies.size(); ++i) {
        const auto e = std::pow(10, spline_energies[i]);
        const auto val = par(e);
        const auto res = std::log10(val);
        spline_att[i] = std::log10(val);
    }

    dxmc::CubicSplineInterpolator spline(spline_energies, spline_att);

    std::vector<dxmc::Material> materials(1);
    materials[0] = mat;

    dxmc::AttenuationLutInterpolator<T> attLut(materials, T { 1 }, T { 150 });




    //printing
    std::cout.precision(std::numeric_limits<T>::digits10 + 1);
    std::cout << "Energy, basic, spline\n";
    for (std::size_t i = 0; i < energies.size(); ++i) {
        std::cout << energies[i] << ", ";
        std::cout << att[i] << ", ";
        std::cout << std::pow(10, spline(std::log10(energies[i])));
        if (i < spline_energies.size()) {
            std::cout << ", " << std::pow(10, spline_energies[i]) << ", ";
            std::cout << std::pow(10, spline_att[i]);
        }
        std::cout << std::endl;
    }

    return false;
}

int main()
{

    // bool splinetest = testSpline();
    bool photo = testSplineAttenuation<float>();
    //assert(splinetest);
    return 1;
}
