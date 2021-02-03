

#include "dxmc/attenuationinterpolator.h"
#include "dxmc/attenuationlut.h"
#include "xraylib.h"

#include <cassert>
#include <iostream>
#include <vector>

using namespace dxmc;

constexpr double ERRF = 1E-4;

template <Floating T>
bool isEqual(T a, T b)
{
    return std::abs(a - b) < ERRF;
}

template <Floating T>
bool testInterpolation()
{
    const T minEnergy = 1;
    const T maxEnergy = 150;

    std::vector<Material> materials;
    materials.emplace_back(("H2O"));
    materials.emplace_back(("H2O"));
    // materials.emplace_back((13));
    //materials.emplace_back((82));

    materials[0].setStandardDensity(1.0);
    materials[1].setStandardDensity(2.0);

    std::vector<T> dens;
    std::transform(materials.cbegin(), materials.cend(), std::back_inserter(dens), [](const auto& m) { return static_cast<T>(m.standardDensity()); });
    std::vector<std::uint8_t> mats(dens.size());
    std::iota(mats.begin(), mats.end(), 0);

    AttenuationLutInterpolator<T> attLut(materials, dens.cbegin(), dens.cend(), mats.cbegin(), T { 1 }, T { 150 });

    std::vector<T> energies;
    const std::size_t resolution = maxEnergy - minEnergy;
    for (std::size_t i = 0; i < resolution; ++i) {
        const T e = minEnergy + (i * (maxEnergy - minEnergy)) / (resolution - 1);
        energies.push_back(e);
    }

    for (const auto& mat : materials) {
        auto binding = mat.getBindingEnergies(minEnergy);
        for (const auto e : binding) {
            energies.push_back(e);
        }
    }
    std::sort(energies.begin(), energies.end());

    const auto material = 0;
    std::cout << "Energy, Total, Photo, Incoher, Coher, Total dxmc, Photo dxmc, Incoher dxmc, Coher dxmc, max total inv, " << materials[material].standardDensity() << "\n ";
    for (const auto e : energies) {
        std::cout << e << ", ";
        std::array<T, 3> xlib {
            static_cast<T>(CS_Photo_CP(materials[material].name().c_str(), e, nullptr)),
            static_cast<T>(CS_Compt_CP(materials[material].name().c_str(), e, nullptr)),
            static_cast<T>(CS_Rayl_CP(materials[material].name().c_str(), e, nullptr))
        };
        std::cout << std::reduce(xlib.cbegin(), xlib.cend()) << ", ";
        for (auto t : xlib) {
            std::cout << t << ", ";
        }
        const auto dx = attLut(material, e);
        const auto d_tot = std::reduce(dx.cbegin(), dx.cend());
        std::cout << d_tot << ", ";
        for (auto t : dx) {
            std::cout << t << ", ";
        }
        for (int i = 0; i < 3; ++i) {
            const auto valid = (dx[i] - xlib[i]) / xlib[i] * 1000;
            if (valid > 0.6)
                return false;
        }
        std::cout << attLut.maxAttenuationInverse(e) << ", ";
        std::cout << "\n";
    }

    return true;
}

int main(int argc, char* argv[])
{
    Material mat('K');

    assert(testInterpolation<float>());

    return EXIT_FAILURE;
}
