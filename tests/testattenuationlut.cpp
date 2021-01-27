

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
bool testInterpolation(const Material& mat)
{

    std::vector<Material> mats { mat };

    AttenuationLut<T> lut;
    lut.generate(mats);

    std::vector<T> energies(lut.energyBegin(), lut.energyEnd());
    for (auto& e : energies) {
        e += (energies[1] - energies[0]) / 2;
    }

    std::cout << "Energy, Photo, Incoher, Coher, Photo dxmc, Incoher dxmc, Coher dxmc\n";

    for (const auto e : energies) {
        std::cout << e << ", ";

        std::cout << CS_Photo_CP(mat.name().c_str(), e, nullptr) << ", ";
        std::cout << CS_Compt_CP(mat.name().c_str(), e, nullptr) << ", ";
        std::cout << CS_Rayl_CP(mat.name().c_str(), e, nullptr) << ", ";
        const auto abc = lut.photoComptRayAttenuation(0, e);
        for (const auto a : abc) {
            std::cout << a << ", ";
        }
        std::cout << "\n";
    }
   
    for (const auto& el : lut.electronShellConfiguration(0)) {
        std::cout << el.bindingEnergy << ", "<< el.numberElectrons << "\n";
    }

    return true;
}

int main(int argc, char* argv[])
{
    Material mat('K');
    assert(testInterpolation<float>(mat));

    return EXIT_FAILURE;
}
