
#include "dxmc/material.h"
#include "dxmc/world.h"

#include <array>
#include <cstdlib>
#include <memory>
#include <vector>

template <typename T>
bool testWorld()
{
    dxmc::World<T> w;
    std::array<T, 3> spacing { 1, 1, 1 };
    std::array<std::size_t, 3> dim { 128, 128, 128 };
    w.setDimensions(dim);
    w.setSpacing(spacing);

    //dxmc::Material water("Water, Liquid");
    //w.addMaterialToMap(water);
    w.addMaterialToMap(dxmc::Material("Water, Liquid"));

    //arrays
    auto dens = std::make_shared<std::vector<T>>(w.size(), T { 1 });
    auto mat = std::make_shared<std::vector<std::uint8_t>>(w.size(), 0);
    w.setDensityArray(dens);
    w.setMaterialIndexArray(mat);
    w.makeValid();
    return w.isValid();
}

template <typename T>
bool testCTDIvol()
{
    dxmc::CTDIPhantom<T> w;
    return w.isValid();
}

int main(int argc, char* argv[])
{
    bool success = true;
    success = success && testWorld<float>();
    success = success && testWorld<double>();
    success = success && testCTDIvol<float>();
    success = success && testCTDIvol<double>();

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
