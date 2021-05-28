

#include "dxmc/material.hpp"
#include "dxmc/transport.hpp"

#include <cassert>
#include <iostream>
#include <numeric>

using namespace dxmc;

bool testMaterialParser()
{
    auto materialNISTList = Material::getNISTCompoundNames();

    std::vector<Material> materials;

    materials.push_back(Material("Urea"));
    materials.push_back(Material(2));
    materials.push_back(Material("H2O"));
    materials[2].setStandardDensity(1.0);
    materials.push_back(Material("Bone, Compact (ICRU)"));

    std::map<std::size_t, Material> mapping;
    for (std::size_t i = 0; i < materials.size(); i++)
        mapping[i] = materials[i];

    std::cout << "Test material parser\n";

    for (auto& m : materials) {
        std::cout << m.name() << " Valid: " << m.isValid() << "\n";        
    }

    bool success = true;
    for (auto& m : materials)
        success = success && m.isValid();
    assert(success);
    return success;
}

int main(int argc, char* argv[])
{
    auto success = testMaterialParser();
    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
