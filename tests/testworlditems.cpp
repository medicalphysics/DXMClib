/*This file is part of DXMClib.

DXMClib is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

DXMClib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with DXMClib. If not, see < https://www.gnu.org/licenses/>.

Copyright 2023 Erlend Andersen
*/

#include "dxmc/world/worlditems/aavoxelgrid.hpp"
#include "dxmc/world/worlditems/ctdiphantom.hpp"

#include <chrono>
#include <fstream>
#include <iostream>

template <typename T>
bool testAAVoxelGrid()
{
    dxmc::AAVoxelGrid<T, 5, 2> item;

    auto air = dxmc::Material<T, 5>::byNistName("Air, Dry (near sea level)").value();
    auto pmma = dxmc::Material<T, 5>::byNistName("Polymethyl Methacralate (Lucite, Perspex)").value();
    const auto air_dens = dxmc::NISTMaterials<T>::density("Air, Dry (near sea level)");
    const auto pmma_dens = dxmc::NISTMaterials<T>::density("Polymethyl Methacralate (Lucite, Perspex)");
    auto lead = dxmc::Material<T, 5>::byZ(82).value();
    const auto lead_dens = dxmc::AtomHandler<T>::Atom(82).standardDensity;

    const std::array<std::size_t, 3> dim = { 13, 13, 13 };
    const auto size = std::reduce(dim.cbegin(), dim.cend(), size_t { 1 }, std::multiplies<>());
    std::array<T, 3> spacing = { 2.0f, 2.0f, 2.0f };

    // material arrays
    std::vector<T> dens(size, air_dens);
    std::vector<std::uint8_t> materialIdx(size, 0);
    std::vector<dxmc::Material<T, 5>> materials;
    materials.push_back(air);
    materials.push_back(pmma);
    materials.push_back(lead);

    materialIdx[1] = 1;
    dens[1] = pmma_dens;
    materialIdx[2] = 2;
    dens[2] = lead_dens;

    item.setData(dim, dens, materialIdx, materials);
    item.setSpacing(spacing);

    bool success = true;

    T energy = dxmc::MIN_ENERGY<T>() * 2;
    while (energy < dxmc::MAX_ENERGY<T>()) {
        const auto att = lead.attenuationValues(energy).sum() * lead_dens;
        const auto attmax = item.maxAttenuationValue(energy);
        success = success && attmax >= att;
        // std::cout << energy << ", " << att << ", " << attmax << std::endl;
        energy += T { 2 };
    }

    return success;
}

template <dxmc::Floating T>
bool testCTDIPhantom() {
    dxmc::CTDIPhantom<T, 5, 2> phantom;

}

int main(int argc, char* argv[])
{
    auto success = true;

    success = success && testAAVoxelGrid<float>();
    success = success && testAAVoxelGrid<double>();

    if (success)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}
