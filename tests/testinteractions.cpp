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

#include "dxmc/interactions.hpp"
#include "dxmc/material/material.hpp"

#include <iostream>

template <typename T>
bool testCoherent(std::size_t Z = 13, T energy = 50)
{
    bool success = true;

    auto material_opt = dxmc::Material2<T>::byZ(Z);
    auto material = material_opt.value();

    dxmc::Particle<T> particle;
    particle.pos = { 0, 0, 0 };
    particle.dir = { 0, 0, 1 };
    particle.energy = energy;
    particle.weight = 1;

    dxmc::RandomState state;

    dxmc::interactions::rayleightScatter<T, 1>(particle, material, state);

    // distribution test goes here, compare to xraylib ? ? ?

    return success;
}

int main()
{
    std::cout << "Testing interactions" << std::endl;
    bool success = true;
    success = success && testCoherent<double>();
    success = success && testCoherent<float>();
    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
