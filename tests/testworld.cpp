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

Copyright 2022 Erlend Andersen
*/

#include "dxmc/vectormath.hpp"
#include "dxmc/world/box.hpp"
#include "dxmc/world/ctdiphantom.hpp"
#include "dxmc/world/sphere.hpp"
#include "dxmc/world/world.hpp"
#include "dxmc/dxmcrandom.hpp"

#include <chrono>
#include <fstream>
#include <iostream>

template <dxmc::Floating T>
bool testWorld()
{
    dxmc::World2<T, dxmc::Sphere<T>, dxmc::CTDIPhantom<T>, dxmc::Box<T>> world;

    dxmc::CTDIPhantom<T> ctdi;
    dxmc::Sphere<T> sphere(32, { 50, 50, 50 });
    dxmc::Box<T> box;
    box.translate({ 50, 0, -50 });

    world.addItem(ctdi);
    world.addItem(sphere);
    world.addItem(box);

    world.build();

    std::array<T, 3> trans { 5, 5, 5 };
    

    dxmc::Particle<T> p {
        .pos = { 0, 0, -100 }, .dir = { 0, 0, 1 }
    };

    auto res = world.intersect(p);

    dxmc::RandomState rand;
    world.transport(p, rand);

    return false;
}

int main(int argc, char* argv[])
{
    auto success = true;
    success = success && testWorld<double>();
    success = success && testWorld<float>();

    if (success)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}
