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

#include "dxmc/dxmcrandom.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/visualization/visualizeworld.hpp"
#include "dxmc/world/world.hpp"
#include "dxmc/world/worlditems/ctdiphantom.hpp"
#include "dxmc/world/worlditems/depthdose.hpp"
#include "dxmc/world/worlditems/worldbox.hpp"
#include "dxmc/world/worlditems/worldsphere.hpp"

#include <iostream>

template <dxmc::Floating T>
bool testWorld()
{
    using Sphere = dxmc::WorldSphere<T, 5, 2>;
    using CTDI = dxmc::CTDIPhantom<T, 5, 2>;
    using World = dxmc::World<T, Sphere, CTDI>;

    World world;
    world.reserveNumberOfItems(2);
    auto& sphere = world.addItem<Sphere>();
    auto& ctdi = world.addItem<CTDI>();

    sphere.translate({ 0, 0, 32 });
    sphere.setRadius(16);

    world.build();
    world.clearDose();
    dxmc::VisualizeWorld<T> viz(world);

    viz.setDistance(150);
    viz.setPolarAngleDeg(180);
    viz.setAzimuthalAngleDeg(0);
    viz.suggestFOV();
    viz.setCameraPosition({ 100, 100, 100 });

    std::size_t size = 512;
    std::vector<std::uint8_t> buffer(size * size * 4, 255);

    viz.generate(world, buffer, size, size);
    viz.savePNG("test.png", buffer, size, size);

    bool test = false;
    return test;
}

template <typename T, typename WO>
bool testWorldItem()
{

    WO obj;

    dxmc::Particle<T> p;
    p.dir = { 0, 0, 1 };
    p.pos = { 0, 0, -2000 };

    dxmc::vectormath::normalize(p.dir);

    auto t = obj.intersect(p);
    if (!t.valid())
        return false;
    p.border_translate(t.intersection);

    dxmc::RandomState state;
    p.energy = 60;
    p.weight = 1;
    obj.transport(p, state);

    auto t_fin = obj.intersect(p);
    if (t_fin.valid())
        return false;

    return true;
}

template <typename T>
bool testWorldItems()
{
    bool success = true;

    success = success && testWorldItem<T, dxmc::CTDIPhantom<T>>();
    if (success)
        std::cout << "SUCCESS Basic test for CTDIPhantom<T> for sizeof(T) = " << sizeof(T) << std::endl;
    else
        std::cout << "FAILURE Basic test for CTDIPhantom<T> for sizeof(T) = " << sizeof(T) << std::endl;

    success = success && testWorldItem<T, dxmc::WorldBox<T>>();
    if (success)
        std::cout << "SUCCESS Basic test for WorldBox<T> for sizeof(T) = " << sizeof(T) << std::endl;
    else
        std::cout << "FAILURE Basic test for WorldBox<T> for sizeof(T) = " << sizeof(T) << std::endl;

    success = success && testWorldItem<T, dxmc::WorldSphere<T>>();
    if (success)
        std::cout << "SUCCESS Basic test for WorldSphere<T> for sizeof(T) = " << sizeof(T) << std::endl;
    else
        std::cout << "FAILURE Basic test for WorldSphere<T> for sizeof(T) = " << sizeof(T) << std::endl;
    return success;
}

template <typename T>
bool testBorderCrossing()
{

    dxmc::World<T, dxmc::WorldBox<T>> w;

    dxmc::WorldBox<T> b1, b2;
    b2.translate({ 0, 0, 2 });
    w.addItem(b1);
    w.addItem(b2);
    w.build();

    dxmc::Particle<T> p;
    p.dir = { 0, 0, 1 };
    p.pos = { 0, 0, -2000 };
    dxmc::vectormath::normalize(p.dir);

    dxmc::RandomState state;
    p.energy = 60;
    p.weight = 1;
    w.transport(p, state);

    return true;
}

int main(int argc, char* argv[])
{
    std::cout << "World Tests\n";
    auto success = true;

    success = success && testWorld<float>();
    success = success && testWorld<double>();

    success = success && testBorderCrossing<float>();
    success = success && testBorderCrossing<double>();

    // success = success && testTransport<float>();
    // success = success && testTransport<double>();

    success = success && testWorldItems<double>();
    success = success && testWorldItems<float>();

    if (success)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}
