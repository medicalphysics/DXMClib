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

#include "dxmc/beams/pencilbeam.hpp"
#include "dxmc/transport.hpp"
#include "dxmc/world/world.hpp"
#include "dxmc/world/worlditems/worldsphere.hpp"

#include <iostream>

template <dxmc::Floating T>
bool testForcedinteractions()
{
    using Sphere = dxmc::WorldSphere<T, 5, 1, false>;
    using SphereF = dxmc::WorldSphere<T, 5, 1, true>;

    using World = dxmc::World<T, Sphere>;
    using WorldF = dxmc::World<T, SphereF>;

    World w;
    WorldF wf;

    constexpr T radii = 5;
    auto& sphere = w.addItem<Sphere>({ radii });
    auto& spheref = wf.addItem<SphereF>({ radii });

    auto material_water = dxmc::Material<T, 5>::byNistName("Water, Liquid").value();
    sphere.setMaterial(material_water, 1);
    spheref.setMaterial(material_water, 1);

    w.build();
    wf.build();

    dxmc::PencilBeam<T> beam({ -100, 0, 0 }, { 1, 0, 0 }, 20);
    beam.setNumberOfExposures(124);
    beam.setNumberOfParticlesPerExposure(1e4);

    dxmc::Transport transport;
    transport(wf, beam);
    transport(w, beam);

    auto dose = sphere.doseScored().dose();
    auto dosef = spheref.doseScored().dose();

    auto diff = std::abs(dose - dosef);

    auto sttd = sphere.doseScored().standardDeviation();
    auto sttdf = spheref.doseScored().standardDeviation();

    std::cout << "Random " << sphere.doseScored().dose() << std::endl;
    std::cout << "Forced " << spheref.doseScored().dose() << std::endl;
    std::cout << "Diff: " << diff << ", sttdev: ";
    std::cout << sttd << ", " << sttdf << std::endl;

    return false;
}

template <dxmc::Floating T>
bool testForcedinteractionsExits()
{
    using Sphere = dxmc::WorldSphere<T, 5, 1, false>;
    using SphereF = dxmc::WorldSphere<T, 5, 1, true>;

    using World = dxmc::World<T, Sphere>;
    using WorldF = dxmc::World<T, SphereF>;

    World w(2);
    WorldF wf(2);


    constexpr T radii = 1;
    auto& sphere1 = w.addItem<Sphere>({ radii });
    auto& sphere2 = w.addItem<Sphere>({ radii, { radii * 2 + 1, 0, 0 } });
    auto& spheref1 = wf.addItem<SphereF>({ radii });
    auto& spheref2 = wf.addItem<SphereF>({ radii, { radii * 2 + 1, 0, 0 } });

    auto material_water = dxmc::Material<T, 5>::byNistName("Water, Liquid").value();
    sphere1.setMaterial(material_water, 1);
    sphere2.setMaterial(material_water, 1);
    spheref1.setMaterial(material_water, 1);
    spheref2.setMaterial(material_water, 1);

    w.build();
    wf.build();


    dxmc::PencilBeam<T> beam({ -100, 0, 0 }, { 1, 0, 0 }, 20);
    beam.setNumberOfExposures(48);
    beam.setNumberOfParticlesPerExposure(1e4);

    dxmc::Transport transport;
    transport(wf, beam);
    transport(w, beam);

    auto dose1 = sphere1.doseScored().dose();
    auto dose2 = sphere2.doseScored().dose();
    auto dosef1 = spheref1.doseScored().dose();
    auto dosef2 = spheref2.doseScored().dose();

    auto diff1 = std::abs(dose1 - dosef1);
    auto diff2 = std::abs(dose2 - dosef2);

    auto sttd1 = sphere1.doseScored().standardDeviation();
    auto sttdf1 = spheref1.doseScored().standardDeviation();

    auto sttd2 = sphere2.doseScored().standardDeviation();
    auto sttdf2 = spheref2.doseScored().standardDeviation();

    std::cout << "Random " << dose1 << ", " << dose2 << " stddev: " << sttd1 << ", " << sttd2 << std::endl;
    std::cout << "Forced " << dosef1 << ", " << dosef2 << " stddev: " << sttdf1 << ", " << sttdf2 << std::endl;
    std::cout << "Diff: " << diff1 << ", " << diff2 << std::endl;

    return false;
}

int main(int argc, char* argv[])
{
    auto success = true;

    success = success && testForcedinteractionsExits<double>();
    //success = success && testForcedinteractions<float>();

    if (success)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}
