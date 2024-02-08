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

    constexpr T radii = 10;
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

    std::cout << sphere.doseScored().dose() << std::endl;
    std::cout << spheref.doseScored().dose() << std::endl;
    std::cout << "Diff: " << diff << ", sttdev: ";
    std::cout << sttd << ", " << sttdf << std::endl;

    return false;
}

int main(int argc, char* argv[])
{
    auto success = true;

    success = success && testForcedinteractions<double>();
    success = success && testForcedinteractions<float>();

    if (success)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}
