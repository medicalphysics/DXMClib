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
#include "dxmc/world/worlditems/ctdiphantom.hpp"
#include "dxmc/world/worlditems/worldbox.hpp"

#include <iostream>

template <typename T>
bool testTransport()
{
    dxmc::CTDIPhantom<T> phantom;
    dxmc::WorldBox<T> box;
    box.translate({ 0, -50, 0 });
    box.setNistMaterial("Water, Liquid");

    dxmc::World2<T, dxmc::CTDIPhantom<T>, dxmc::WorldBox<T>> world;
    world.addItem(phantom);
    world.addItem(box);
    world.build();

    dxmc::PencilBeam<T> beam;
    beam.setPosition({ 0, -1000, 0 });
    beam.setDirection({ 0, 1, 0 });
    beam.setNumberOfParticlesPerExposure(100);
    beam.setNumberOfExposures(4);

    dxmc::Transport<T> transport;
    transport.setNumberOfThreads(1);

    transport(world, beam);

    const auto& ctdi_vec = world.template getItems<dxmc::CTDIPhantom<T>>();
    auto dose_ctdi = ctdi_vec[0].dose(0);
    std::cout << "CTDI: " << dose_ctdi.energyImparted() << std::endl;

    const auto& box_vec = world.template getItems<dxmc::WorldBox<T>>();
    auto dose_box = box_vec[0].dose(0);
    std::cout << "Box: " << dose_box.energyImparted() << std::endl;

    return true;
}

int main()
{
    bool success = true;
    success = success && testTransport<float>();

    if (success)
        return 0;

    return 1;
}