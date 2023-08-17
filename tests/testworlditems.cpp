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
#include "dxmc/world/worlditems/aavoxelgrid.hpp"
#include "dxmc/world/worlditems/ctdiphantom.hpp"
#include "dxmc/world/worlditems/depthdose.hpp"
#include "dxmc/world/worlditems/fluencescore.hpp"
#include "dxmc/world/worlditems/tetrahedalmesh.hpp"
#include "dxmc/world/worlditems/triangulatedmesh.hpp"
#include "dxmc/world/worlditems/worldbox.hpp"
#include "dxmc/world/worlditems/worldboxgrid.hpp"
#include "dxmc/world/worlditems/worldcylinder.hpp"
#include "dxmc/world/worlditems/worldsphere.hpp"

template <dxmc::Floating T, typename U>
    requires std::is_base_of<dxmc::WorldItemBase<T>, U>::value
bool testItem()
{
    dxmc::World<T, U> world;
    world.reserveNumberOfItems(1);
    auto& item = world.template addItem<U>({});
    world.build();

    dxmc::PencilBeam<T> beam;
    beam.setNumberOfExposures(1);
    beam.setNumberOfParticlesPerExposure(8);

    dxmc::Transport transport;
    transport.setNumberOfThreads(1);

    transport(world, beam);

    return true;
}

template <dxmc::Floating T>
bool basicTestAllItems()
{
    auto success = true;
    success = success && testItem<T, dxmc::AAVoxelGrid<T>>();
    success = success && testItem<T, dxmc::CTDIPhantom<T>>();
    success = success && testItem<T, dxmc::DepthDose<T>>();
    success = success && testItem<T, dxmc::FluenceScore<T>>();
    success = success && testItem<T, dxmc::TetrahedalMesh<T>>();
    success = success && testItem<T, dxmc::TriangulatedMesh<T>>();
    success = success && testItem<T, dxmc::WorldBox<T>>();
    success = success && testItem<T, dxmc::WorldBoxGrid<T>>();
    success = success && testItem<T, dxmc::WorldCylinder<T>>();
    success = success && testItem<T, dxmc::WorldSphere<T>>();

    return success;
}

int main(int argc, char* argv[])
{
    auto success = true;

    success = success && basicTestAllItems<float>();
    success = success && basicTestAllItems<double>();

    if (success)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}
