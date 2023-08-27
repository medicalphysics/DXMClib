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

#include "dxmc/beams/dxbeam.hpp"
#include "dxmc/transport.hpp"
#include "dxmc/world/visualization/visualizeworld.hpp"
#include "dxmc/world/world.hpp"
#include "dxmc/world/worlditems/aavoxelgrid.hpp"
#include "dxmc/world/worlditems/ctdiphantom.hpp"
#include "dxmc/world/worlditems/triangulatedmesh.hpp"
#include "dxmc/world/worlditems/worldsphere.hpp"

#include <vector>

int main()
{
    using CTDIPhantom = dxmc::CTDIPhantom<double, 5, 1>;
    using Mesh = dxmc::TriangulatedMesh<double, 5, 1>;
    using Sphere = dxmc::WorldSphere<double, 5, 1>;
    using VGrid = dxmc::AAVoxelGrid<double, 5, 1>;
    using World = dxmc::World<double, CTDIPhantom, Mesh, Sphere, VGrid>;

    using Viz = dxmc::VisualizeWorld<double>;

    World world {};
    world.reserveNumberOfItems(4);
    auto& carm = world.addItem<Mesh>({ "carm.stl" });
    auto& table = world.addItem<Mesh>({ "table.stl" });
    auto& ctdi = world.addItem<CTDIPhantom>({});
    ctdi.translate({ 16, 0, 9 });

    table.translate({ 0, 0, -17 });
    ctdi.translate({ 0, 0, -17 });
    auto ctdi_aabb = ctdi.AABB();

    world.build();

    // adding beam
    using Beam = dxmc::DXBeam<double>;
    const std::array<double, 3> source_pos = { 16, 0, -70 };
    Beam beam(source_pos);
    // beam.setCollimationAnglesDeg(7, 7);
    beam.setBeamSize(20, 20, 100);

    Viz viz(world);
    auto buffer = viz.generateBuffer(1024, 1024);
    viz.addLineProp(beam, 150, 0.1);

    viz.setDistance(500);
    viz.setAzimuthalAngleDeg(90);
    std::vector<double> angles;
    for (std::size_t i = 0; i < 5; ++i)
        angles.push_back(i * 30);

    for (auto a : angles) {
        viz.setPolarAngleDeg(a);
        viz.suggestFOV();
        viz.generate(world, buffer);
        std::string name = "test" + std::to_string(int(a)) + ".png";
        viz.savePNG(name, buffer);
    }

    return 0;
}