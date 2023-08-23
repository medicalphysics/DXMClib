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
#include "dxmc/world/worlditems/ctdiphantom.hpp"
#include "dxmc/world/worlditems/triangulatedmesh.hpp"

#include <vector>

int main()
{
    using CTDIPhantom = dxmc::CTDIPhantom<double>;
    using Mesh = dxmc::TriangulatedMesh<double>;
    using World = dxmc::World<double, CTDIPhantom, Mesh>;
    using Viz = dxmc::VisualizeWorld<double>;

    World world;
    world.reserveNumberOfItems(3);
    auto& carm = world.addItem<Mesh>({ "carm.stl" });
    auto& table = world.addItem<Mesh>({ "table.stl" });
    auto& ctdi = world.addItem<CTDIPhantom>({});

   // carm.scale(100);

    world.build();

    std::size_t N = 512;
    double d = 501;
    std::vector<std::uint8_t> image_buffer(N * N * 4, 0);
    Viz viz(world);
    viz.suggestFOV();
    viz.setCameraPosition({ d, 0, 0 });
    viz.generate(world, image_buffer, N, N);
    viz.savePNG("testx.png", image_buffer, N, N);

    viz.setCameraPosition({ 0, d, 0 });
    viz.suggestFOV();
    viz.generate(world, image_buffer, N, N);
    viz.savePNG("testy.png", image_buffer, N, N);

    viz.setCameraPosition({ 0, 0, d });
    viz.suggestFOV();
    viz.generate(world, image_buffer, N, N);
    viz.savePNG("testz.png", image_buffer, N, N);

    viz.setCameraPosition({ 0, d/2, d });
    viz.suggestFOV();
    viz.generate(world, image_buffer, N, N);
    viz.savePNG("test.png", image_buffer, N, N);
    return 0;
}