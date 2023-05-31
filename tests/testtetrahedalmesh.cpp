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

#include "dxmc/world/visualization/visualizeworld.hpp"
#include "dxmc/world/world.hpp"
#include "dxmc/world/worlditems/tetrahedalmesh.hpp"
#include "dxmc/world/worlditems/tetrahedalmesh/tetrahedalmeshreader.hpp"

#include <iostream>
#include <string>
#include <vector>

template <typename T>
void writeImage(const std::vector<T>& buffer, const std::string& name)
{
    std::ofstream file;
    file.open(name, std::ios::out | std::ios::binary);
    file.write((char*)buffer.data(), buffer.size() * sizeof(T));
    file.close();
}

template <typename T>
bool testReader()
{
    dxmc::TetrahedalmeshReader<T> reader;
    reader.readICRP145Phantom("MRCP_AM.node", "MRCP_AM.ele", "MRCP_AM_media.dat");
    // mesh.readICRP145Phantom("MRCP_AM.node", "MRCP_AM.ele");
    return false;
}

template <typename T>
std::tuple<std::vector<std::array<std::size_t, 4>>, std::vector<std::array<T, 3>>> tetrahedron()
{
    std::vector<std::array<T, 3>> vertices;
    vertices.push_back({ -5, 0, -5 }); // 0
    vertices.push_back({ 0, -5, -5 }); // 1
    vertices.push_back({ 5, 5, -5 }); // 2
    vertices.push_back({ 0, 0, 5 }); // 3

    std::vector<std::array<std::size_t, 4>> nodes;
    nodes.push_back({ 0, 1, 2, 3 });

    return std::make_pair(nodes, vertices);
}

template <typename T>
bool testMeshVisualization()
{

    using Mesh = dxmc::TetrahedalMesh<T, 5, 2>;
    using World = dxmc::World2<T, Mesh>;

    World world;
    auto& mesh = world.addItem<Mesh>({});

    // const auto [nodes, vertices] = tetrahedron<T>();
    // mesh.setData(nodes, vertices);

    world.build(T { 0 });

    dxmc::VisualizeWorld<T> viz(world);
    viz.setDistance(500);
    viz.setPolarAngle(std::numbers::pi_v<T> / 3);
    viz.setAzimuthalAngle(std::numbers::pi_v<T> / 2);
    int height = 512;
    int width = 512;
    std::vector<T> buffer(height * width * 4, T { 1 });

    // viz.setCameraPosition({ 0, -800, -800 });
    viz.suggestFOV();
    viz.generate(world, buffer, width, height);

    writeImage(buffer, "color.bin");
    return false;
}
int main()
{
    std::cout << "Testing ray intersection on tetrahedal mesh\n";

    bool success = true;
    //success = success && testMeshVisualization<double>();
    success = success && testReader<double>();
    success = success && testReader<float>();

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}