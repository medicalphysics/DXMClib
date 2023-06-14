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

#include "dxmc/beams/pencilbeam.hpp"
#include "dxmc/transport.hpp"
#include "dxmc/world/visualization/visualizeworld.hpp"
#include "dxmc/world/world.hpp"
#include "dxmc/world/worlditems/tetrahedalmesh.hpp"
#include "dxmc/world/worlditems/tetrahedalmesh/tetrahedalmeshreader.hpp"
#include "dxmc/world/worlditems/worldbox.hpp"

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
    reader.readICRP145Phantom("MRCP_AM.node", "MRCP_AM.ele", "MRCP_AM_media.dat", "icrp145organs.csv");
    // mesh.readICRP145Phantom("MRCP_AM.node", "MRCP_AM.ele");
    return false;
}

template <typename T>
std::vector<dxmc::Tetrahedron<T>> tetCube()
{
    std::vector<std::array<T, 3>> v(8);
    v[0] = { -1, 1, 1 };
    v[1] = { 1, -1, 1 };
    v[2] = { 1, 1, 1 };
    v[3] = { -1, -1, -1 };
    v[4] = { 1, -1, -1 };
    v[5] = { -1, -1, 1 };
    v[6] = { -1, 1, -1 };
    v[7] = { 1, 1, -1 };

    for (auto& i : v)
        for (auto& n : i)
            n *= 10;

    std::vector<dxmc::Tetrahedron<T>> t(6);
    t[0] = { v[1], v[7], v[0], v[2], 0, 0 }; //*
    t[1] = { v[7], v[3], v[0], v[6], 0, 0 }; //*
    t[2] = { v[1], v[3], v[0], v[4], 0, 0 }; //*
    t[3] = { v[1], v[7], v[4], v[0], 0, 0 }; //*
    t[4] = { v[7], v[3], v[4], v[0], 0, 0 }; //*
    t[5] = { v[1], v[3], v[5], v[0], 0, 0 }; //*

    for (auto& tet : t)
        tet.validVerticeOrientation();

    // std::vector<dxmc::Tetrahedron<T>> t(1);
    // t[0] = { v[2], v[6], v[3], v[0], 0, 0 };

    return t;
}

template <typename T, std::size_t N = 5, int L = 2>
dxmc::TetrahedalMesh<T, N, L> simpletetrahedron()
{
    auto tets = tetCube<T>();

    std::vector<dxmc::Material2<T, N>> mats;
    mats.push_back(dxmc::Material2<T, N>::byNistName("Water, Liquid").value());
    std::vector<T> dens(tets.size(), 1);

    std::vector<std::string> names(1);
    dxmc::TetrahedalMesh<T, N, L> mesh(std::move(tets), dens, mats, names);

    return mesh;
}

template <typename T, std::size_t N = 5, int L = 2>
dxmc::TetrahedalMesh<T, N, L> simpletetrahedron2()
{
    std::vector<dxmc::Tetrahedron<T>> tets;

    std::vector<std::array<T, 3>> p;
    p.push_back({ -1, 0, -1 });
    p.push_back({ 1, 0, -1 });
    p.push_back({ 1, 2, -1 });
    p.push_back({ 1, 0, 1 });

    tets.push_back({ p[0], p[1], p[2], p[3] });

    bool valid = true;
    for (const auto& t : tets)
        valid = valid && t.validVerticeOrientation();

    std::vector<dxmc::Material2<T, N>> mats;
    mats.push_back(dxmc::Material2<T, N>::byNistName("Water, Liquid").value());
    std::vector<T> dens(1, 1);

    std::vector<std::string> names(1);
    dxmc::TetrahedalMesh<T, N, L> mesh(std::move(tets), dens, mats, names);

    return mesh;
}

template <typename T, std::size_t N = 5, int L = 2>
void testMeshCubeVisualization()
{

    using Mesh = dxmc::TetrahedalMesh<T, N, L>;
    using World = dxmc::World2<T, Mesh>;

    World world;

    // dxmc::TetrahedalmeshReader<T> reader;
    // auto& mesh = world.template addItem<Mesh>(reader.readICRP145Phantom("MRCP_AM.node", "MRCP_AM.ele", "MRCP_AM_media.dat", "icrp145organs.csv"));

    auto& mesh = world.template addItem<Mesh>(simpletetrahedron<T, N, L>());

    // const auto [nodes, vertices] = tetrahedron<T>();
    // mesh.setData(nodes, vertices);

    world.build(T { 0 });

    dxmc::VisualizeWorld<T> viz(world);

    int height = 512;
    int width = 512;
    std::vector<T> buffer(height * width * 4, T { 1 });

    for (std::size_t i = 0; i < 12; ++i) {
        viz.setDistance(500);
        viz.setPolarAngle(std::numbers::pi_v<T> / 3.0);
        viz.setAzimuthalAngle(std::numbers::pi_v<T> * i / 6.0);
        // viz.setCameraPosition({ -60, -30, -10 });
        viz.suggestFOV();
        viz.generate(world, buffer, width, height);
        std::string name = "color_" + std::to_string(i) + ".bin";
        writeImage(buffer, name);
    }
}

template <dxmc::Floating T, std::size_t N = 5, int L = 2, bool BOX = false>
T testDoseScoring(std::array<T, 3> pos = { 0, -100, 0 }, std::array<T, 3> dir = { 0, 1, 0 })
{
    using Mesh = dxmc::TetrahedalMesh<T, N, L>;
    using Box = dxmc::WorldBox<T, N, L>;
    using World = dxmc::World2<T, Mesh, Box>;
    using Material = dxmc::Material2<T, N>;
    World world;

    if constexpr (BOX) {
        auto& box = world.template addItem<Box>({ 10 });
        auto water = Material::byNistName("Water, Liquid").value();
        const T density = 1;
        box.setMaterial(water, density);
    } else {
        world.template addItem<Mesh>(simpletetrahedron<T, N, L>());
    }
    world.build();

    auto c = world.center();
    auto d = dxmc::vectormath::subtract({ 0, 0, 0 }, c);
    world.translate(d);
    auto c1 = world.center();

    dxmc::PencilBeam<T> beam(pos, dir);
    beam.setNumberOfParticlesPerExposure(1E6);
    beam.setNumberOfExposures(4);

    dxmc::Transport transport;
    // transport.setNumberOfThreads(1);
    transport(world, beam);
    if constexpr (BOX) {
        auto& boxes = world.template getItems<Box>();
        return boxes[0].dose(0).energyImparted();
    } else {
        auto& meshes = world.template getItems<Mesh>();
        auto& mesh = meshes[0];
        T dose = 0;
        for (std::size_t i = 0; i < mesh.numberOfCollections(); ++i)
            dose += mesh.dose(i).energyImparted();
        return dose;
    }
}

int main()
{
    std::cout << "Testing ray intersection on tetrahedal mesh\n";

    testMeshCubeVisualization<double>();

    bool success = true;
    std::cout << "Test tetrahedalmesh dose scoring\n";

    std::vector<std::pair<std::array<double, 3>, std::array<double, 3>>> bp;
    bp.push_back(std::make_pair(std::array<double, 3> { -100, 0, 0 }, std::array<double, 3> { 1, 0, 0 }));
    bp.push_back(std::make_pair(std::array<double, 3> { 0, -100, 0 }, std::array<double, 3> { 0, 1, 0 }));
    bp.push_back(std::make_pair(std::array<double, 3> { 0, 0, -100 }, std::array<double, 3> { 0, 0, 1 }));
    bp.push_back(std::make_pair(std::array<double, 3> { 100, 0, 0 }, std::array<double, 3> { -1, 0, 0 }));
    bp.push_back(std::make_pair(std::array<double, 3> { 0, 100, 0 }, std::array<double, 3> { 0, -1, 0 }));
    bp.push_back(std::make_pair(std::array<double, 3> { 0, 0, 100 }, std::array<double, 3> { 0, 0, -1 }));

    for (const auto& [pos, dir] : bp) {
        const auto tm = testDoseScoring<double, 5, 1, false>({ -100, 0, 0 }, { 1, 0, 0 });
        const auto tc = testDoseScoring<double, 5, 1, true>({ -100, 0, 0 }, { 1, 0, 0 });
        std::cout << "Dose to mesh: " << tm << ", dose cube: " << tc << ", difference[%]: " << (tm / tc - 1) * 100 << "\n";
    }

    if (success)
        std::cout << "SUCCESS ";
    else
        std::cout << "FAILURE ";

    // success = success && testReader<double>();
    // success = success && testReader<float>();

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}