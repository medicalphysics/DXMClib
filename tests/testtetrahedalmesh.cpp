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

    bool success = true;
    for (auto& tet : t)
        success = success && tet.validVerticeOrientation();

    if (!success)
        t.clear();
    return t;
}

template <typename T, std::size_t N = 5, int L = 2, bool Fluence = true>
dxmc::TetrahedalMesh<T, N, L, Fluence> simpletetrahedron()
{
    auto tets = tetCube<T>();

    std::vector<dxmc::Material<T, N>> mats;
    mats.push_back(dxmc::Material<T, N>::byNistName("Water, Liquid").value());
    std::vector<T> dens(tets.size(), 1);

    std::vector<std::string> names(1);
    dxmc::TetrahedalMesh<T, N, L, Fluence> mesh(tets, dens, mats, names, 1);

    return mesh;
}

template <typename T, std::size_t N = 5, int L = 2, bool Fluence = true>
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

    std::vector<dxmc::Material<T, N>> mats;
    mats.push_back(dxmc::Material<T, N>::byNistName("Water, Liquid").value());
    std::vector<T> dens(1, 1);

    std::vector<std::string> names(1);
    dxmc::TetrahedalMesh<T, N, L, Fluence> mesh(tets, dens, mats, names, 1);

    return mesh;
}

bool testintersection()
{
    auto mesh = simpletetrahedron<double, 5, 1, false>();
    dxmc::Particle<double> p { .pos = { 0, 0, -100 }, .dir = { 0, 0, 1 } };
    auto res = mesh.intersect(p);
    bool success = res.valid() && std::abs(res.intersection - 90) < 0.001;
    return success;
}

bool testTransport()
{

    dxmc::World<double, dxmc::TetrahedalMesh<double, 5, 1, true>> w;
    w.reserveNumberOfItems(1);
    auto& mesh = w.addItem(simpletetrahedron<double, 5, 1, true>());
    w.build();

    dxmc::PencilBeam<double> beam({ 0, 0.0001, -100 }, { 0, 0, 1 }, 60);
    beam.setNumberOfExposures(32);
    beam.setNumberOfParticlesPerExposure(10000);
    dxmc::Transport transport;
    transport.setNumberOfThreads(1);
    transport(w, beam);

    for (const auto& tet : mesh.tetrahedrons()) {
        std::cout << tet.doseScored().dose() << std::endl;
    }

    return false;
}

int main()
{

    bool success = true;
    success = success && testintersection();
    success = success && testTransport();
    if (success)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}