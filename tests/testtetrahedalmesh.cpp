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
#include "dxmc/transportprogress.hpp"
#include "dxmc/world/visualization/visualizeworld.hpp"
#include "dxmc/world/world.hpp"
#include "dxmc/world/worlditems/aavoxelgrid.hpp"
#include "dxmc/world/worlditems/tetrahedalmesh.hpp"
#include "dxmc/world/worlditems/tetrahedalmesh/tetrahedalmeshreader.hpp"
#include "dxmc/world/worlditems/worldbox.hpp"

#include <iostream>
#include <string>
#include <vector>

std::vector<dxmc::Tetrahedron> tetCube()
{
    std::vector<std::array<double, 3>> v(8);
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
            n *= 100;

    std::vector<dxmc::Tetrahedron> t(6);

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

template <std::size_t N = 5, int L = 2, bool Fluence = true>
dxmc::TetrahedalMesh<N, L, Fluence> simpletetrahedron()
{
    auto tets = tetCube();

    std::vector<dxmc::Material<N>> mats;
    mats.push_back(dxmc::Material<N>::byNistName("Water, Liquid").value());
    std::vector<double> dens(tets.size(), 1);

    std::vector<std::string> names(1);
    dxmc::TetrahedalMesh<N, L, Fluence> mesh(tets, dens, mats, names, 1);

    return mesh;
}

template <std::size_t N = 5, int L = 2>
dxmc::WorldBox<N, L> simplebox()
{
    dxmc::WorldBox<N, L> box(100);
    auto water = dxmc::Material<N>::byNistName("Water, Liquid").value();
    box.setMaterial(water, 1);
    return box;
}

template <std::size_t N = 5, int L = 2>
dxmc::AAVoxelGrid<N, L, 255> simplegrid()
{

    dxmc::AAVoxelGrid<N, L, 255> grid;
    std::vector<std::uint8_t> m(1000, 0);
    std::vector<double> d(1000, 1);
    std::vector<dxmc::Material<N>> mats;
    mats.push_back(dxmc::Material<N>::byNistName("Water, Liquid").value());
    grid.setData({ 10, 10, 10 }, d, m, mats);
    grid.setSpacing({ 20, 20, 20 });
    return grid;
}

template <std::size_t N = 5, int L = 2, bool Fluence = true>
dxmc::TetrahedalMesh<N, L> simpletetrahedron2()
{
    std::vector<dxmc::Tetrahedron> tets;

    std::vector<std::array<double, 3>> p;
    p.push_back({ -1, 0, -1 });
    p.push_back({ 1, 0, -1 });
    p.push_back({ 1, 2, -1 });
    p.push_back({ 1, 0, 1 });

    tets.push_back({ p[0], p[1], p[2], p[3] });

    bool valid = true;
    for (const auto& t : tets)
        valid = valid && t.validVerticeOrientation();

    std::vector<dxmc::Material<N>> mats;
    mats.push_back(dxmc::Material<N>::byNistName("Water, Liquid").value());
    std::vector<double> dens(1, 1);

    std::vector<std::string> names(1);
    dxmc::TetrahedalMesh<N, L, Fluence> mesh(tets, dens, mats, names, 1);

    return mesh;
}

bool testIntersectionMesh()
{
    auto mesh = simpletetrahedron<5, 1, false>();
    dxmc::Particle p { .pos = { 0, 0, -1000 }, .dir = { 0, 0, 1 } };
    auto res = mesh.intersect(p);
    bool success = res.valid() && std::abs(res.intersection - 90) < 0.001;
    return success;
}

bool testTransport()
{

    using M1 = dxmc::TetrahedalMesh<5, 1, true>;
    using M2 = dxmc::TetrahedalMesh<5, 1, false>;
    using B = dxmc::WorldBox<5, 1>;
    using G = dxmc::AAVoxelGrid<5, 1, 255>;

    constexpr std::size_t N_HIST = 1000000;
    constexpr std::size_t N_EXP = 16;

    dxmc::PencilBeam<> beam({ 50, 50, -1000 }, { 0, 0, 1 }, 60);
    beam.setNumberOfExposures(N_EXP);
    beam.setNumberOfParticlesPerExposure(N_HIST);

    dxmc::TransportProgress progress;

    {
        dxmc::World<M1> w;
        w.reserveNumberOfItems(1);
        auto& mesh = w.addItem(simpletetrahedron<5, 1, true>());
        w.build();

        double dose = 0;
        dxmc::Transport transport;

        transport(w, beam, &progress);
        for (const auto& t : mesh.tetrahedrons())
            dose += t.energyScored().energyImparted();
        std::cout << dose << " " << progress.humanTotalTime() << std::endl;
    }

    {
        dxmc::World<M2> w;
        w.reserveNumberOfItems(1);
        auto& mesh = w.addItem(simpletetrahedron<5, 1, false>());
        w.build();

        double dose = 0;
        dxmc::Transport transport;
        // transport.setNumberOfThreads(1);
        transport(w, beam, &progress);
        for (const auto& t : mesh.tetrahedrons())
            dose += t.energyScored().energyImparted();
        std::cout << dose << " " << progress.humanTotalTime() << std::endl;
    }

    {
        dxmc::World<G> w;
        w.reserveNumberOfItems(1);
        auto& mesh = w.addItem(simplegrid<5, 1>());
        w.build();

        dxmc::Transport transport;
        transport(w, beam, &progress);
        auto s = mesh.size();
        double dose = 0;
        for (std::size_t i = 0; i < s; ++i)
            dose += mesh.energyScored(i).energyImparted();
        std::cout << dose << " " << progress.humanTotalTime() << std::endl;
    }

    {
        dxmc::World<B> w;
        w.reserveNumberOfItems(1);
        auto& box = w.addItem(simplebox<5, 1>());
        w.build();

        dxmc::Transport transport;
        transport(w, beam, &progress);
        auto dose = box.energyScored().energyImparted();
        std::cout << dose << " " << progress.humanTotalTime() << std::endl;
    }

    return false;
}

bool testIntersection()
{

    dxmc::Particle p { .pos = { 0, 0, -100 }, .dir = { 0, 0, 1 } };
    const auto tets = tetCube();

    bool success = true;
    for (const auto& tet : tets) {
        auto i1 = tet.intersect(p);
        auto p2 = p;
        if (i1.valid()) {
            p2.border_translate(i1.intersection);
            bool is_inside = tet.pointInside(p2.pos);
            auto i2 = tet.intersect(p2);
            success = success && is_inside == i2.rayOriginIsInsideItem;
        }
    }
    return success;
}

int main()
{

    bool success = true;
    success = success && testIntersection();
    success = success && testTransport();
    if (success)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}