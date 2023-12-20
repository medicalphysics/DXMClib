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
#include "dxmc/world/worlditems/enclosedroom.hpp"
#include "dxmc/world/worlditems/fluencescore.hpp"
#include "dxmc/world/worlditems/tetrahedalmesh.hpp"
#include "dxmc/world/worlditems/triangulatedmesh.hpp"
#include "dxmc/world/worlditems/triangulatedopensurface.hpp"
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

    if constexpr (std::is_same_v<U, dxmc::TriangulatedMesh<T>> || std::is_same_v<U, dxmc::TriangulatedOpenSurface<T>>) {
        std::vector<std::array<T, 3>> p;
        constexpr T d = 30;
        p.push_back({ 1, 1, 0 }); // 0
        p.push_back({ 1, -1, 0 }); // 1
        p.push_back({ -1, -1, 0 }); // 2
        p.push_back({ -1, 1, 0 }); // 3
        p.push_back({ 0, 0, 1 });
        for (auto& i : p)
            for (auto& j : i)
                j *= d;

        std::vector<dxmc::Triangle<T>> t;
        t.push_back({ p[0], p[1], p[4] });
        t.push_back({ p[1], p[2], p[4] });
        t.push_back({ p[2], p[3], p[4] });
        t.push_back({ p[3], p[0], p[4] });

        if constexpr (std::is_same_v<U, dxmc::TriangulatedMesh<T>>) {
            // underside
            t.push_back({ p[0], p[3], p[2] });
            t.push_back({ p[2], p[1], p[0] });
        }
        item.setData(t);

    } else if constexpr (std::is_same_v<U, dxmc::TetrahedalMesh<T>>) {
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

        std::vector<T> dens(1, 1);
        std::vector<dxmc::Material<T>> mats;
        mats.push_back(dxmc::Material<T>::byChemicalFormula("H2O").value());

        item.setData(t, dens, mats);
    }

    item.translate({ 1, 1, 1 });
    item.translate({ -1, -1, -1 });
    auto center = item.center();
    auto aabb = item.AABB();
    dxmc::Particle<T> p { .pos = { 0, 0, 0 }, .dir = { 0, 0, 1 } };
    auto intersection = item.intersect(p);
    auto intersectionViz = item.intersectVisualization(p);

    world.build();
    dxmc::PencilBeam<T> beam;
    beam.setNumberOfExposures(1);
    beam.setNumberOfParticlesPerExposure(8);
    dxmc::Transport transport;
    transport.setNumberOfThreads(1);
    transport(world, beam);

    auto energy = item.energyScored();
    item.addEnergyScoredToDoseScore();
    item.clearEnergyScored();
    auto dose = item.doseScored(0);
    item.clearDoseScored();
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
    success = success && testItem<T, dxmc::TriangulatedOpenSurface<T>>();
    success = success && testItem<T, dxmc::WorldBox<T>>();
    success = success && testItem<T, dxmc::WorldBoxGrid<T>>();
    success = success && testItem<T, dxmc::WorldCylinder<T>>();
    success = success && testItem<T, dxmc::WorldSphere<T>>();
    success = success && testItem<T, dxmc::EnclosedRoom<T>>();

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
