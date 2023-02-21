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

#include "dxmc/dxmcrandom.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/world.hpp"
#include "dxmc/world/worlditems/worldbox.hpp"

#include <chrono>
#include <fstream>
#include <iostream>

template <typename W, typename T>
auto create_image(W& world, const std::array<T, 3>& camera_pos, bool print = true)
{
    const std::int64_t Nx = 512;
    const std::int64_t Ny = 512;

    std::vector<T> buffer(Nx * Ny, 0);
    std::vector<std::size_t> idx(buffer.size());
    std::iota(idx.begin(), idx.end(), 0);

    const auto& aabb = world.AABB();
    const std::array<T, 3> aabb_center = { (aabb[0] + aabb[3]) / 2, (aabb[1] + aabb[4]) / 2, (aabb[2] + aabb[5]) / 2 };
    const std::array<T, 3> aabb_extent = { (aabb[3] - aabb[0]), (aabb[4] - aabb[1]), (aabb[5] - aabb[3]) };

    auto dir = dxmc::vectormath::subtract(aabb_center, camera_pos);
    dxmc::vectormath::normalize(dir);
    auto collapse_dim = dxmc::vectormath::argmax3(dir);

    std::array<T, 3> dx = { 1, 0, 0 };
    std::array<T, 3> dy = { 0, 1, 0 };
    std::array<T, 3> dz = { 0, 0, 1 };
    std::array<T, 3> cm_canditates { dxmc::vectormath::dot(dx, dir), dxmc::vectormath::dot(dy, dir), dxmc::vectormath::dot(dz, dir) };
    const std::array<T, 3> c1 = dxmc::vectormath::argmax3(cm_canditates) == 0 ? dxmc::vectormath::cross(dir, dz) : dxmc::vectormath::cross(dir, dx);
    const std::array<T, 3> c2 = dxmc::vectormath::cross(dir, c1);

    const std::array<T, 3> aabb_p1 { aabb[0], aabb[1], aabb[2] };
    const std::array<T, 3> aabb_p2 { aabb[3], aabb[4], aabb[5] };

    const T FOV = dxmc::vectormath::lenght(dxmc::vectormath::subtract(aabb_p1, aabb_p2));

    const auto t0 = std::chrono::high_resolution_clock::now();
    std::transform(std::execution::par_unseq, idx.cbegin(), idx.cend(), buffer.begin(), [&](const auto i) {
        const std::int64_t y = i / Nx;
        const std::int64_t x = i - y * Nx;
        const T dx = (FOV / Nx) * (Nx / 2 - x);
        const T dy = (FOV / Nx) * (Ny / 2 - y);

        std::array<T, 3> target;
        for (std::size_t j = 0; j < 3; ++j) {
            target[j] = aabb_center[j] + c1[j] * dx + c2[j] * dy;
        }

        dxmc::Particle<T> p;
        p.pos = camera_pos;
        p.dir = dxmc::vectormath::subtract(target, camera_pos);

        dxmc::vectormath::normalize(p.dir);

        const auto intersection = world.intersect(p);
        if (intersection.item)
            return intersection.intersection;
        return T { 0 };
    });
    const auto t1 = std::chrono::high_resolution_clock::now();
    if (print) {
        std::cout << "Intersect time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << std::endl;
        std::ofstream file;
        file.open("world.bin", std::ios::out | std::ios::binary);
        file.write((char*)buffer.data(), buffer.size() * sizeof(T));
        file.close();
    }
    return t1 - t0;
}

template <typename T, typename WO>
bool testWorldItem()
{

    WO obj;

    dxmc::Particle<T> p;
    p.dir = { 0, 0, 1 };
    p.pos = { 0, 0, -2000 };

    auto t = obj.intersectForward(p);
    if (!t)
        return false;
    p.translate(*t);

    dxmc::RandomState state;
    p.energy = 60;
    p.weight = 1;
    obj.transport(p, state);

    auto t_fin = obj.intersectForward(p);
    if (t_fin)
        return false;

    return true;
}

template <typename T>
bool testWorldItems()
{
    bool success = true;

    success = success && testWorldItem<T, dxmc::WorldBox<T>>();

    return success;
}

template <typename T>
bool testCreateWorld()
{

    dxmc::World2<T, dxmc::WorldBox<T>> world;

    for (std::size_t i = 0; i < 2; ++i) {
        dxmc::WorldBox<T> box;
        std::array<T, 3> dir = { 0, 0, 5.0f * i };
        box.translate(dir);
        auto material_set = box.setNistMaterial("Water, Liquid");
        world.addItem(std::move(box));
    }

    world.build();

    dxmc::Particle<T> p;
    p.pos = { 0, 0, -5 };
    p.dir = { 0, 0, 1 };
    p.energy = 60;
    p.weight = 1;

    dxmc::RandomState state;

    world.transport(p, state);

    return true;
}

int main(int argc, char* argv[])
{
    auto success = true;

    success = success && testWorldItems<double>();
    success = success && testWorldItems<float>();

    success = success && testCreateWorld<double>();
    success = success && testCreateWorld<float>();

    if (success)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}
