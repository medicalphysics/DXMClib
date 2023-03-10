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
#include "dxmc/world/worlditems/ctdiphantom.hpp"
#include "dxmc/world/worlditems/depthdose.hpp"
#include "dxmc/world/worlditems/worldbox.hpp"
#include "dxmc/world/worlditems/worldsphere.hpp"

#include <chrono>
#include <fstream>
#include <iostream>
#include <thread>

template <typename W, typename T>
auto create_image(W& world, const std::array<T, 3>& camera_pos, bool print = true)
{
    const std::int64_t Nx = 1024;
    const std::int64_t Ny = 1024;

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

template <typename T>
void testTransportWorker(dxmc::World2<T, dxmc::DepthDose<T>>& w, dxmc::Particle<T> p_temp, const std::size_t N)
{
    dxmc::RandomState state;
    for (std::size_t i = 0; i < N; ++i) {
        auto p = p_temp;
        w.transport(p, state);
    }
}

template <typename T>
bool testTransport(bool print = false)
{
    constexpr std::size_t N = 1E6;

    dxmc::World2<T, dxmc::DepthDose<T>> w;

    {
        std::array<T, 3> pos = { 0, 0, 0 };
        dxmc::DepthDose<T> sp(1);
        sp.setNistMaterial("Water, Liquid");
        sp.setMaterialDensity(1);
        w.addItem(std::move(sp));
    }

    w.build();

    const dxmc::Particle<T> p_temp = {
        .pos = { 0, 0, -1000 },
        .dir = { 0, 0, 1 },
        .energy = 60,
        .weight = 1
    };

    const auto tstart = std::chrono::high_resolution_clock::now();

    std::vector<std::jthread> threads;
    for (std::size_t i = 0; i < std::jthread::hardware_concurrency() - 1; ++i) {
        threads.emplace_back(testTransportWorker<T>, std::ref(w), p_temp, N);
    }

    testTransportWorker<T>(w, p_temp, N);
    for (auto& t : threads) {
        t.join();
    }
    const auto tend = std::chrono::high_resolution_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(tend - tstart);
    auto items = w.template getItems<dxmc::DepthDose<T>>();

    std::vector<T> diff;
    for (const auto& it : items) {
        auto v = it.doseVector();
        auto d = it.depthVector();

        auto mat = it.material();
        const auto att_tot = mat.attenuationValues(p_temp.energy).sum();
        auto dens = it.density();

        for (std::size_t i = 0; i < d.size(); ++i) {
            if (print)
                std::cout << d[i] << ", " << v[i].energyImparted() << std::endl;

            const T ana = std::exp(-(d[i] - d[0]) * att_tot * dens);
            const T expr = v[i].energyImparted() / v[0].energyImparted();
            diff.push_back(std::abs(ana - expr));
        }
    }

    auto max_diff = std::max_element(diff.cbegin(), diff.cend());

    auto success = *max_diff < 0.1;
    if (success)
        std::cout << "SUCCESS ";
    else
        std::cout << "FAILURE ";
    std::cout << "Depth dose test with sizeof(T) = " << sizeof(T);
    std::cout << " [time: " << time.count() << "ms] for " << N * (threads.size() + 1) << " particles\n";
    return success;
}

template <typename T, typename WO>
bool testWorldItem()
{

    WO obj;

    dxmc::Particle<T> p;
    p.dir = { 0, 0, 1 };
    p.pos = { 0, 0, -2000 };

    dxmc::vectormath::normalize(p.dir);

    auto t = obj.intersect(p);
    if (!t.valid())
        return false;
    p.border_translate(t.intersection);

    dxmc::RandomState state;
    p.energy = 60;
    p.weight = 1;
    obj.transport(p, state);

    auto t_fin = obj.intersect(p);
    if (t_fin.valid())
        return false;

    return true;
}

template <typename T>
bool testWorldItems()
{
    bool success = true;

    success = success && testWorldItem<T, dxmc::CTDIPhantom<T>>();
    if (success)
        std::cout << "SUCCESS Basic test for CTDIPhantom<T> for sizeof(T) = " << sizeof(T) << std::endl;
    else
        std::cout << "FAILURE Basic test for CTDIPhantom<T> for sizeof(T) = " << sizeof(T) << std::endl;

    success = success && testWorldItem<T, dxmc::WorldBox<T>>();
    if (success)
        std::cout << "SUCCESS Basic test for WorldBox<T> for sizeof(T) = " << sizeof(T) << std::endl;
    else
        std::cout << "FAILURE Basic test for WorldBox<T> for sizeof(T) = " << sizeof(T) << std::endl;

    success = success && testWorldItem<T, dxmc::WorldSphere<T>>();
    if (success)
        std::cout << "SUCCESS Basic test for WorldSphere<T> for sizeof(T) = " << sizeof(T) << std::endl;
    else
        std::cout << "FAILURE Basic test for WorldSphere<T> for sizeof(T) = " << sizeof(T) << std::endl;
    return success;
}

template <typename T>
bool testBorderCrossing()
{

    dxmc::World2<T, dxmc::WorldBox<T>> w;

    dxmc::WorldBox<T> b1, b2;
    b2.translate({ 0, 0, 2 });
    w.addItem(b1);
    w.addItem(b2);
    w.build();

    dxmc::Particle<T> p;
    p.dir = { 0, 0, 1 };
    p.pos = { 0, 0, -2000 };
    dxmc::vectormath::normalize(p.dir);

    dxmc::RandomState state;
    p.energy = 60;
    p.weight = 1;
    w.transport(p, state);

    return true;
}

int main(int argc, char* argv[])
{
    std::cout << "World Tests\n";
    auto success = true;

    success = success && testBorderCrossing<float>();
    success = success && testBorderCrossing<double>();

    // success = success && testTransport<float>();
    // success = success && testTransport<double>();

    success = success && testWorldItems<double>();
    success = success && testWorldItems<float>();

    if (success)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}
