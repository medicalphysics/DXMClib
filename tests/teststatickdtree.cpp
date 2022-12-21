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

#include "dxmc/world/ctdiphantom.hpp"
#include "dxmc/world/statickdtree.hpp"
#include "dxmc/world/triangle.hpp"

#include <chrono>
#include <fstream>
#include <iostream>

template <typename T, typename U>
auto create_image(const U& tree, const std::array<T, 3>& camera_pos, bool print = true)
{
    const std::int64_t Nx = 512;
    const std::int64_t Ny = 512;

    std::vector<T> buffer(Nx * Ny, 0);
    std::vector<std::size_t> idx(buffer.size());
    std::iota(idx.begin(), idx.end(), 0);

    const auto aabb = tree.AABB();
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

        const auto intersection = tree.intersect(p);

        return intersection ? *intersection : 0;
    });
    const auto t1 = std::chrono::high_resolution_clock::now();
    if (print) {
        std::cout << "Intersect time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << std::endl;
        std::ofstream file;
        file.open("intersect.bin", std::ios::out | std::ios::binary);
        file.write((char*)buffer.data(), buffer.size() * sizeof(T));
        file.close();
    }
    return t1 - t0;
}
template <dxmc::Floating T>
std::array<T, 3> neg_center(const std::array<T, 6>& aabb)
{
    std::array<T, 3> c;
    for (int i = 0; i < 3; ++i) {
        c[i] = -(aabb[i] + aabb[i + 3]) / 2;
    }
    return c;
}

template <typename T>
bool teststatickdtree()
{

    // dxmc::KDTreeNode<T, dxmc::Triangle<T>> node;

    bool success = true;

    dxmc::KDTree<T, dxmc::Triangle<T>, dxmc::CTDIPhantom<T>> tree;

    std::array<T, 9> vertices = { 0, 0, 0, 10, 0, 0, 0, 0, 10 };
    dxmc::Triangle<T> tri(vertices.data());

    for (int i = -1; i < 2; ++i) {
        for (int j = -1; j < 2; ++j) {
            for (int k = -1; k < 2; ++k) {
                constexpr T step = 80;
                std::array<T, 3> pos { step * i, step * j, step * k };
                dxmc::CTDIPhantom<T> ph { 16, pos, 15 };
                tree.insert(ph);
            }
        }
    }

    // dxmc::CTDIPhantom<T> ph { 16, { 0, 0, 0 }, 32 };
    // tree.insert(ph);

    // translate all objects
    std::array<T, 3> trans { 1, 1, 1 };
    // tree.translate(trans);

    // calculate aabb
    auto aabb = tree.AABB();

    // calculate center
    auto center = tree.center();
    tree.build();
    auto depth = tree.depth();

    // std::array<T, 3> camera { 900, 1000, 800 };
    std::array<T, 3> camera { 1000, 1000, 1000 };
    auto time = create_image(tree, camera);
    std::cout << " Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(time).count() << std::endl;
    return true;
}

int main()
{

    auto success = teststatickdtree<float>();
    // success = success && teststatickdtree<double>();
    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}