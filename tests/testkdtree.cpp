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

#include "dxmc/world/kdtree.hpp"
#include "dxmc/world/statickdtree.hpp"
#include "dxmc/world/triangle.hpp"
#include "dxmc/world/triangulatedmesh.hpp"
#include "dxmc/world/ctdiphantom.hpp"

#include <chrono>
#include <fstream>
#include <iostream>

template <typename U>
void test_ptr(U& kdtree, const std::size_t depth = 8)
{
    using T = typename U::Type;
    const std::size_t Nx = 512;
    const std::size_t Ny = 512;
    const std::size_t axis = 0;

    std::vector<T> buffer(Nx * Ny, 0);
    std::vector<std::size_t> idx(buffer.size());
    std::iota(idx.begin(), idx.end(), 0);

    const auto aabb = kdtree.AABB();

    const T dx = (aabb[1 + 3] - aabb[1]) / Nx;
    const T dy = (aabb[2 + 3] - aabb[2]) / Ny;

    const auto t0 = std::chrono::high_resolution_clock::now();

    std::transform(std::execution::par, idx.cbegin(), idx.cend(), buffer.begin(), [&](const auto i) {
        const std::size_t y = i / Ny;
        const std::size_t x = i - y * Ny;

        const std::array<T, 3> plane { 0, aabb[1] + dx * x, aabb[2] + dy * y };

        dxmc::Particle<T> p;
        p.pos = { -1000, 0, 0 };
        for (std::size_t j = 0; j < 3; ++j) {
            p.dir[j] = plane[j] - p.pos[j];
        }
        dxmc::vectormath::normalize(p.dir);

        const auto intersection = kdtree.intersect(p, aabb);

        return intersection ? *intersection : 0;
    });
    const auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "KDTree depth: " << kdtree.depth();
    std::cout << "    Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << std::endl;
    std::ofstream file;
    file.open("intersect.bin", std::ios::out | std::ios::binary);
    file.write((char*)buffer.data(), buffer.size() * sizeof(T));
    file.close();
}

template <typename T>
void testGeom(const dxmc::TriangulatedMesh<T>& mesh, const std::size_t depth = 6)
{
    const std::size_t Nx = 64;
    const std::size_t Ny = 64;
    const std::size_t axis = 0;

    std::vector<T> buffer(Nx * Ny, 0);
    std::vector<std::size_t> idx(buffer.size());
    std::iota(idx.begin(), idx.end(), 0);

    auto triangles = mesh.getTriangles();
    dxmc::KDTree<dxmc::Triangle<T>> kdtree(triangles, depth);
    const auto aabb = kdtree.AABB();

    const T dx = (aabb[1 + 3] - aabb[1]) / Nx;
    const T dy = (aabb[2 + 3] - aabb[2]) / Ny;

    const auto t0 = std::chrono::high_resolution_clock::now();

    std::transform(std::execution::par, idx.cbegin(), idx.cend(), buffer.begin(), [&](const auto i) {
        const std::size_t y = i / Ny;
        const std::size_t x = i - y * Ny;

        const std::array<T, 3> plane { 0, aabb[1] + dx * x, aabb[2] + dy * y };

        dxmc::Particle<T> p;
        p.pos = { -1000, 0, 0 };
        for (std::size_t j = 0; j < 3; ++j) {
            p.dir[j] = plane[j] - p.pos[j];
        }
        dxmc::vectormath::normalize(p.dir);

        const auto intersection = kdtree.intersect(p, aabb);
        const auto intersection2 = mesh.intersect(p);

        return intersection ? *intersection : 0;
    });
    const auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "KDTree depth: " << kdtree.depth();
    std::cout << "    Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << std::endl;
    std::ofstream file;

    file.open("intersect.bin", std::ios::out | std::ios::binary);
    file.write((char*)buffer.data(), buffer.size() * sizeof(T));
    file.close();
}


int main(int argc, char* argv[])
{

    dxmc::TriangulatedMesh<double> mesh("bunny_low.stl");

    dxmc::Particle<double> p;
    p.pos = { 0, 0, 0 };
    p.dir = { 1, 0, 0 };
    dxmc::vectormath::normalize(p.dir);

    const auto aabb_p = mesh.AABB();
    std::array<double, 3> dist;
    for (std::size_t i = 0; i < 3; ++i)
        dist[i] = -(aabb_p[i] + aabb_p[i + 3]) * 0.5;
    mesh.translate(dist);
    const auto aabb = mesh.AABB();

    auto triangles = mesh.getTriangles();

    std::vector<dxmc::Triangle<double>*> triangles_ptr;
    for (auto& t : triangles)
        triangles_ptr.push_back(&t);

    auto kdtree = dxmc::KDTree(triangles);
    auto kdtree_ptr = dxmc::KDTree(triangles_ptr);
    const auto& kdtree_m = mesh.kdtree();

    // test_ptr(kdtree);
    // test_ptr(kdtree_ptr);

    std::cout << "KDTree depth : " << kdtree.depth() << std::endl;

    const auto aabb_k = kdtree.AABB();

    auto res = kdtree.intersect(p, aabb);
    if (res)
        std::cout << "Intersect KDTree: " << *res << std::endl;

    // brute
    double bh_min = std::numeric_limits<double>::max();
    double bh_max = std::numeric_limits<double>::lowest();
    for (const auto& t : triangles) {
        auto h = t.intersect<0>(p);
        if (h) {
            bh_min = std::min(*h, bh_min);
            bh_max = std::max(*h, bh_max);
        }
    }
    std::cout << "Brute force hits: " << bh_min << ", " << bh_max << std::endl;

    testGeom(mesh, 5);

    return EXIT_SUCCESS;
}
