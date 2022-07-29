

#include "dxmc/world/kdtree.hpp"
#include "dxmc/world/triangulatedmesh.hpp"

#include <fstream>
#include <iostream>

template <typename T, typename U>
void testGeom(const dxmc::KDTree<T, U>& kdtree)
{
    const std::size_t Nx = 512;
    const std::size_t Ny = 512;
    const std::size_t axis = 0;

    std::vector<T> buffer(Nx * Ny, 0);
    std::vector<std::size_t> idx(buffer.size());
    std::iota(idx.begin(), idx.end(), 0);

    const auto aabb = kdtree.AABB();

    const T dx = (aabb[1 + 3] - aabb[1]) / Nx;
    const T dy = (aabb[2 + 3] - aabb[2]) / Ny;

    std::transform(std::execution::par, idx.cbegin(), idx.cend(), buffer.begin(), [&](const auto i) {
        const std::size_t y = i / Ny;
        const std::size_t x = i - y * Ny;

        const std::array<T, 3> plane { 0, aabb[1] + dx * x, aabb[2] + dy * y };

        dxmc::Particle<T> p;
        p.pos = { -100, 0, 0 };
        for (std::size_t j = 0; j < 3; ++j) {
            p.dir[j] = plane[j] - p.pos[j];
        }
        dxmc::vectormath::normalize(p.dir);

        const auto intersection = kdtree.intersect(p, aabb);

        return intersection ? *intersection : 0;
    });

    std::ofstream file;
    file.open("intersect.bin", std::ios::out | std::ios::binary);
    file.write((char*)buffer.data(), buffer.size() * sizeof(T));
    file.close();
}

int main(int argc, char* argv[])
{
    dxmc::STLReader<double> reader;

    auto mesh = reader("duck.stl");
    std::cout << "Message: " << reader.message() << std::endl;

    dxmc::Particle<double> p;
    p.pos = { -100, 0, 0 };
    p.dir = { 1, 0, 0 };
    dxmc::vectormath::normalize(p.dir);

    const auto aabb_p = mesh.AABB();
    std::array<double, 3> dist;
    for (std::size_t i = 0; i < 3; ++i)
        dist[i] = -(aabb_p[i] + aabb_p[i + 3]) * 0.5;
    mesh.translate(dist);
    const auto aabb = mesh.AABB();

    auto triangles = mesh.getTriangles();
    dxmc::KDTree<double, dxmc::Triangle<double>> kdtree(triangles);
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
    testGeom(kdtree);
    return EXIT_SUCCESS;
}
