

#include "dxmc/world/triangulatedmesh.hpp"

#include <chrono>
#include <fstream>
#include <iostream>

template <dxmc::KDTreeType U>
void create_image(U& object)
{
    using T = typename U::Type;
    const std::size_t Nx = 64;
    const std::size_t Ny = 64;
    const std::size_t axis = 0;

    std::vector<T> buffer(Nx * Ny, 0);
    std::vector<std::size_t> idx(buffer.size());
    std::iota(idx.begin(), idx.end(), 0);

    const auto aabb = object.AABB();

    const T dx = (aabb[1 + 3] - aabb[1]) / Nx;
    const T dy = (aabb[2 + 3] - aabb[2]) / Ny;

    const auto t0 = std::chrono::high_resolution_clock::now();

    std::transform(std::execution::par, idx.cbegin(), idx.cend(), buffer.begin(), [&](const auto i) {
        const std::size_t y = i / Ny;
        const std::size_t x = i - y * Ny;

        const std::array<T, 3> plane { 0, aabb[1] + dx * x, aabb[2] + dy * y };

        dxmc::Particle<T> p;
        p.pos = { -200, 0, 0 };
        for (std::size_t j = 0; j < 3; ++j) {
            p.dir[j] = plane[j] - p.pos[j];
        }
        dxmc::vectormath::normalize(p.dir);

        const auto intersection = object.intersect(p);

        return intersection ? *intersection : 0;
    });
    const auto t1 = std::chrono::high_resolution_clock::now();
    std::cout << "Intersect time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count() << std::endl;
    std::ofstream file;
    file.open("intersect.bin", std::ios::out | std::ios::binary);
    file.write((char*)buffer.data(), buffer.size() * sizeof(T));
    file.close();
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

int main(int argc, char* argv[])
{
    dxmc::TriangulatedMesh<double> mesh("bunny_low.stl", 5);
    //dxmc::TriangulatedMesh<double> mesh("duck.stl", 5);

    const auto& aabb_pre = mesh.AABB();

    auto triangles = mesh.getTriangles();
    std::cout << "Number of triangles: " << triangles.size() << std::endl;

    std::cout << "AABB before translate: ";
    for (auto v : aabb_pre) {
        std::cout << v << " ";
    }
    std::cout << std::endl;
    auto t = neg_center(aabb_pre);
    mesh.translate(t);
    const auto& aabb = mesh.AABB();
    std::cout << "AABB after translate: ";
    for (auto v : aabb) {
        std::cout << v << " ";
    }
    std::cout << std::endl;

    const auto& kdtree = mesh.kdtree();
    std::cout << "Tree depth: " << kdtree.depth() << std::endl;

    create_image(mesh);

    return EXIT_SUCCESS;
}
