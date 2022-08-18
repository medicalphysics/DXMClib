

#include "dxmc/world/triangulatedmesh.hpp"

#include <chrono>
#include <fstream>
#include <iostream>

template <typename T>
auto create_image(dxmc::TriangulatedMesh<T>& object, const std::array<T, 3>& camera_pos, bool print = true)
{
    const std::int64_t Nx = 2048;
    const std::int64_t Ny = 2048;

    std::vector<T> buffer(Nx * Ny, 0);
    std::vector<std::size_t> idx(buffer.size());
    std::iota(idx.begin(), idx.end(), 0);

    const auto aabb = object.AABB();
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

        const auto intersection = object.intersect(p);

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
void benchmark(const std::vector<dxmc::Triangle<T>>& tri)
{
    std::array<T, 3> pos { 1000, 1000, 1000 };
    for (std::size_t i = 4; i <= 14; ++i) {
        dxmc::TriangulatedMesh<T> mesh(tri, i);
        std::cout << "Max depth: " << i << " Depth: " << mesh.kdtree().depth();
        auto time = create_image(mesh, pos, false);
        std::cout << " Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(time).count() << std::endl;
    }
}


template <typename T>
std::vector<dxmc::Triangle<T>> getPyramid() {
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
    return t;
}

template <typename T>
std::vector<dxmc::Triangle<T>> getBox()
{
    std::vector<std::array<T, 3>> p;
    constexpr T d = 30;
    p.push_back({ 1, 1, 1 }); // 0
    p.push_back({ 1, 1, -1 }); // 1
    p.push_back({ 1, -1, 1 }); // 2
    p.push_back({ -1, 1, 1 }); // 3
    p.push_back({ -1, -1, 1 }); // 4
    p.push_back({ -1, 1, -1 }); // 5
    p.push_back({ 1, -1, -1 }); // 6
    p.push_back({ -1, -1, -1 }); // 7
    for (auto& i : p)
        for (auto& j : i)
            j *= d;

    std::vector<dxmc::Triangle<T>> t;
    t.push_back({ p[0], p[3], p[4] });
    t.push_back({ p[0], p[4], p[2] });
    t.push_back({ p[6], p[2], p[4] });
    t.push_back({ p[6], p[4], p[7] });
    t.push_back({ p[7], p[4], p[3] });
    t.push_back({ p[7], p[3], p[5] });
    t.push_back({ p[5], p[1], p[6] });
    t.push_back({ p[5], p[6], p[7] });
    t.push_back({ p[1], p[0], p[2] });
    t.push_back({ p[1], p[2], p[6] });
    t.push_back({ p[5], p[3], p[0] });
    t.push_back({ p[5], p[0], p[1] });
    return t;
}
int main(int argc, char* argv[])
{
    //auto reader = dxmc::STLReader<double>("bunny.stl");
    // auto reader = dxmc::STLReader<double>("bunny_low.stl");
    // auto reader = dxmc::STLReader<double>("duck.stl");
    //const auto triangles = reader();
    const auto triangles = getBox<double>();
    //const auto triangles = getPyramid<double>();
    
    dxmc::TriangulatedMesh<double> mesh(triangles, 9);

    // centering aabb
    auto aabb_pre = mesh.AABB();
    auto dist = neg_center(aabb_pre);
    mesh.translate(dist);
    auto aabb = mesh.AABB();

    std::array<double, 3> pos { -1000.0, 500.0, 500.0 };
    create_image(mesh, pos);
    // benchmark(triangles);
    return EXIT_SUCCESS;
}
