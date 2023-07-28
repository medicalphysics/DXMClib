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

#include "dxmc/world/visualization/visualizeworld.hpp"
#include "dxmc/world/world.hpp"
#include "dxmc/world/worlditems/aavoxelgrid.hpp"
#include "dxmc/world/worlditems/worldcylinder.hpp"

#include <iostream>
#include <string>
#include <vector>

template <typename T>
void writeImage(const std::vector<T>& buffer, const std::string& name)
{
    std::ofstream file;
    file.open(name, std::ios::out | std::ios::binary);
    file.write((char*)buffer.data(), buffer.size() * sizeof(T));
    file.close();
}

template <typename T = int, typename U>
std::vector<T> generateDonut(const std::array<std::size_t, 3>& dim, const std::array<U, 3>& spacing = { 1, 1, 1 })
{
    const auto s = std::reduce(dim.cbegin(), dim.cend(), std::size_t { 1 }, std::multiplies<>());
    std::vector<T> d(s, T { 0 });

    const U R = U { 0.25 } * (dim[0] * spacing[0] + dim[1] * spacing[1]) / 2;
    const U r = U { 0.1 } * (dim[0] * spacing[0] + dim[1] * spacing[1]) / 2;

    for (std::size_t z = 0; z < dim[2]; ++z)
        for (std::size_t y = 0; y < dim[1]; ++y)
            for (std::size_t x = 0; x < dim[0]; ++x) {
                const auto flat_ind = x + y * dim[0] + z * dim[0] * dim[1];
                const auto xc = x * spacing[0] - (dim[0] * spacing[0]) / 2 + spacing[0] / 2;
                const auto yc = y * spacing[1] - (dim[1] * spacing[1]) / 2 + spacing[1] / 2;
                const auto zc = z * spacing[2] - (dim[2] * spacing[2]) / 2 + spacing[2] / 2;

                const auto p1 = R - std::sqrt(xc * xc + yc * yc);
                if (p1 * p1 + zc * zc < r * r)
                    d[flat_ind] = 1;
            }
    return d;
}

template <typename T = int>
std::vector<T> generateEdges(const std::array<std::size_t, 3>& dim, const std::array<double, 3>& spacing = { 1, 1, 1 })
{
    const auto s = std::reduce(dim.cbegin(), dim.cend(), std::size_t { 1 }, std::multiplies<>());
    std::vector<T> d(s, T { 0 });

    const auto dim_max = std::max(dim[0], std::max(dim[1], dim[2]));

    const std::size_t c0 = 1; // dim_max * 1 / 4;
    const std::size_t c1 = dim_max - 2; // dim_max * 3 / 4;

    for (std::size_t z = 0; z < dim[2]; ++z)
        for (std::size_t y = 0; y < dim[1]; ++y)
            for (std::size_t x = 0; x < dim[0]; ++x) {
                const auto flat_ind = x + y * dim[0] + z * dim[0] * dim[1];
                if (c0 <= x && x <= c1 && (y == c0 || y == c1) && (z == c0 || z == c1)) {
                    d[flat_ind] = 1;
                }
                if (c0 <= y && y <= c1 && (x == c0 || x == c1) && (z == c0 || z == c1)) {
                    d[flat_ind] = 1;
                }
                if (c0 <= z && z <= c1 && (y == c0 || y == c1) && (x == c0 || x == c1)) {
                    d[flat_ind] = 1;
                }
            }
    return d;
}

template <typename T>
bool testGeometryColor()
{
    std::array<std::size_t, 3> dim = { 128, 128, 128 };
    std::array<T, 3> spacing = { .1, .1, .1 };

    using Grid = dxmc::AAVoxelGrid<T, 5, 2, 0>;
    using Cylinder = dxmc::WorldCylinder<T, 5, 2>;
    using World = dxmc::World2<T, Grid, Cylinder>;

    World world;
    auto& grid = world.addItem<Grid>({});
    auto& cylinder = world.addItem<Cylinder>({});
    cylinder.setRadius(1);
    cylinder.setHeight(5);

    auto air = dxmc::Material<T, 5>::byNistName("Air, Dry (near sea level)").value();
    auto pmma = dxmc::Material<T, 5>::byNistName("Polymethyl Methacralate (Lucite, Perspex)").value();
    const auto air_dens = dxmc::NISTMaterials<T>::density("Air, Dry (near sea level)");
    const auto pmma_dens = dxmc::NISTMaterials<T>::density("Polymethyl Methacralate (Lucite, Perspex)");

    const auto matIdx = generateDonut<std::uint8_t>(dim, spacing);
    // const auto matIdx = generateEdges<std::uint8_t>(dim, spacing);
    std::vector<T> dens(matIdx.size(), 0);
    std::transform(std::execution::par_unseq, matIdx.cbegin(), matIdx.cend(), dens.begin(), [=](const auto i) { return i == 0 ? air_dens : pmma_dens; });

    std::vector<dxmc::Material<T>> materials;
    materials.push_back(air);
    materials.push_back(pmma);

    grid.setData(dim, dens, matIdx, materials);
    grid.setSpacing(spacing);

    world.build(T { 0 });

    dxmc::VisualizeWorld<T> viz(world);
    viz.setPolarAngle(std::numbers::pi_v<T> / 4);
    viz.setAzimuthalAngle((std::numbers::pi_v<T> * 2) / 4 + T {0.1});
    int height = 1024;
    int width = 1024;
    std::vector<T> buffer(height * width * 4, T { 1 });
    viz.generate(world, buffer, width, height);

    writeImage(buffer, "color.bin");
    return false;
}

int main()
{

    bool success = false;
    testGeometryColor<double>();
    // testGeometryColor<float>();
    // testGeometryDistance<double>();

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}