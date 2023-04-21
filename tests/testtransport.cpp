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

#include "dxmc/beams/isotropicmonoenergybeam.hpp"
#include "dxmc/beams/pencilbeam.hpp"
#include "dxmc/transport.hpp"
#include "dxmc/world/world.hpp"
#include "dxmc/world/worlditems/aavoxelgrid.hpp"
#include "dxmc/world/worlditems/ctdiphantom.hpp"
#include "dxmc/world/worlditems/depthdose.hpp"
#include "dxmc/world/worlditems/worldbox.hpp"
#include "dxmc/world/worlditems/worldcylinder.hpp"
#include "dxmc/world/worlditems/worldsphere.hpp"

#include <iostream>

template <typename T, typename W, typename B>
auto runDispatcher(T& transport, W& world, const B& beam)
{
    dxmc::TransportProgress progress;

    bool running = true;
    std::thread job([&]() {
        transport(world, beam, &progress);
        running = false;
    });
    std::string message;
    while (running) {
        std::this_thread::sleep_for(std::chrono::milliseconds(1000));
        std::cout << std::string(message.length(), ' ') << "\r";
        message = progress.message();
        std::cout << message << "\r";
    }
    job.join();
    std::cout << std::string(message.length(), ' ') << "\r";
    return progress.totalTime();
}

template <typename T>
bool testCylinder()
{
    using Cylinder = dxmc::WorldCylinder<T, 5, 1>;
    //using Cylinder=dxmc::WorldSphere<T,5,1>;
    dxmc::World2<T, Cylinder> world;
    
    //auto& cylinder = world.addItem<Cylinder>({ 4 });
    auto& cylinder = world.addItem<Cylinder>({ 4, 60 });

    world.build(60);

    dxmc::Transport transport;

    std::uint64_t nPart = 0;
    std::chrono::milliseconds time;

    dxmc::IsotropicMonoEnergyBeam<T> beam;

    const auto collAngleY = std::atan(2.0f / 60);
    const auto collAngleZ = std::atan(2.0f / 60);
    beam.setCollimationAngles({ -collAngleY, -collAngleZ, collAngleY, collAngleZ });
    beam.setNumberOfParticlesPerExposure(1e5);
    beam.setNumberOfExposures(32);

    std::vector<dxmc::DoseScore<T>> res(8);
    const T angStep = 360 / res.size();
    for (int ang = 0; ang < res.size(); ++ang) {
        const auto a = dxmc::DEG_TO_RAD<T>() * ang * angStep;
        constexpr std::array<T, 3> pos = { -60, 0, 0 };
        constexpr std::array<T, 3> cosy = { 0, 1, 0 };
        constexpr std::array<T, 3> cosz = { 0, 0, 1 };

        auto rpos = dxmc::vectormath::rotate(pos, { 0, 0, 1 }, a);
        auto rcosy = dxmc::vectormath::rotate(cosy, { 0, 0, 1 }, a);
        auto rcosz = dxmc::vectormath::rotate(cosz, { 0, 0, 1 }, a);
        beam.setPosition(rpos);
        beam.setDirectionCosines(rcosy, rcosz);

        time = runDispatcher(transport, world, beam);
        res[ang] = cylinder.dose();

        cylinder.clearDose();

        std::cout << ang * angStep << ", " << res[ang].energyImparted() << ", " << res[ang].relativeUncertainty() << std::endl;
    }

    return true;
}

template <typename T>
bool testCTDI()
{
    using CTDI = dxmc::CTDIPhantom<T>;
    using Box = dxmc::WorldBox<T, 4, 1>;

    dxmc::World2<T, CTDI, Box> world;
    auto& phantom = world.addItem<CTDI>({});

    phantom.setHoleMaterial("Polymethyl Methacralate (Lucite, Perspex)", dxmc::NISTMaterials<T>::density("Polymethyl Methacralate (Lucite, Perspex)"));

    world.build(30);

    dxmc::Transport transport;

    std::uint64_t nPart = 0;
    std::chrono::milliseconds time;

    dxmc::IsotropicMonoEnergyBeam<T> beam;

    const auto collAngleY = std::atan(16.0f / 60);
    const auto collAngleZ = 0; //    std::atan(4.0f / 60);
    beam.setCollimationAngles({ -collAngleY, -collAngleZ, collAngleY, collAngleZ });
    beam.setNumberOfParticlesPerExposure(1e6);
    beam.setNumberOfExposures(32);

    std::vector<std::array<T, 6>> res(8);
    constexpr T angStep = 45;
    for (int ang = 0; ang < res.size(); ++ang) {
        const auto a = dxmc::DEG_TO_RAD<T>() * ang * angStep;
        constexpr std::array<T, 3> pos = { -60, 0, 0 };
        constexpr std::array<T, 3> cosy = { 0, 1, 0 };
        constexpr std::array<T, 3> cosz = { 0, 0, 1 };

        auto rpos = dxmc::vectormath::rotate(pos, { 0, 0, 1 }, a);
        auto rcosy = dxmc::vectormath::rotate(cosy, { 0, 0, 1 }, a);
        auto rcosz = dxmc::vectormath::rotate(cosz, { 0, 0, 1 }, a);
        beam.setPosition(rpos);
        beam.setDirectionCosines(rcosy, rcosz);

        time = runDispatcher(transport, world, beam);
        for (int i = 0; i < 6; ++i) {
            res[ang][i] = (phantom.dose(i).energyImparted() * 1000) / beam.numberOfParticles();
        }

        std::cout << ang * angStep << ", ";
        for (auto a : res[ang]) {
            std::cout << a << ", ";
        }
        std::cout << std::endl;
        phantom.clearDose();
    }

    for (int i = 0; i < res.size(); ++i) {
        std::cout << i * angStep << ", ";
        for (auto a : res[i]) {
            std::cout << a << ", ";
        }
        std::cout << std::endl;
    }

    return true;
}

template <typename T>
bool testDepth(bool print = false)
{
    using Cylinder = dxmc::DepthDose<T, 5, 2>;
    using World = dxmc::World2<T, Cylinder>;
    using Beam = dxmc::PencilBeam<T>;

    World world;
    auto& cylinder = world.addItem<Cylinder>({ T { 0.1 }, T { 20 }, 20 });
    auto material = dxmc::Material2<T, 5>::byNistName("Polymethyl Methacralate (Lucite, Perspex)");
    if (material) {
        cylinder.setMaterial(material.value());
        cylinder.setMaterialDensity(T { 1.19 });
    }
    world.build();

    const T energy = 60;

    Beam beam;
    beam.setEnergy(energy);
    beam.setPosition({ 0, 0, -100 });
    beam.setDirection({ 0, 0, 1 });

    beam.setNumberOfExposures(40);
    beam.setNumberOfParticlesPerExposure(1e4);

    // beam.setParticleWeight(T { 0.5 });

    dxmc::Transport transport;
    auto time = runDispatcher(transport, world, beam);

    if (print)
        std::cout << "pos, energy, std, nevents\n";

    T dose0 = -1;
    T pos0 = -1;
    const auto att = material.value().attenuationValues(beam.energy()).sum();

    bool success = true;

    for (const auto& [pos, d] : cylinder.depthDose()) {
        if (dose0 < 0) {
            dose0 = d.energyImparted();
            pos0 = pos;
        }
        if (print) {
            std::cout << pos << ", ";
            std::cout << d.energyImparted() << ", ";
            std::cout << d.stdEnergyImparted() << ", ";
            std::cout << d.numberOfEvents() << "\n";
        }
        success = success && d.energyImparted() / dose0 - std::exp(-(pos - pos0) * att) < 0.01;
    }
    return success;
}

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
                const auto xc = x * spacing[0] + spacing[0] / 2 - (dim[0] * spacing[0]) / 2;
                const auto yc = y * spacing[1] + spacing[0] / 2 - (dim[1] * spacing[1]) / 2;
                const auto zc = z * spacing[2] + spacing[0] / 2 - (dim[2] * spacing[2]) / 2;

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

    const std::size_t c0 = 3; // dim_max * 1 / 4;
    const std::size_t c1 = dim_max - 4; // dim_max * 3 / 4;

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

template <typename T, std::uint_fast8_t TRANSPARENT = 0>
bool testAAVoxelGridTransport()
{
    bool success = true;

    using AAVoxelGrid = dxmc::AAVoxelGrid<T, 5, 1, TRANSPARENT>;
    using Cylinder = dxmc::WorldCylinder<T, 5, 2>;
    using Sphere = dxmc::WorldSphere<T, 5, 2>;
    using World = dxmc::World2<T, AAVoxelGrid, Sphere>;

    World world;
    auto& grid = world.addItem(AAVoxelGrid());
    auto& sphere = world.addItem(Sphere(3));

    auto air = dxmc::Material2<T, 5>::byNistName("Air, Dry (near sea level)").value();
    auto pmma = dxmc::Material2<T, 5>::byNistName("Polymethyl Methacralate (Lucite, Perspex)").value();
    const auto air_dens = dxmc::NISTMaterials<T>::density("Air, Dry (near sea level)");
    const auto pmma_dens = dxmc::NISTMaterials<T>::density("Polymethyl Methacralate (Lucite, Perspex)");

    sphere.setMaterial(pmma);
    sphere.setMaterialDensity(pmma_dens);

    const std::array<std::size_t, 3> dim = { 64, 64, 64 };
    const auto size = std::reduce(dim.cbegin(), dim.cend(), size_t { 1 }, std::multiplies<>());
    std::array<T, 3> spacing = { .5f, .5f, .5f };

    // setting up grid
    auto matInd = generateDonut<std::uint8_t>(dim, spacing);
    // auto matInd = generateEdges<std::uint8_t>(dim, spacing);
    std::vector<T> dens(size, air_dens);
    std::transform(matInd.cbegin(), matInd.cend(), dens.begin(), [=](const auto i) { return i == 0 ? air_dens : pmma_dens; });
    std::vector<dxmc::Material2<T, 5>> materials;
    materials.push_back(air);
    materials.push_back(pmma);
    grid.setData(dim, dens, matInd, materials);
    grid.setSpacing(spacing);

    dxmc::PencilBeam<T> beam({ 0, 0, -100 }, { 0, 0, 1 }, 60);
    beam.setNumberOfExposures(64);
    beam.setNumberOfParticlesPerExposure(100000);
    dxmc::Transport transport;
    world.build();
    auto time = runDispatcher(transport, world, beam);
    std::cout << std::format("Total time: {}", time) << std::endl;

    std::vector<T> doseArray(dens.size());
    for (std::size_t i = 0; i < size; ++i) {
        const auto& dose = grid.dose(i);
        doseArray[i] = dose.energyImparted();
    }
    writeImage(doseArray, "dose.bin");

    return success;
}

template <typename T>
bool testAAVoxelGrid()
{

    dxmc::AAVoxelGrid<T, 5> item;

    auto air = dxmc::Material2<T, 5>::byNistName("Air, Dry (near sea level)").value();
    auto pmma = dxmc::Material2<T, 5>::byNistName("Polymethyl Methacralate (Lucite, Perspex)").value();
    const auto air_dens = dxmc::NISTMaterials<T>::density("Air, Dry (near sea level)");
    const auto pmma_dens = dxmc::NISTMaterials<T>::density("Polymethyl Methacralate (Lucite, Perspex)");

    const std::array<std::size_t, 3> dim = { 64, 64, 64 };
    const auto size = std::reduce(dim.cbegin(), dim.cend(), size_t { 1 }, std::multiplies<>());
    std::array<T, 3> spacing = { .20f, .20f, .20f };

    // material arrays
    std::vector<T> dens(size, air_dens);
    std::vector<std::uint8_t> materialIdx(size, 0);
    std::vector<dxmc::Material2<T, 5>> materials;
    materials.push_back(air);
    materials.push_back(pmma);

    // test indices and assign material
    bool success = true;
    item.setData(dim, dens, materialIdx, materials);
    item.setSpacing(spacing);

    for (std::size_t z = 0; z < dim[2]; ++z)
        for (std::size_t y = 0; y < dim[1]; ++y)
            for (std::size_t x = 0; x < dim[0]; ++x) {
                const std::array tind = { x, y, z };
                const auto find = item.flatIndex(tind);
                const auto ind = item.index(find);
                success = success && ind == tind;
                /* if (x > dim[0] / 2 && y > dim[1] / 2 && z > dim[2] / 2) {
                    dens[find] = pmma_dens;
                    materialIdx[find] = 1;
                } else if (x < dim[0] / 2 && y < dim[1] / 2 && z < dim[2] / 2) {
                    dens[find] = pmma_dens;
                    materialIdx[find] = 1;
                }*/
                if (x == 1 && y == 1 && z == 1) {
                    dens[find] = pmma_dens;
                    materialIdx[find] = 1;
                }
            }

    item.setData(dim, dens, materialIdx, materials);
    item.setSpacing(spacing);

    dxmc::Particle<T> p;

    p.pos = { -100, 0, 0 };
    p.dir = { 1, 0, 0 };
    auto res = item.intersect(p);
    success = success && res.valid() && res.intersection == 99;

    p.pos = { 100, 0, 0 };
    p.dir = { -1, 0, 0 };
    res = item.intersect(p);
    success = success && res.valid() && res.intersection == 99;

    p.pos = { -100, -100, -100 };
    p.dir = { 1, 1, 1 };
    dxmc::vectormath::normalize(p.dir);
    res = item.intersect(p);
    auto val = std::sqrt(3 * 100 * T { 100 }) - std::sqrt(3 * spacing[0] * spacing[0] / 4);
    success = success && res.valid() && val - res.intersection < 1E-6;

    return success;
}

int main()
{
    bool success = true;
    success = success && testCylinder<double>();
    // success = success && testCTDI<double>();

    // success = success && testAAVoxelGridTransport<double, 0>();
    //  success = success && testAAVoxelGridTransport<float, 0>();
    /* success = success && testAAVoxelGridTransport<float, 255>();

    success = success && testAAVoxelGrid<float>();
    success = success && testAAVoxelGrid<double>();

    success = success && testDepth<float>();
    success = success && testDepth<double>();

    success = success && testTransport<float>();
    success = success && testTransport<float>();
    success = success && testTransport<float>();
    success = success && testTransport<double>();
    success = success && testTransport<double>();
    success = success && testTransport<double>();
    */
    if (success)
        return EXIT_SUCCESS;

    return EXIT_FAILURE;
}