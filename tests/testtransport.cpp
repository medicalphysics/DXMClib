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
bool testTransport()
{
    using CTDI = dxmc::CTDIPhantom<T>;
    using Box = dxmc::WorldBox<T, 4, 1>;

    dxmc::World2<T, CTDI, Box> world;
    auto& phantom = world.addItem<CTDI>({});
    auto& box = world.addItem<Box>({ T { 10 } });

    box.translate({ 0, -50, 0 });
    box.setNistMaterial("Water, Liquid");

    world.build();

    dxmc::Transport transport;
    // transport.setNumberOfThreads(4);

    std::uint64_t nPart = 0;
    std::chrono::milliseconds time;
    if constexpr (false) {
        dxmc::PencilBeam<T> beam;
        beam.setPosition({ 0, -1000, 0 });
        beam.setDirection({ 0, 1, 0 });
        beam.setNumberOfParticlesPerExposure(1e3);
        beam.setNumberOfExposures(4);
        nPart = beam.numberOfParticles();
        time = runDispatcher(transport, world, beam);
    } else {
        dxmc::IsotropicMonoEnergyBeam<T> beam;
        beam.setPosition({ 0, -1000, 0 });
        beam.setDirectionCosines({ 0, 0, 1, 1, 0, 0 });
        beam.setCollimationAngles({ .01, .01 });
        beam.setNumberOfParticlesPerExposure(1e5);
        beam.setNumberOfExposures(4);
        nPart = beam.numberOfParticles();
        time = runDispatcher(transport, world, beam);
    }

    std::cout << "SUCCESS Simple transport test ";
    std::cout << nPart << " histories in " << dxmc::TransportProgress::human_time(time) << " ";
    std::cout << (nPart * 1000) / std::chrono::duration_cast<std::chrono::milliseconds>(time).count() << " histories/sec for sizeof(T) = " << sizeof(T) << std::endl;

    const auto& ctdi_vec = world.template getItems<dxmc::CTDIPhantom<T>>();
    auto dose_ctdi = ctdi_vec[0].dose(0);
    std::cout << "CTDI: " << dose_ctdi.energyImparted();
    std::cout << ", " << dose_ctdi.numberOfEvents();
    std::cout << ", " << dose_ctdi.energyImpartedSquared();
    std::cout << std::endl;

    const auto& box_vec = world.template getItems<dxmc::WorldBox<T, 4, 1>>();
    auto dose_box = box_vec[0].dose(0);
    std::cout << "Box: " << dose_box.energyImparted();
    std::cout << ", " << dose_box.numberOfEvents();
    std::cout << ", " << dose_box.energyImpartedSquared();
    std::cout << std::endl;

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

template <typename T, std::uint_fast8_t TRANSPARENT = 255>
bool testAAVoxelGridTransport()
{
    bool success = true;

    using AAVoxelGrid = dxmc::AAVoxelGrid<T, 5, 2, TRANSPARENT>;

    dxmc::World2<T, AAVoxelGrid> world;
    auto& item = world.addItem(AAVoxelGrid());

    auto air = dxmc::Material2<T, 5>::byNistName("Air, Dry (near sea level)").value();
    auto pmma = dxmc::Material2<T, 5>::byNistName("Polymethyl Methacralate (Lucite, Perspex)").value();
    const auto air_dens = dxmc::NISTMaterials<T>::density("Air, Dry (near sea level)");
    //const auto pmma_dens = dxmc::NISTMaterials<T>::density("Polymethyl Methacralate (Lucite, Perspex)");
    const auto pmma_dens = 1;


    const std::array<std::size_t, 3> dim = { 100, 13, 13 };
    const auto size = std::reduce(dim.cbegin(), dim.cend(), size_t { 1 }, std::multiplies<>());
    std::array<T, 3> spacing = { .20f, .20f, .20f };

    // material arrays
    std::vector<T> dens(size, air_dens);
    std::vector<std::uint8_t> materialIdx(size, 0);
    std::vector<dxmc::Material2<T, 5>> materials;
    materials.push_back(air);
    materials.push_back(pmma);

    // test indices and assign material

    item.setSpacing(spacing);
    std::size_t i = 0;
    for (std::size_t z = 0; z < dim[2]; ++z) {
        for (std::size_t y = 0; y < dim[1]; ++y) {
            for (std::size_t x = 0; x < dim[0]; ++x) {

                dens[i] = pmma_dens;
                materialIdx[i] = 1;
                i++;
            }
        }
    }
    materialIdx[0] = 0;

    item.setData(dim, dens, materialIdx, materials);
    item.setSpacing(spacing);

    dxmc::PencilBeam<T> beam({ -100,0,0 }, { 1,0,0 }, 60);
    beam.setNumberOfExposures(40);
    beam.setNumberOfParticlesPerExposure(1000);
    dxmc::Transport transport;
    //transport.setNumberOfThreads(1);
    world.build();
    auto time = runDispatcher(transport, world, beam);
    std::cout << std::format("Total time: {}", time) << std::endl;

    for (std::size_t x = 0; x < dim[0]; ++x) {
        const auto y = dim[1] / 2;
        const auto z = dim[2] / 2;
        const auto ind = item.flatIndex({ x, y, z });
        const auto d = item.dose(ind).energyImparted();
        std::cout << (x + T { 0.5 }) * spacing[0] << ", " << d << "," << static_cast<int>(materialIdx[ind]) << std::endl;
    }

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

    const std::array<std::size_t, 3> dim = { 3, 3, 3 };
    const auto size = std::reduce(dim.cbegin(), dim.cend(), size_t { 1 }, std::multiplies<>());
    std::array<T, 3> spacing = { 2.0f, 2.0f, 2.0f };

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

    success = success && testAAVoxelGridTransport<double, 0>();
    success = success && testAAVoxelGridTransport<double, 255>();
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
        return 0;

    return 1;
}