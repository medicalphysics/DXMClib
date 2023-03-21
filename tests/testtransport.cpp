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

    dxmc::Transport<T> transport;
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

    dxmc::Transport<T> transport;
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

int main()
{
    bool success = true;

    success = success && testDepth<float>();
    success = success && testDepth<double>();

    success = success && testTransport<float>();
    success = success && testTransport<float>();
    success = success && testTransport<float>();
    success = success && testTransport<double>();
    success = success && testTransport<double>();
    success = success && testTransport<double>();

    if (success)
        return 0;

    return 1;
}