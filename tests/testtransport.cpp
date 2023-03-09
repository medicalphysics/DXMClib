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
    dxmc::CTDIPhantom<T> phantom;
    dxmc::WorldBox<T, 4, 1> box(10);
    box.translate({ 0, -50, 0 });
    box.setNistMaterial("Water, Liquid");

    dxmc::World2<T, dxmc::CTDIPhantom<T>, dxmc::WorldBox<T, 4, 1>> world;
    world.addItem(phantom);

    world.addItem(box);
    world.build();

    dxmc::Transport<T> transport;
    transport.setNumberOfThreads(4);

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
        beam.setNumberOfExposures(400);
        nPart = beam.numberOfParticles();
        time = runDispatcher(transport, world, beam);
        // transport(world, beam);
    }

    std::cout << "SUCCESS Simple transport test ";
    std::cout << nPart << " histories in " << dxmc::TransportProgress::human_time(time) << " ";
    std::cout << (nPart*1000) / std::chrono::duration_cast<std::chrono::milliseconds>(time).count() << " histories/sec for sizeof(T) = " << sizeof(T) << std::endl;

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

int main()
{
    bool success = true;
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