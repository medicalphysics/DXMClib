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

#include "dxmc/beams/isotropicbeam.hpp"
#include "dxmc/beams/pencilbeam.hpp"
#include "dxmc/transport.hpp"
#include "dxmc/tube.hpp"
#include "dxmc/world/visualization/visualizeworld.hpp"
#include "dxmc/world/world.hpp"
#include "dxmc/world/worlditems/fluencescore.hpp"
#include "dxmc/world/worlditems/worldbox.hpp"

#include <iostream>
#include <string>
#include <vector>


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
        std::this_thread::sleep_for(std::chrono::milliseconds(500));
        std::cout << std::string(message.length(), ' ') << "\r";
        message = progress.message();
        std::cout << message << "\r";
    }
    job.join();
    std::cout << std::string(message.length(), ' ') << "\r";
    return progress.totalTime();
}

template <typename T, std::size_t N = 5, int L = 2>
void testfluencescore()
{
    using Box = dxmc::WorldBox<T, N, L>;
    using FluenceScore = dxmc::FluenceScore<T>;
    using World = dxmc::World2<T, Box, FluenceScore>;
    using Material = dxmc::Material2<T, N>;

    World world;
    {
        auto& box = world.template addItem<Box>({ 10 });
        auto water = Material::byNistName("Water, Liquid").value();
        const T density = 1;
        box.setMaterial(water, density);
        world.template addItem<FluenceScore>({ 20, { 0, 0, 15 }, { 0, 0, 1 } });
    }
    world.build();

    dxmc::VisualizeWorld<T> viz(world);
    int height = 512;
    int width = 512;
    std::vector<T> buffer(height * width * 4, T { 1 });
    for (std::size_t i = 0; i < 12; ++i) {
        viz.setDistance(500);
        viz.setPolarAngle(std::numbers::pi_v<T> / 3.0);
        viz.setAzimuthalAngle(std::numbers::pi_v<T> * i / 6.0);
        // viz.setCameraPosition({ -60, -30, -10 });
        viz.suggestFOV();
        viz.generate(world, buffer, width, height);
        std::string name = "fluence_score" + std::to_string(i) + ".png";
        viz.savePNG(name, buffer, width, height);
    }

    dxmc::Tube<T> tube;
    tube.setVoltage(140);
    tube.setAlFiltration(8);
    tube.setEnergyResolution(1);
    const auto ts = tube.getSpecter();

    dxmc::IsotropicBeam<T> beam({ -1000, 0, 0 }, { 0, 1, 0, 0, 0, 1 });
    beam.setEnergySpecter(ts);
    beam.setCollimationAngles(0, 0, 0, 0);

    beam.setNumberOfParticlesPerExposure(1E6);
    beam.setNumberOfExposures(100);

    dxmc::Transport transport;
    runDispatcher(transport, world, beam);
    // transport(world, beam);
    // transport.setNumberOfThreads(4);

    const auto& fluence_vec = world.template getItems<FluenceScore>();
    const auto& fluence = fluence_vec[0];
    const auto specter = fluence.getSpecter();
    const auto N = beam.numberOfExposures() * beam.numberOfParticles();

    std::cout << "Energy, IntensityIn, Energy, IntensityOut\n";
    for (std::size_t i = 0; i < ts.size(); ++i) {
        std::cout << ts[i].first << ", " << N * ts[i].second << ", ";
        std::cout << specter[i + 1].first << ", " << specter[i + 1].second << std::endl;
    }
}

int main()
{
    std::cout << "Testing fluence scoring\n";
    bool success = true;

    testfluencescore<double, 5, 2>();
    if (success)
        std::cout << "SUCCESS ";
    else
        std::cout << "FAILURE ";

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}