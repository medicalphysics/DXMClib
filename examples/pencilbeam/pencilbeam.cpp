

#include "dxmc/beams/pencilbeam.hpp"
#include "dxmc/material/material.hpp"
#include "dxmc/transport.hpp"
#include "dxmc/world/world.hpp"
#include "dxmc/world/worlditems/depthdose.hpp"

#include <algorithm>
#include <iostream>
#include <thread>

auto calculate()
{
    // Lets create a world that describes our voxelized model to score dose in.

    // Start creating a depthscore cylinder
    using Cylinder = dxmc::DepthDose<>;

    dxmc::World<Cylinder> world;
    world.reserveNumberOfItems(8);

    auto& cylinder = world.template addItem<Cylinder>({ 1, 50 });

    // We need to specify dimensions and voxel spacing (in millimeters)

    // Now we fill the world box with some materials
    // We specify three materials
    auto aluminium = dxmc::Material<>::byZ(13).value(); // air("Air, Dry (near sea level)"); // Material air("N0.76O0.23Ar0.01") is equivalent
    // Material water("Water, Liquid"); // Material water("H2O") is equivalent
    // Material aluminium(13); // Material aluminum

    cylinder.setMaterial(aluminium, 2.27 /*g/cm3*/);

    world.build(100);

    dxmc::PencilBeam<> beam({ 0, 0, -100 }, { 0, 0, 1 });
    beam.setNumberOfExposures(64);
    beam.setNumberOfParticlesPerExposure(10000);

    auto nThreads = std::max(std::thread::hardware_concurrency(), std::uint32_t { 1 });

    auto time_elapsed = dxmc::Transport::runConsole(world, beam, nThreads);

    return time_elapsed;
}

int main(int argc, char* argv[])
{
    auto n_secs_d = calculate();
    std::cout << "Simulation time: ";
    std::cout << n_secs_d << " seconds\n";
    return 1;
}