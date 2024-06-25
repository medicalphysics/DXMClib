

#include "dxmc/beams/pencilbeam.hpp"
#include "dxmc/material/material.hpp"
#include "dxmc/transport.hpp"
#include "dxmc/world/visualization/visualizeworld.hpp"
#include "dxmc/world/world.hpp"
#include "dxmc/world/worlditems/depthdose.hpp"

#include <algorithm>
#include <iostream>
#include <thread>

auto calculate()
{
    // Create a thin aluminium cylinder and .

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

    world.build(1);

    /*dxmc::VisualizeWorld viz(world);
    auto buffer = viz.template createBuffer<double>(1024 , 1024 );
    viz.setDistance(400);
    viz.setAzimuthalAngleDeg(90);
    viz.setPolarAngleDeg(30);
    viz.suggestFOV();
    viz.generate(world, buffer);
    viz.savePNG("dose.png", buffer);
*/

    dxmc::PencilBeam<> beam({ 0, 0, -100 }, { 0, 0, 1 });
    beam.setNumberOfExposures(24);
    beam.setNumberOfParticlesPerExposure(100000);

    auto nThreads = std::max(std::thread::hardware_concurrency(), std::uint32_t { 1 });
    auto time_elapsed = dxmc::Transport::runConsole(world, beam, nThreads);

    // get max dose
    double max_dose = 0;
    for (std::size_t i = 0; i < cylinder.resolution(); ++i) {
        auto dose_val = cylinder.doseScored(i).dose();
        std::cout << dose_val << std::endl;
        max_dose = std::max(max_dose, dose_val);
    }

    dxmc::VisualizeWorld viz(world);
    auto buffer = viz.template createBuffer<double>(1024, 1024);
    viz.setDistance(400);
    viz.setAzimuthalAngleDeg(90);
    viz.setPolarAngleDeg(30);
    viz.suggestFOV();
    viz.setColorByValueMinMax(0, 50000.0);
    viz.addColorByValueItem(world.getItemPointers()[0]);
    //  viz.addLineProp(beam, 50, .2);
    viz.generate(world, buffer);
    viz.savePNG("dose.png", buffer);

    return time_elapsed;
}

int main(int argc, char* argv[])
{
    auto n_secs_d = calculate();
    std::cout << "Simulation time: ";
    std::cout << n_secs_d << " seconds\n";
    return 1;
}