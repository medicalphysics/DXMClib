

#include "dxmc/beams/pencilbeam.hpp"
#include "dxmc/material/material.hpp"
#include "dxmc/transport.hpp"
#include "dxmc/world/visualization/visualizeworld.hpp"
#include "dxmc/world/world.hpp"
#include "dxmc/world/worlditems/depthdose.hpp"
#include "dxmc/world/worlditems/enclosedroom.hpp"

#include <algorithm>
#include <iostream>
#include <thread>

auto calculate()
{
    // Create a thin aluminium cylinder and .

    // Start creating a depthscore cylinder
    using Cylinder = dxmc::DepthDose<5, 1>;
    using Room = dxmc::EnclosedRoom<5, 1>;

    dxmc::World<Cylinder, Room> world;
    world.reserveNumberOfItems(2);

    auto& cylinder = world.template addItem<Cylinder>({ 1, 10 });
    auto aluminium = dxmc::Material<5>::byZ(13).value(); // air("Air, Dry (near sea level)"); // Material air("N0.76O0.23Ar0.01") is equivalent
    cylinder.setMaterial(aluminium, 2.27 /*g/cm3*/);

    auto concrete = dxmc::Material<5>::byNistName("Concrete, Ordinary").value();
    auto concrete_density = dxmc::NISTMaterials::density("Concrete, Ordinary");
    auto& room = world.template addItem<Room>({ 10, 300 });
    room.setMaterial(concrete, concrete_density);

    world.build();

    dxmc::PencilBeam<> beam({ 0, 0, -10 }, { 0, 0, 1 });
    beam.setNumberOfExposures(12);
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
    viz.setDistance(60);
    viz.setAzimuthalAngleDeg(90);
    viz.setPolarAngleDeg(30);
    viz.suggestFOV(5);

    viz.addLineProp(beam, 50, .2);
    viz.generate(world, buffer);
    viz.savePNG("cylinder.png", buffer);

    viz.setColorByValueMinMax(-0.01, max_dose);
    viz.addColorByValueItem(world.getItemPointers()[0]);
    viz.generate(world, buffer);
    viz.savePNG("cylinder_dose.png", buffer);

    return time_elapsed;
}

int main(int argc, char* argv[])
{
    auto n_secs_d = calculate();
    std::cout << "Simulation time: ";
    std::cout << n_secs_d << " seconds\n";
    return 1;
}