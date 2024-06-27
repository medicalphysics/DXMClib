

#include "dxmc/beams/pencilbeam.hpp"
#include "dxmc/material/material.hpp"
#include "dxmc/transport.hpp"
#include "dxmc/world/visualization/visualizeworld.hpp"
#include "dxmc/world/world.hpp"
#include "dxmc/world/worlditems/depthdose.hpp"

#include <algorithm>
#include <iostream>
#include <thread>

void example()
{
    // Create a thin aluminium cylinder and print depth dose
    std::cout << "Example Pencilbeam\nTransport of monoenergetic photons in a thin long aluminium cylinder\n";

    // Start creating a depthscore cylinder
    constexpr int N_ATOMIC_SHELLS = 5; // Number of atomic shells to consider binding energies.
    constexpr int LOW_ENERGY_CORRECTION = 1; /* 0: No binding energy correction, 1: Livermore corection, 2; impulse approx. correction*/

    using Cylinder = dxmc::DepthDose<N_ATOMIC_SHELLS, LOW_ENERGY_CORRECTION>;

    // Create a world that can consist of one or more depthdose objects
    dxmc::World<Cylinder> world;
    // Reserve number of items in world, not nessecary here but it's good practice.
    world.reserveNumberOfItems(1);

    // Adding a depthdose object
    auto& cylinder = world.template addItem<Cylinder>({ 1 /* radius */, 10 /* lenght */ });
    // Set material and density
    auto aluminium = dxmc::Material<N_ATOMIC_SHELLS>::byZ(13).value();
    cylinder.setMaterial(aluminium, 2.27 /* g/cm3 */);
    // Example for constructing other materials
    //    auto concrete = dxmc::Material<5>::byNistName("Concrete, Ordinary").value();
    //    auto concrete_density = dxmc::NISTMaterials::density("Concrete, Ordinary");
    //    cylinder.setMaterial(concrete, concrete_density);

    std::cout << "with radius " << cylinder.radius() << " cm and height " << cylinder.length() << " cm\n";

    // Building world
    world.build();

    // Define a radiation source
    dxmc::PencilBeam<> beam({ 0, 0, -10 } /* position */, { 0, 0, 1 } /* direction */);
    beam.setNumberOfExposures(64); // number of jobs
    beam.setNumberOfParticlesPerExposure(1000000); // histories per job

    // Lunch simulation
    auto nThreads = std::max(std::thread::hardware_concurrency(), std::uint32_t { 1 });
    auto time_elapsed = dxmc::Transport::runConsole(world, beam, nThreads);

    // Get max dose and print some values
    std::cout << "Depth dose in cylinder for " << beam.numberOfParticles() << " photons of " << beam.energy() << " keV\n";
    std::cout << "Simulation time: " << time_elapsed << std::endl;
    std::cout << "Depth [cm], Energy scored [keV], Dose per Air Kerma [mGy/mGy], Uncertainty [%], #Events\n";
    double max_dose = 0;
    for (std::size_t i = 0; i < cylinder.resolution(); ++i) {
        auto dose_val = cylinder.doseScored(i).dose();
        max_dose = std::max(max_dose, dose_val);
        std::cout << (cylinder.length() * (i + 0.5)) / cylinder.resolution() << ", ";
        std::cout << cylinder.energyScored(i).energyImparted() << ", ";
        std::cout << cylinder.doseScored(i).dose() << ", ";
        std::cout << cylinder.doseScored(i).relativeUncertainty() * 100 << ", ";
        std::cout << cylinder.doseScored(i).numberOfEvents() << std::endl;
    }

    // Save some images
    dxmc::VisualizeWorld viz(world);
    auto buffer = viz.template createBuffer<double>(1024, 1024);
    viz.setDistance(60);
    viz.setAzimuthalAngleDeg(90);
    viz.setPolarAngleDeg(30);
    viz.suggestFOV(2); // zoom = 2
    viz.setColorOfItem<std::uint8_t>(world.getItemPointers()[0], { 255, 192, 203 }); // making it pink
    viz.generate(world, buffer);
    viz.savePNG("cylinder.png", buffer);

    viz.setColorByValueMinMax(-0.01, max_dose); // color by dose value
    viz.addColorByValueItem(world.getItemPointers()[0]);
    viz.generate(world, buffer);
    viz.savePNG("cylinder_dose.png", buffer);
}

int main(int argc, char* argv[])
{
    example();
    return 1;
}