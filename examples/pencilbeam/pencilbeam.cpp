

#include "dxmc/material.hpp"
#include "dxmc/source.hpp"
#include "dxmc/transport.hpp"
#include "dxmc/world.hpp"

#include <array>
#include <iostream>
#include <memory>
#include <numeric>

using namespace dxmc;

template <typename T>
double calculate()
{
    // Lets create a world that describes our voxelized model to score dose in.
    World<T> world;
    // We need to specify dimensions and voxel spacing (in millimeters)
    std::array<std::size_t, 3> dimensions = { 56, 56, 56 };
    world.setDimensions(dimensions);
    std::array<T, 3> spacing = { 1.0, 1.0, 1.0 };
    world.setSpacing(spacing);

    // Now we fill the world box with some materials
    // We specify three materials
    Material air("Air, Dry (near sea level)"); // Material air("N0.76O0.23Ar0.01") is equivalent
    Material water("Water, Liquid"); // Material water("H2O") is equivalent
    Material aluminium(13); // Material aluminum

    // Lets add the materials to the world material map
    world.addMaterialToMap(air);
    world.addMaterialToMap(water);
    world.addMaterialToMap(aluminium);

    // now we create the density array and material index array
    auto arraySize = world.size();
    // the type of density array is double
    auto densityArray = std::make_shared<std::vector<T>>(arraySize);
    // the type of material index array is std::uint8_t or old school unsigned char (up to 255 materials are supported)
    auto materialIndexArray = std::make_shared<std::vector<std::uint8_t>>(arraySize);
    // Now fill the arrays
    for (std::size_t i = 0; i < dimensions[0]; ++i)
        for (std::size_t j = 0; j < dimensions[1]; ++j)
            for (std::size_t k = 0; k < dimensions[2]; ++k) {
                auto index = i + j * dimensions[0] + k * dimensions[0] * dimensions[1];
                if (k < dimensions[2] / 3) { // The first 1/3 part of the box is air.
                    densityArray->data()[index] = air.standardDensity();
                    materialIndexArray->data()[index] = 0;
                } else if (k < dimensions[2] * 2 / 3) { // The second 1/3 part of the box is water.
                    densityArray->data()[index] = water.standardDensity();
                    materialIndexArray->data()[index] = 1;
                } else { // The third 1/3 part of the box is aluminum.
                    densityArray->data()[index] = aluminium.standardDensity();
                    materialIndexArray->data()[index] = 2;
                }
            }
    // Add the density array and the material index array to world.
    world.setDensityArray(densityArray);
    world.setMaterialIndexArray(materialIndexArray);
    // The world is now complete, to test if it is valid and all material indices are correct
    // and arrays have correct dimensions we can call:
    world.makeValid();

    // We need a beam source, in this case a simple pencil beam
    PencilSource<T> pen;

    // PencilSource only support monochromatic intensity, in this case we use 60 keV initial photon energy
    pen.setPhotonEnergy(60.0); // keV

    // Position the source above the world box
    std::array<T, 3> source_position = { 0, 0, -(dimensions[2] * spacing[2]) };
    pen.setPosition(source_position);
    // We want the source plane normal aka the beam direction to point towards the box
    // To set the direction of the beam we must specify the beam direction cosines.
    // This is basically two orthonormal vectors on the beam plane [x, y].
    std::array<T, 6> beam_cosines = { 1, 0, 0, 0, 1, 0 };
    pen.setDirectionCosines(beam_cosines);
    // The beam direction is then x cross y
    std::array<T, 3> beam_direction;
    vectormath::cross(beam_cosines.data(), beam_direction.data());
    // The beam direction should be [0, 0, 1]

    // Specify number of exposures, each exposure is run on a single thread
    // so make sure there are more exposures than cores on your CPU for optimal performance.
    // For CT sources, number of exposures are determined by the geometry.
    pen.setTotalExposures(20);
    // Set total number of histories per exposure
    pen.setHistoriesPerExposure(1000000);

    // Create a transport object
    Transport<T> transport;
    // Do the simulation, this might take some time depending on total number of histories.
    auto res = transport(world, &pen);

    // Print the dose along z direction of box
    const auto& materials = world.materialMap();
    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(res.simulationTime).count() / 1000.0;
    std::cout << "Simulation time, " << time << ", seconds\n";
    std::cout << "Depth [mm],  Dose, nEvents, Material index, Material\n";
    auto dose = res.dose;
    for (std::size_t k = 0; k < dimensions[2]; ++k) {
        std::cout << spacing[2] * k << ", ";
        const auto offset = dimensions[0] * dimensions[1];
        const auto k_offset = dimensions[0] * dimensions[1] * k;
        auto sum_start = dose.begin() + k_offset;
        auto sum_end = dose.begin() + k_offset + offset;
        auto sum = std::accumulate(sum_start, sum_end, 0.0);
        std::cout << sum / static_cast<T>(offset) << ", ";
        auto n_events = std::accumulate(res.nEvents.cbegin() + k_offset, res.nEvents.cbegin() + k_offset + offset, 0);
        std::cout << n_events / static_cast<T>(offset) << ", ";
        auto material_index = static_cast<std::size_t>(materialIndexArray->data()[k_offset]);
        std::cout << material_index << ", " << materials[material_index].prettyName() << std::endl;
    }
    return time;
}

int main(int argc, char* argv[])
{
    auto n_secs_d = calculate<double>();
    auto n_secs_f = calculate<float>();
    std::cout << "Simulation time for double precision: ";
    std::cout << n_secs_d << " seconds\n";
    std::cout << "Simulation time for single precision: ";
    std::cout << n_secs_f << " seconds\n";
    std::cout << "Single precision was " << std::ceil(100 - n_secs_f / n_secs_d * 100);
    std::cout << " percent faster\n";
    return 1;
}