

#include "dxmc/material.h"
#include "dxmc/source.h"
#include "dxmc/transport.h"
#include "dxmc/world.h"

#include <array>
#include <iostream>
#include <memory>
#include <numeric>

using namespace dxmc;

int main(int argc, char* argv[])
{
    // Lets create a world that describes our voxelized model to score dose in.
    World world;
    // We need to specify dimensions and voxel spacing (in millimeters)
    std::array<std::size_t, 3> dimensions = { 56, 56, 56 };
    world.setDimensions(dimensions);
    std::array<double, 3> spacing = { 1.0, 1.0, 1.0 };
    world.setSpacing(spacing);

    // Now we fill the world box with som materials
    // We specify three materials
    Material air("Air, Dry (near sea level)"); // Material air("N0.76O0.23Ar0.01") is eqvivalent
    Material water("Water, Liquid"); // Material water("H2O") is eqvivalent
    Material aluminium(13); // aluminium

    // Lets add the materials to the world material map
    world.addMaterialToMap(air);
    world.addMaterialToMap(water);
    world.addMaterialToMap(aluminium);

    //now we create the density array and material index array
    auto arraySize = dimensions[0] * dimensions[1] * dimensions[2];
    // the type of density array is double
    auto densityArray = std::make_shared<std::vector<double>>(arraySize);
    // the type of material index array is unsigned char (up to 255 materials are supported)
    auto materialIndexArray = std::make_shared<std::vector<unsigned char>>(arraySize);
    // Now fill the arrays
    for (std::size_t i = 0; i < dimensions[0]; ++i)
        for (std::size_t j = 0; j < dimensions[1]; ++j)
            for (std::size_t k = 0; k < dimensions[2]; ++k) {
                auto index = i + j * dimensions[0] + k * dimensions[0] * dimensions[1];
                if (k < dimensions[2] / 3) {
                    densityArray->data()[index] = air.standardDensity();
                    ;
                    materialIndexArray->data()[index] = 0;
                } else if (k < dimensions[2] * 2 / 3) {
                    densityArray->data()[index] = water.standardDensity();
                    materialIndexArray->data()[index] = 1;
                } else {
                    densityArray->data()[index] = aluminium.standardDensity();
                    materialIndexArray->data()[index] = 2;
                }
            }

    world.setDensityArray(densityArray);
    world.setMaterialIndexArray(materialIndexArray);

    // we need a beam source
    PencilSource pen;

    //pencilsource only support monochromatic intensity
    pen.setPhotonEnergy(60.0); // keV
    world.setAttenuationLutMaxEnergy(60.0);
    //position the source above the world box
    //std::array<double, 3> source_position = { 0,0, -dimensions[2] * spacing[2] };
    std::array<double, 3> source_position = { 0, 0, 0 };
    pen.setPosition(source_position);
    // we want the cross product aka the beam direction to point towards the box
    std::array<double, 6> beam_cosines = { -1, 0, 0, 0, 1, 0 };
    pen.setDirectionCosines(beam_cosines);
    pen.setTotalExposures(100); // number of batches
    pen.setHistoriesPerExposure(10000); // histories per batch

    //setting normalizing DAP value
    pen.setAirDose(1.0);
    //run simulation
    auto test = world.validate();
    auto res = transport::run(world, &pen);
    auto dose = res.dose;

    //printing dose dose along z direction of box
    auto materials = world.materialMap();
    std::cout << "Depth [mm],  Dose, Material index, Material\n";
    for (std::size_t k = 0; k < dimensions[2]; ++k) {
        std::cout << spacing[2] * k << ", ";
        const auto offset = dimensions[0] * dimensions[1];
        const auto k_offset = dimensions[0] * dimensions[1] * k;
        auto sum_start = dose.begin() + k_offset;
        auto sum_end = dose.begin() + k_offset + offset;
        auto sum = std::accumulate(sum_start, sum_end, 0.0);
        std::cout << sum / static_cast<double>(offset) << ", ";
        auto material_index = static_cast<std::size_t>(materialIndexArray->data()[k_offset]);
        std::cout << material_index << ", " << materials[material_index].prettyName() << std::endl;
    }

    return 1;
}