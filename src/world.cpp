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

Copyright 2019 Erlend Andersen
*/

/*
#include "dxmc/world.h"
#include "dxmc/material.h"
#include "dxmc/vectormath.h"
#include <algorithm>
#include <execution>
#include <future>

namespace dxmc {

constexpr double AIRDENSITY = 0.001205; // g/cm3


std::vector<std::size_t> circleIndices2D(const std::array<std::size_t, 2>& dim, const std::array<double, 2>& spacing, const std::array<double, 2>& center, double radius)
{
    std::vector<std::size_t> indices;
    const auto r_int = static_cast<std::size_t>(std::ceil(radius));
    indices.reserve(r_int * r_int);

    std::array<int, 4> flimits;
    flimits[0] = static_cast<int>((center[0] - radius) / spacing[0]);
    flimits[1] = static_cast<int>((center[0] + radius) / spacing[0]) + 1;
    flimits[2] = static_cast<int>((center[1] - radius) / spacing[1]);
    flimits[3] = static_cast<int>((center[1] + radius) / spacing[1]) + 1;
    std::array<std::size_t, 4> limits;
    limits[0] = static_cast<std::size_t>(std::max(flimits[0], 0));
    limits[1] = static_cast<std::size_t>(std::min(flimits[1], static_cast<int>(dim[0])));
    limits[2] = static_cast<std::size_t>(std::max(flimits[2], 0));
    limits[3] = static_cast<std::size_t>(std::min(flimits[3], static_cast<int>(dim[1])));

    const double r2 = radius * radius;
    for (std::size_t i = limits[0]; i < limits[1]; ++i) {
        const double posX = center[0] - i * spacing[0];
        for (std::size_t j = limits[2]; j < limits[3]; ++j) {
            const double posY = center[1] - j * spacing[1];
            if ((posX * posX + posY * posY) <= r2)
                indices.push_back(i + j * dim[0]);
        }
    }
    return indices;
}

CTDIPhantom::CTDIPhantom(std::size_t diameter)
    : World()
{
    setDirectionCosines(1.0, 0.0, 0.0, 0.0, 1.0, 0.0);
    setSpacing(1.0, 1.0, 2.5);
    setOrigin(0.0, 0.0, 0.0);
    if (diameter % 2 == 0) // make sure odd dimensions in xy
        setDimensions(diameter + 3, diameter + 3, 60);
    else
        setDimensions(diameter + 2, diameter + 2, 60);

    Material air("Air, Dry (near sea level)");
    Material pmma("Polymethyl Methacralate (Lucite, Perspex)");
    addMaterialToMap(air);
    addMaterialToMap(air); // we want two to differentiate between air around phantom and air inside rods
    addMaterialToMap(pmma);

    constexpr double holeDiameter = 13.1;
    constexpr double holeRadii = holeDiameter / 2.0;
    const double holeDisplacement = diameter % 2 == 0 ? 10 + 3 : 10 + 2;
    const double radii = static_cast<double>(diameter) / 2.0;
    const auto dim = dimensions();
    const auto sp = spacing();

    //making phantom
    auto dBufferPtr = std::make_shared<std::vector<double>>(size(), airDensity());
    auto mBufferPtr = std::make_shared<std::vector<unsigned char>>(size(), 0);
    auto measurementPtr = std::make_shared<std::vector<std::uint8_t>>(size(), 0);
    setDensityArray(dBufferPtr);
    setMaterialIndexArray(mBufferPtr);
    setMeasurementMapArray(measurementPtr);

    //air holes indices
    std::array<std::size_t, 2> fdim = { dim[0], dim[1] };
    std::array<double, 2> fspacing = { sp[0], sp[1] };
    std::array<double, 2> fcenter = { dim[0] * sp[0] * 0.5, dim[1] * sp[1] * 0.5 };

    std::array<std::array<double, 2>, 5> hcenters;
    hcenters[0][0] = 0.0;
    hcenters[0][1] = 0.0;
    hcenters[1][0] = 0.0;
    hcenters[1][1] = radii - holeDisplacement;
    hcenters[2][0] = 0.0;
    hcenters[2][1] = -radii + holeDisplacement;
    hcenters[3][0] = radii - holeDisplacement;
    hcenters[3][1] = 0.0;
    hcenters[4][0] = -radii + holeDisplacement;
    hcenters[4][1] = 0.0;
    for (auto& cvec : hcenters)
        for (std::size_t i = 0; i < 2; ++i)
            cvec[i] += fcenter[i];

    //setting up measurement indices
    for (std::size_t k = 0; k < dim[2]; ++k) {
        const double slicePosCenter = sp[2] * (k - dim[2] * 0.5 + 0.5);
        if (std::abs(slicePosCenter) <= 50.0) // from -50 to +50 mm into center of the phantom
        {
            const std::size_t offset = k * dim[0] * dim[1];
            for (std::size_t i = 0; i < 5; ++i) {
                auto idx = circleIndices2D(fdim, fspacing, hcenters[i], holeRadii);
                std::transform(idx.begin(), idx.end(), idx.begin(), [=](std::size_t ind) { return ind + offset; }); //adding offset
                auto& indices = m_holePositions[i];
                indices.insert(indices.end(), idx.begin(), idx.end());
            }
        }
    }

    //making phantom
    auto dBuffer = dBufferPtr->data();
    auto mBuffer = mBufferPtr->data();

    for (std::size_t k = 0; k < dim[2]; ++k) {
        const std::size_t offset = k * dim[0] * dim[1];
        auto indices = circleIndices2D(fdim, fspacing, fcenter, radii);
        for (auto idx : indices) {
            dBuffer[idx + offset] = pmma.standardDensity();
            mBuffer[idx + offset] = 2;
        }
        //air holes
        for (std::size_t i = 0; i < 5; ++i) {
            auto indicesHoles = circleIndices2D(fdim, fspacing, hcenters[i], holeRadii);
            for (auto idx : indicesHoles) {
                dBuffer[idx + offset] = air.standardDensity();
                mBuffer[idx + offset] = 1;
            }
        }
    }

    //filling measurementArray
    auto measBuffer = measurementPtr->data();
    for (const auto vIdx : m_holePositions)
        for (const auto idx : vIdx)
            measBuffer[idx] = 1;
    validate();
}
constexpr double CTDIPhantom::airDensity()
{
    return AIRDENSITY;
}

const std::vector<std::size_t>& CTDIPhantom::holeIndices(CTDIPhantom::HolePosition position)
{
    if (position == HolePosition::West)
        return m_holePositions[4];
    else if (position == HolePosition::East)
        return m_holePositions[3];
    else if (position == HolePosition::North)
        return m_holePositions[2];
    else if (position == HolePosition::South)
        return m_holePositions[1];
    return m_holePositions[0];
}
}
*/