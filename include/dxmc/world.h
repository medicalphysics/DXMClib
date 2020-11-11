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

#pragma once // include guard

#include "dxmc/floating.h"
#include "dxmc/material.h"
#include "dxmc/tube.h"
#include "dxmc/vectormath.h"

#include <algorithm>
#include <array>
#include <execution>
#include <map>
#include <memory>
#include <ranges>
#include <string>
#include <vector>

namespace dxmc {

template <Floating T = double>
class World {
private:
    std::array<T, 3> m_spacing = { 1.0, 1.0, 1.0 };
    std::array<T, 3> m_origin = { 0.0, 0.0, 0.0 };
    std::array<T, 6> m_directionCosines = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0 };
    std::array<std::size_t, 3> m_dimensions = { 0L, 0L, 0L };
    std ::array<T, 6> m_matrixExtent { 0, 0, 0, 0, 0, 0 };

    // m_density and m_materialIndex can outlive member variables of this class, i.e shared pointers
    std::shared_ptr<std::vector<T>> m_density = nullptr;
    std::shared_ptr<std::vector<std::uint8_t>> m_materialIndex = nullptr;
    std::shared_ptr<std::vector<std::uint8_t>> m_measurementMap = nullptr;
    std::vector<Material> m_materialMap;
    bool m_valid = false;

protected:
    void updateMatrixExtent(void)
    {
        for (std::size_t i = 0; i < 3; i++) {
            const T halfDist = (m_dimensions[i] * m_spacing[i]) * T { 0.5 };
            const T lower = m_origin[i] - halfDist;
            const T upper = m_origin[i] + halfDist;
            // Preventing rounding errors when using extent to determine if a
            // particle is inside the world
            m_matrixExtent[i * 2] = std::nextafter(lower, upper);
            m_matrixExtent[i * 2 + 1] = std::nextafter(upper, lower);
        }
    }

    bool validate();

public:
    World()
    {
        updateMatrixExtent();
    }
    void setDimensions(const std::array<std::size_t, 3>& dimensions)
    {
        m_valid = false;
        m_dimensions = dimensions;
        updateMatrixExtent();
    }
    void setSpacing(const std::array<T, 3>& spacing)
    {
        m_valid = false;
        m_spacing = spacing;
        updateMatrixExtent();
    }
    void setOrigin(const std::array<T, 3>& origin)
    {
        m_valid = false;
        m_origin = origin;
        updateMatrixExtent();
    }

    void setDirectionCosines(const std::array<T, 6>& cosines)
    {
        m_valid = false;
        m_directionCosines = cosines;
        vectormath::normalize(m_directionCosines.data());
        vectormath::normalize(&m_directionCosines[3]);
    }

    std::size_t size(void) const { return m_dimensions[0] * m_dimensions[1] * m_dimensions[2]; }
    const std::array<std::size_t, 3>& dimensions(void) const { return m_dimensions; }
    const std::array<T, 3>& spacing(void) const { return m_spacing; }
    const std::array<T, 3>& origin(void) const { return m_origin; }
    const std::array<T, 6>& directionCosines(void) const { return m_directionCosines; }
    std::array<T, 3> depthDirection(void) const
    {
        std::array<T, 3> dir;
        vectormath::cross(m_directionCosines.data(), dir.data());
        return dir;
    }
    void setDensityArray(std::shared_ptr<std::vector<T>> densityArray)
    {
        m_valid = false;
        m_density = densityArray;
    }
    std::shared_ptr<std::vector<T>> densityArray(void) { return m_density; }
    const std::shared_ptr<std::vector<T>> densityArray(void) const { return m_density; }

    void setMaterialIndexArray(std::shared_ptr<std::vector<std::uint8_t>> materialIndex)
    {
        m_valid = false;
        m_materialIndex = materialIndex;
    }
    std::shared_ptr<std::vector<std::uint8_t>> materialIndexArray(void) { return m_materialIndex; }
    const std::shared_ptr<std::vector<std::uint8_t>> materialIndexArray(void) const { return m_materialIndex; }

    void setMeasurementMapArray(std::shared_ptr<std::vector<std::uint8_t>> measurementMap)
    {
        m_valid = false;
        m_measurementMap = measurementMap;
    }
    std::shared_ptr<std::vector<std::uint8_t>> measurementMapArray(void) { return m_measurementMap; }
    const std::shared_ptr<std::vector<std::uint8_t>> measurementMapArray(void) const { return m_measurementMap; }

    const std::vector<Material>& materialMap(void) const { return m_materialMap; }
    bool addMaterialToMap(const Material& material)
    {
        m_valid = false;
        if (material.isValid()) {
            m_materialMap.push_back(material);
            return true;
        }
        return false;
    }
    bool addMaterialToMap(Material&& material)
    {
        m_valid = false;
        if (material.isValid()) {
            m_materialMap.push_back(material);
            return true;
        }
        return false;
    }

    void clearMaterialMap(void)
    {
        m_valid = false;
        m_materialMap.clear();
    }
    const std::array<T, 6>& matrixExtent(void) const { return m_matrixExtent; }

    void makeValid(void)
    {
        if (!m_valid) {
            m_valid = validate();
        }
    }

    bool isValid(void) const { return m_valid; }
    [[nodiscard]] bool isValid(void)
    {
        makeValid();
        return m_valid;
    }
};

template <Floating T>
bool World<T>::validate()
{
    if ((m_dimensions[0] * m_dimensions[1] * m_dimensions[2]) <= 0) {
        return false;
    }

    if ((m_spacing[0] * m_spacing[1] * m_spacing[2]) <= T { 0 }) {
        return false;
    }

    if (!m_density || !m_materialIndex) {
        return false;
    }
    if (!m_measurementMap) {
        m_measurementMap = std::make_shared<std::vector<std::uint8_t>>(size(), 0);
    }

    const auto elements = size();
    if (elements == 0) {
        return false;
    }

    if (m_density->size() != elements || m_materialIndex->size() != elements || m_measurementMap->size() != elements) {
        return false;
    }

    //testValidMaterials
    for (auto const& mat : m_materialMap) {
        const bool validMaterial = mat.isValid();
        if (!validMaterial) {
            return false;
        }
    }

    auto depthDirectionCosine = depthDirection();
    const T testOrtogonality = vectormath::dot(depthDirectionCosine.data(), m_directionCosines.data()) + vectormath::dot(depthDirectionCosine.data(), &m_directionCosines[3]) + vectormath::dot(&m_directionCosines[0], &m_directionCosines[3]);

    if (std::abs(testOrtogonality) > 0.001) {
        return false;
    }

    auto [min, max] = std::minmax_element(std::execution::par_unseq, m_materialIndex->cbegin(), m_materialIndex->cend());
    const auto n = m_materialMap.size();
    return (*min >= 0) && (*min < n) && (*max >= 0) && (*max < n);
}

template <Floating T = double>
class CTDIPhantom final : public World<T> {

private:
    std::array<std::vector<std::size_t>, 5> m_holePositions;

protected:
    static std::vector<std::size_t> circleIndices2D(const std::array<std::size_t, 2>& dim, const std::array<T, 2>& spacing, const std::array<T, 2>& center, T radius);

public:
    enum class HolePosition { Center,
        West,
        East,
        South,
        North };
    CTDIPhantom(std::size_t size = 320); // size im mm
    const std::vector<std::size_t>& holeIndices(CTDIPhantom<T>::HolePosition position)
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
    static constexpr T airDensity() { return T { 0.001205 }; } // g/cm3
    static constexpr std::uint64_t ctdiMinHistories() { return 100E6; }
};

template <Floating T>
CTDIPhantom<T>::CTDIPhantom(std::size_t diameter)
    : World<T>()
{
    const std::array<T, 6> directionCosines { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0 };
    const std::array<T, 3> sp { 1, 1, 2.5 };
    const std::array<T, 3> origin { 0, 0, 0 };
    this->setDirectionCosines(directionCosines);
    this->setSpacing(sp);
    this->setOrigin(origin);

    if (diameter % 2 == 0) { // make sure odd dimensions in xy
        std::array<std::size_t, 3> dim { diameter + 3, diameter + 3, 60 };
        this->setDimensions(dim);
    } else {
        std::array<std::size_t, 3> dim { diameter + 2, diameter + 2, 60 };
        this->setDimensions(dim);
    }
    Material air("Air, Dry (near sea level)");
    Material pmma("Polymethyl Methacralate (Lucite, Perspex)");
    this->addMaterialToMap(air);
    this->addMaterialToMap(air); // we want two to differentiate between air around phantom and air inside rods
    this->addMaterialToMap(pmma);

    constexpr T holeDiameter { 13.1 };
    constexpr T holeRadii { holeDiameter / 2.0 };
    const T holeDisplacement = diameter % 2 == 0 ? T { 13 } : T { 13 };
    const T radii = static_cast<T>(diameter) / T { 2 };
    const auto& dim = this->dimensions();

    //making phantom
    auto dBufferPtr = std::make_shared<std::vector<T>>(this->size(), airDensity());
    auto mBufferPtr = std::make_shared<std::vector<std::uint8_t>>(this->size(), 0);
    auto measurementPtr = std::make_shared<std::vector<std::uint8_t>>(this->size(), 0);
    this->setDensityArray(dBufferPtr);
    this->setMaterialIndexArray(mBufferPtr);
    this->setMeasurementMapArray(measurementPtr);

    //air holes indices
    std::array<std::size_t, 2> fdim = { dim[0], dim[1] };
    std::array<T, 2> fspacing = { sp[0], sp[1] };
    std::array<T, 2> fcenter = { dim[0] * sp[0] * T { 0.5 }, dim[1] * sp[1] * T { 0.5 } };

    std::array<std::array<T, 2>, 5> hcenters;
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
        const T slicePosCenter = sp[2] * (k - dim[2] * 0.5 + 0.5);
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
    for (const auto& vIdx : m_holePositions)
        for (const auto idx : vIdx)
            measBuffer[idx] = 1;
    this->makeValid();
}
template <Floating T>
std::vector<std::size_t> CTDIPhantom<T>::circleIndices2D(const std::array<std::size_t, 2>& dim, const std::array<T, 2>& spacing, const std::array<T, 2>& center, T radius)
{
    std::vector<std::size_t> indices;
    const auto r_int = static_cast<std::size_t>(std::ceil(std::abs(radius)));
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

    const auto r2 = radius * radius;
    for (std::size_t i = limits[0]; i < limits[1]; ++i) {
        const auto posX = center[0] - i * spacing[0] + spacing[0] / 2;
        for (std::size_t j = limits[2]; j < limits[3]; ++j) {
            const auto posY = center[1] - j * spacing[1] + spacing[1] / 2;
            if ((posX * posX + posY * posY) <= r2)
                indices.push_back(i + j * dim[0]);
        }
    }
    return indices;
}
}