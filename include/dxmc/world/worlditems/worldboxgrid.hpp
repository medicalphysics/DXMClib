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

#pragma once

#include "dxmc/dxmcrandom.hpp"
#include "dxmc/floating.hpp"
#include "dxmc/interactions.hpp"
#include "dxmc/material/material.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/basicshapes/aabb.hpp"
#include "dxmc/world/worlditems/worlditembase.hpp"

#include <limits>
#include <optional>

namespace dxmc {

template <Floating T, std::size_t NMaterialShells = 5, int Lowenergycorrection = 2>
class WorldBoxGrid final : public WorldItemBase<T> {
public:
    WorldBoxGrid(const std::array<T, 6>& aabb = { -1, -1, -1, 1, 1, 1 })
        : WorldItemBase<T>()
        , m_aabb(aabb)
        , m_material(Material2<T, NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        m_materialDensity = NISTMaterials<T>::density("Air, Dry (near sea level)");
        setVoxelDimensions({ 1, 1, 1 });
    }

    WorldBoxGrid(T aabb_size, std::array<T, 3> pos = { 0, 0, 0 })
        : WorldItemBase<T>()
        , m_material(Material2<T, NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] = -std::abs(aabb_size) + pos[i];
            m_aabb[i + 3] = std::abs(aabb_size) + pos[i];
        }
        m_materialDensity = NISTMaterials<T>::density("Air, Dry (near sea level)");
        setVoxelDimensions({ 1, 1, 1 });
    }

    void setMaterial(const Material2<T, NMaterialShells>& material)
    {
        m_material = material;
    }
    void setVoxelDimensions(const std::array<std::size_t, 3>& dim)
    {
        m_voxelDim = dim;
        for (std::size_t i = 0; i < 3; ++i) {
            m_voxelSize[i] = (m_aabb[i + 3] - m_aabb[i]) / m_voxelDim[i];
        }
        const auto ndim = m_voxelDim[0] * m_voxelDim[1] * m_voxelDim[2];
        m_dose.resize(ndim);
    }
    std::size_t totalNumberOfVoxels() const { return m_voxelDim[0] * m_voxelDim[1] * m_voxelDim[2]; }

    const std::array<std::size_t, 3>& voxelDimensions() const { return m_voxelDim; }

    const std::array<T, 3>& voxelSpacing() const { return m_voxelSize; }

    std::size_t gridIndex(const std::array<T, 3>& pos) const noexcept
    {
        const auto x = static_cast<std::size_t>((pos[0] - m_aabb[0]) / m_voxelSize[0]);
        const auto y = static_cast<std::size_t>((pos[1] - m_aabb[1]) / m_voxelSize[1]);
        const auto z = static_cast<std::size_t>((pos[2] - m_aabb[2]) / m_voxelSize[2]);
        if (x < m_voxelDim[0] && y < m_voxelDim[1] && z < m_voxelDim[2]) {
            return x + y * m_voxelDim[0] + z * m_voxelDim[0] * m_voxelDim[1];
        }
        return m_voxelDim[0] * m_voxelDim[1] * m_voxelDim[2];
    }

    void setMaterialDensity(T density) { m_materialDensity = density; }

    bool setNistMaterial(const std::string& nist_name)
    {
        const auto mat = Material2<T, NMaterialShells>::byNistName(nist_name);
        if (mat) {
            m_material = mat.value();
            m_materialDensity = NISTMaterials<T>::density(nist_name);
            return true;
        }
        return false;
    }

    void translate(const std::array<T, 3>& dist) noexcept override
    {
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] += dist[i];
            m_aabb[i + 3] += dist[i];
        }
    }

    void clearDose() override
    {
        for (auto& d : m_dose) {
            d.clear();
        }
    }

    std::array<T, 3> center() const noexcept override
    {
        std::array<T, 3> c {
            (m_aabb[0] + m_aabb[3]) * T { 0.5 },
            (m_aabb[1] + m_aabb[4]) * T { 0.5 },
            (m_aabb[2] + m_aabb[5]) * T { 0.5 },
        };
        return c;
    }

    std::array<T, 6> AABB() const noexcept override
    {
        return m_aabb;
    }

    WorldIntersectionResult<T> intersect(const Particle<T>& p) const noexcept override
    {
        return basicshape::AABB::intersect(p, m_aabb);
    }

    void transport(Particle<T>& p, RandomState& state) noexcept override
    {
        bool cont = basicshape::AABB::pointInside(p.pos, m_aabb);
        bool updateAtt = true;
        AttenuationValues<T> att;
        T attSumInv;
        while (cont) {
            if (updateAtt) {
                att = m_material.attenuationValues(p.energy);
                attSumInv = 1 / (att.sum() * m_materialDensity);
                updateAtt = false;
            }
            const auto stepLen = -std::log(state.randomUniform<T>()) * attSumInv; // cm
            const auto intLen = intersect(p).intersection; // this can not be nullopt

            if (stepLen < intLen) {
                // interaction happends
                p.translate(stepLen);
                const auto intRes = interactions::template interact<T, NMaterialShells, Lowenergycorrection>(att, p, m_material, state);
                const auto doseIdx = gridIndex(p.pos);
                m_dose[doseIdx].scoreEnergy(intRes.energyImparted);
                cont = intRes.particleAlive;
                updateAtt = intRes.particleEnergyChanged;

            } else {
                // transport to border
                p.border_translate(intLen);
                cont = false;
            }
        }
    }

    const DoseScore<T>& dose(std::size_t index = 0) const override
    {
        return m_dose.at(index);
    }

protected:
private:
    std::array<T, 6> m_aabb;
    Material2<T, NMaterialShells> m_material;
    std::vector<DoseScore<T>> m_dose;
    T m_materialDensity = 1;
    std::array<T, 3> m_voxelSize = { 1, 1, 1 };
    std::array<std::size_t, 3> m_voxelDim = { 1, 1, 1 };
};

}
