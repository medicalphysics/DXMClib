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
class WorldBox final : public WorldItemBase<T> {
public:
    WorldBox(const std::array<T, 6>& aabb = { -1, -1, -1, 1, 1, 1 })
        : WorldItemBase<T>()
        , m_aabb(aabb)
        , m_material(Material<T, NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        m_materialDensity = NISTMaterials<T>::density("Air, Dry (near sea level)");
    }

    WorldBox(T aabb_size, std::array<T, 3> pos = { 0, 0, 0 })
        : WorldItemBase<T>()
        , m_material(Material<T, NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] = -std::abs(aabb_size) + pos[i];
            m_aabb[i + 3] = std::abs(aabb_size) + pos[i];
        }
        m_materialDensity = NISTMaterials<T>::density("Air, Dry (near sea level)");
    }

    void setMaterial(const Material<T, NMaterialShells>& material)
    {
        m_material = material;
    }
    void setMaterial(const Material<T, NMaterialShells>& material, T density)
    {
        m_material = material;
        m_materialDensity = std::abs(density);
    }
    void setMaterialDensity(T density) { m_materialDensity = std::abs(density); }

    bool setNistMaterial(const std::string& nist_name)
    {
        const auto mat = Material<T, NMaterialShells>::byNistName(nist_name);
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

    void clearEnergyScored() override
    {
        m_energyScored.clear();
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

    VisualizationIntersectionResult<T, WorldItemBase<T>> intersectVisualization(const Particle<T>& p) const noexcept override
    {
        return basicshape::AABB::template intersectVisualization<T, WorldItemBase<T>>(p, m_aabb);
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
                m_energyScored.scoreEnergy(intRes.energyImparted);
                cont = intRes.particleAlive;
                updateAtt = intRes.particleEnergyChanged;

            } else {
                // transport to border
                p.border_translate(intLen);
                cont = false;
            }
        }
    }

    const EnergyScore<T>& energyScored(std::size_t index = 0) const override
    {
        return m_energyScored;
    }

protected:
private:
    std::array<T, 6> m_aabb;
    Material<T, NMaterialShells> m_material;
    T m_materialDensity = 1;
    EnergyScore<T> m_energyScored;
};

}
