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

Copyright 2022 Erlend Andersen
*/

#pragma once

#include "dxmc/dxmcrandom.hpp"
#include "dxmc/floating.hpp"
#include "dxmc/interactions.hpp"
#include "dxmc/material/material.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/basicshapes/aabb.hpp"
#include "dxmc/world/basicshapes/cylinder.hpp"
#include "dxmc/world/worlditems/worlditembase.hpp"

#include <limits>
#include <optional>

namespace dxmc {

template <Floating T, std::size_t NMaterialShells = 5, int LOWENERGYCORRECTION = 2>
class TG195World42 final : public WorldItemBase<T> {
protected:
    static bool insideChild(const Particle<T>& p, const basicshape::cylinder::Cylinder<T>& child, const std::array<T, 6>& aabb)
    {
        if (basicshape::AABB::pointInside(p.pos, aabb)) {
            return basicshape::cylinder::pointInside(p.pos, child);
        }
        return false;
    }

public:
    TG195World42(T radius = T { 16 }, T height = T { 10 }, const std::array<T, 3>& pos = { 0, 0, 0 })
        : WorldItemBase<T>()
        , m_material(Material<T, NMaterialShells>::byNistName("Air, Dry (near sea level)").value())

    {
        m_cylinder.center = pos;
        m_cylinder.direction = { 0, 0, 1 };
        m_cylinder.radius = radius;
        m_cylinder.half_height = height / 2;

        m_materialDensity = NISTMaterials<T>::density("Air, Dry (near sea level)");

        m_centerChild.center = m_cylinder.center;
        m_centerChild.direction = m_cylinder.direction;
        m_centerChild.radius = T { 0.5 };
        m_centerChild.half_height = 5;
        m_periferyChild.center = { m_cylinder.center[0] - (m_cylinder.radius - 1), m_cylinder.center[1], m_cylinder.center[2] };
        m_periferyChild.direction = m_cylinder.direction;
        m_periferyChild.radius = T { 0.5 };
        m_periferyChild.half_height = 5;
        for (std::size_t i = 0; i < 3; ++i) {
            m_centerChild_aabb[i] = m_centerChild.center[i] - m_centerChild.radius;
            m_centerChild_aabb[i + 3] = m_centerChild.center[i] + m_centerChild.radius;
            m_periferyChild_aabb[i] = m_periferyChild.center[i] - m_periferyChild.radius;
            m_periferyChild_aabb[i + 3] = m_periferyChild.center[i] + m_periferyChild.radius;
            if (i == 2) {
                m_centerChild_aabb[i] = m_centerChild.center[i] - m_centerChild.half_height;
                m_centerChild_aabb[i + 3] = m_centerChild.center[i] + m_centerChild.half_height;
                m_periferyChild_aabb[i] = m_periferyChild.center[i] - m_periferyChild.half_height;
                m_periferyChild_aabb[i + 3] = m_periferyChild.center[i] + m_periferyChild.half_height;
            }
        }
    }

    void setMaterial(const Material<T, NMaterialShells>& material)
    {
        m_material = material;
    }

    void setMaterialDensity(T density) { m_materialDensity = density; }

    void translate(const std::array<T, 3>& dist) override
    {
        m_cylinder.center = vectormath::add(m_cylinder.center, dist);
        m_centerChild.center = vectormath::add(m_centerChild.center, dist);
        m_periferyChild.center = vectormath::add(m_periferyChild.center, dist);
        for (std::size_t i = 0; i < 3; ++i) {
            m_periferyChild_aabb[i] += dist[i];
            m_periferyChild_aabb[i + 3] += dist[i];
            m_centerChild_aabb[i] += dist[i];
            m_centerChild_aabb[i + 3] += dist[i];
        }
    }

    std::array<T, 3> center() const override
    {
        return m_cylinder.center;
    }

    std::array<T, 6> AABB() const override
    {
        std::array<T, 6> aabb {
            m_cylinder.center[0] - m_cylinder.radius,
            m_cylinder.center[1] - m_cylinder.radius,
            m_cylinder.center[2] - m_cylinder.radius,
            m_cylinder.center[0] + m_cylinder.radius,
            m_cylinder.center[1] + m_cylinder.radius,
            m_cylinder.center[2] + m_cylinder.radius
        };
        return aabb;
    }

    WorldIntersectionResult<T> intersect(const Particle<T>& p) const noexcept override
    {
        return basicshape::cylinder::intersect(p, m_cylinder);
    }

    VisualizationIntersectionResult<T, WorldItemBase<T>> intersectVisualization(const Particle<T>& p) const noexcept override
    {
        return basicshape::cylinder::template intersectVisualization<T, WorldItemBase<T>>(p, m_cylinder);
    }

    void transport(Particle<T>& p, RandomState& state) noexcept override
    {
        bool cont = basicshape::cylinder::pointInside(p.pos, m_cylinder);
        bool updateAtt = true;
        AttenuationValues<T> att;
        T attSumInv;
        while (cont) {
            if (updateAtt) {
                att = m_material.attenuationValues(p.energy);
                attSumInv = 1 / (att.sum() * m_materialDensity);
                updateAtt = false;
            }

            constexpr T alpha = T { 0.4 };

            const auto stepLen = -std::log(state.randomUniform<T>()) * attSumInv * alpha;

            const auto intLen = intersect(p);
            if (stepLen < intLen.intersection) {
                // interaction happends
                p.translate(stepLen);
                if (insideChild(p, m_centerChild, m_centerChild_aabb)) {
                    const auto intRes = interactions::template interactForced<T, NMaterialShells, LOWENERGYCORRECTION>(alpha, att, p, m_material, state);
                    m_centerChild_energyScored.scoreEnergy(intRes.energyImparted);
                    updateAtt = intRes.particleEnergyChanged;
                    cont = intRes.particleAlive;
                } else if (insideChild(p, m_periferyChild, m_periferyChild_aabb)) {
                    const auto intRes = interactions::template interactForced<T, NMaterialShells, LOWENERGYCORRECTION>(alpha, att, p, m_material, state);
                    m_periferyChild_energyScored.scoreEnergy(intRes.energyImparted);
                    updateAtt = intRes.particleEnergyChanged;
                    cont = intRes.particleAlive;
                } else if (state.randomUniform<T>() < alpha) {
                    const auto intRes = interactions::template interact<T, NMaterialShells, LOWENERGYCORRECTION>(att, p, m_material, state);
                    updateAtt = intRes.particleEnergyChanged;
                    cont = intRes.particleAlive;
                }

            } else {
                // transport to border
                p.border_translate(intLen.intersection);
                cont = false;
            }
        }
    }

    void clearEnergyScored() override
    {
        m_energyScored.clear();
        m_centerChild_energyScored.clear();
        m_periferyChild_energyScored.clear();
    }
    const EnergyScore<T>& energyScoredCenterCylinder() const
    {
        return m_centerChild_energyScored;
    }
    const EnergyScore<T>& energyScoredPeriferyCylinder() const
    {
        return m_periferyChild_energyScored;
    }
    const EnergyScore<T>& energyScored(std::size_t index = 0) const override
    {
        if (index == 0)
            return m_centerChild_energyScored;
        else if (index == 1)
            return m_periferyChild_energyScored;
        return m_energyScored;
    }
    void addEnergyScoredToDoseScore(T calibration_factor = 1)
    {
        const T c_vol = m_centerChild.volume();
        const T p_vol = m_periferyChild.volume();
        const T vol = m_cylinder.volume() - c_vol - p_vol;
        m_doseScored.addScoredEnergy(m_energyScored, vol, m_materialDensity, calibration_factor);
        m_centerChild_doseScored.addScoredEnergy(m_centerChild_energyScored, c_vol, m_materialDensity, calibration_factor);
        m_periferyChild_doseScored.addScoredEnergy(m_periferyChild_energyScored, p_vol, m_materialDensity, calibration_factor);
    };
    const DoseScore<T>& doseScored(std::size_t index = 0) const
    {
        if (index == 0)
            return m_centerChild_doseScored;
        else if (index == 1)
            return m_periferyChild_doseScored;
        return m_doseScored;
    }
    void clearDoseScored()
    {
        m_doseScored.clear();
        m_centerChild_doseScored.clear();
        m_periferyChild_doseScored.clear();
    }

private:
    basicshape::cylinder::Cylinder<T> m_cylinder;
    T m_materialDensity = 1;
    Material<T, NMaterialShells> m_material;
    EnergyScore<T> m_energyScored;
    EnergyScore<T> m_periferyChild_energyScored;
    EnergyScore<T> m_centerChild_energyScored;
    DoseScore<T> m_doseScored;
    DoseScore<T> m_periferyChild_doseScored;
    DoseScore<T> m_centerChild_doseScored;
    basicshape::cylinder::Cylinder<T> m_centerChild;
    basicshape::cylinder::Cylinder<T> m_periferyChild;
    std::array<T, 6> m_centerChild_aabb = { 0, 0, 0, 0, 0, 0 };
    std::array<T, 6> m_periferyChild_aabb = { 0, 0, 0, 0, 0, 0 };
};
}
