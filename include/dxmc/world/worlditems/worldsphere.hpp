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
#include "dxmc/world/basicshapes/sphere.hpp"
#include "dxmc/world/worlditems/worlditembase.hpp"

#include <limits>
#include <optional>

namespace dxmc {

template <Floating T, std::size_t NMaterialShells = 5, int LOWENERGYCORRECTION = 2, bool FORCEINTERACTIONS = false>
class WorldSphere final : public WorldItemBase<T> {
public:
    WorldSphere(T radius = T { 16 }, const std::array<T, 3>& pos = { 0, 0, 0 })
        : WorldItemBase<T>()
        , m_radius(std::abs(radius))
        , m_center(pos)
        , m_material(Material<T, NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        m_materialDensity = NISTMaterials<T>::density("Air, Dry (near sea level)");
    }

    void setRadius(T r)
    {
        m_radius = std::abs(r);
    }

    T radius() const { return m_radius; }

    void setMaterial(const Material<T, NMaterialShells>& material)
    {
        m_material = material;
    }

    void setMaterial(const Material<T, NMaterialShells>& material, T density)
    {
        m_material = material;
        setMaterialDensity(density);
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

    void translate(const std::array<T, 3>& dist) override
    {
        m_center = vectormath::add(m_center, dist);
    }

    std::array<T, 3> center() const override
    {
        return m_center;
    }

    std::array<T, 6> AABB() const override
    {
        std::array<T, 6> aabb {
            m_center[0] - m_radius,
            m_center[1] - m_radius,
            m_center[2] - m_radius,
            m_center[0] + m_radius,
            m_center[1] + m_radius,
            m_center[2] + m_radius
        };
        return aabb;
    }

    WorldIntersectionResult<T> intersect(const Particle<T>& p) const noexcept override
    {
        return basicshape::sphere::intersect(p, m_center, m_radius);
    }

    VisualizationIntersectionResult<T, WorldItemBase<T>> intersectVisualization(const Particle<T>& p) const noexcept override
    {
        auto inter = basicshape::sphere::template intersectVisualization<T, WorldItemBase<T>>(p, m_center, m_radius);
        if (inter.valid())
            inter.value = m_dose.dose();
        return inter;
    }

    void transport(Particle<T>& p, RandomState& state) noexcept override
    {
        if constexpr (FORCEINTERACTIONS)
            transportForced(p, state);
        else
            transportRandom(p, state);
    }

    const EnergyScore<T>& energyScored(std::size_t index = 0) const override
    {
        return m_energyScored;
    }

    void clearEnergyScored() override
    {
        m_energyScored.clear();
    }

    void addEnergyScoredToDoseScore(T calibration_factor = 1) override
    {
        const auto volume = (4 * std::numbers::pi_v<T> * m_radius * m_radius * m_radius) / 3;
        m_dose.addScoredEnergy(m_energyScored, volume, m_materialDensity, calibration_factor);
    }

    const DoseScore<T>& doseScored(std::size_t index = 0) const override
    {
        return m_dose;
    }

    void clearDoseScored() override
    {
        m_dose.clear();
    }

protected:
    void transportRandom(Particle<T>& p, RandomState& state) noexcept
    {
        bool cont = basicshape::sphere::pointInside(p.pos, m_center, m_radius);
        bool updateAtt = false;
        auto att = m_material.attenuationValues(p.energy);
        auto attSumInv = 1 / (att.sum() * m_materialDensity);
        while (cont) {
            if (updateAtt) {
                att = m_material.attenuationValues(p.energy);
                attSumInv = 1 / (att.sum() * m_materialDensity);
                updateAtt = false;
            }
            const auto stepLen = -std::log(state.randomUniform<T>()) * attSumInv; // cm
            const auto intLen = intersect(p).intersection; // this must be valid

            if (stepLen < intLen) {
                // interaction happends
                p.translate(stepLen);
                const auto intRes = interactions::template interact<T, NMaterialShells, LOWENERGYCORRECTION>(att, p, m_material, state);
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

    void transportForced(Particle<T>& p, RandomState& state) noexcept
    {
        bool cont = basicshape::sphere::pointInside(p.pos, m_center, m_radius);
        bool updateAtt = false;
        auto att = m_material.attenuationValues(p.energy);
        auto attSumInv = 1 / (att.sum() * m_materialDensity);
        while (cont) {
            if (updateAtt) {
                att = m_material.attenuationValues(p.energy);
                attSumInv = 1 / (att.sum() * m_materialDensity);
                updateAtt = false;
            }
            const auto stepLen = -std::log(state.randomUniform<T>()) * attSumInv; // cm
            const auto intLen = intersect(p).intersection; // this must be valid

            if (stepLen < intLen) {
                // interaction happends
                p.translate(stepLen);
                const auto photoProb = att.photoelectric / att.sum();
                if (state.randomUniform<T>() < photoProb) {
                    m_energyScored.scoreEnergy(p.energy * p.weight);
                    cont = false;
                } else {
                    const auto intRes = interactions::template interactScatter<T, NMaterialShells, LOWENERGYCORRECTION>(att, p, m_material, state);
                    m_energyScored.scoreEnergy(intRes.energyImparted);
                    cont = intRes.particleAlive;
                    updateAtt = intRes.particleEnergyChanged;
                }
            } else {
                // transport to border
                p.border_translate(intLen);
                cont = false;
            }
        }
    }

private:
    T m_radius = 0;
    std::array<T, 3> m_center;
    Material<T, NMaterialShells> m_material;
    T m_materialDensity = 1;
    EnergyScore<T> m_energyScored;
    DoseScore<T> m_dose;
};
}
