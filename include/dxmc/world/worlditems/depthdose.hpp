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
#include "dxmc/interactions.hpp"
#include "dxmc/material/material.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/particletracker.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/basicshapes/aabb.hpp"
#include "dxmc/world/basicshapes/cylinder.hpp"

#include <limits>
#include <optional>

namespace dxmc {

template <std::size_t NMaterialShells = 5, int LOWENERGYCORRECTION = 2, bool FORCEINTERACTIONS = false>
class DepthDose {
public:
    DepthDose(double radius = 16, double height = 10, std::size_t resolution = 100, const std::array<double, 3>& pos = { 0, 0, 0 })
        : m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        m_cylinder.center = pos;
        m_cylinder.direction = { 0, 0, 1 };
        m_cylinder.radius = radius;
        m_cylinder.half_height = height / 2;

        m_materialDensity = NISTMaterials::density("Air, Dry (near sea level)");
        m_energyScored.resize(resolution);
        m_dose.resize(resolution);
    }

    double length() const { return m_cylinder.half_height * 2; }
    double radius() const { return m_cylinder.radius; }

    std::size_t resolution() const
    {
        return m_energyScored.size();
    }

    void setMaterial(const Material<NMaterialShells>& material)
    {
        m_material = material;
    }

    void setMaterialDensity(double density)
    {
        m_materialDensity = std::abs(density);
    }

    void setMaterial(const Material<NMaterialShells>& material, double density)
    {
        setMaterial(material);
        setMaterialDensity(density);
    }

    bool setNistMaterial(const std::string& nist_name)
    {
        const auto mat = Material<NMaterialShells>::byNistName(nist_name);
        if (mat) {
            m_material = mat.value();
            m_materialDensity = NISTMaterials::density(nist_name);
            return true;
        }
        return false;
    }

    const Material<NMaterialShells>& material() const
    {
        return m_material;
    }

    double density() const { return m_materialDensity; }

    void translate(const std::array<double, 3>& dist)
    {
        m_cylinder.center = vectormath::add(m_cylinder.center, dist);
    }

    const std::array<double, 3>& center() const
    {
        return m_cylinder.center;
    }

    std::array<double, 6> AABB() const
    {
        std::array aabb {
            m_cylinder.center[0] - m_cylinder.radius,
            m_cylinder.center[1] - m_cylinder.radius,
            m_cylinder.center[2] - m_cylinder.half_height,
            m_cylinder.center[0] + m_cylinder.radius,
            m_cylinder.center[1] + m_cylinder.radius,
            m_cylinder.center[2] + m_cylinder.half_height
        };
        return aabb;
    }

    WorldIntersectionResult intersect(const ParticleType auto& p) const
    {
        return basicshape::cylinder::intersect(p, m_cylinder);
    }

    template <typename U>
    VisualizationIntersectionResult<U> intersectVisualization(const ParticleType auto& p) const
    {
        auto inter = basicshape::cylinder::template intersectVisualization<U>(p, m_cylinder);
        if (inter.valid()) {
            auto p_int = p;
            p_int.translate(inter.intersection);
            const auto dose_ind_f = (p_int.pos[2] - (m_cylinder.center[2] - m_cylinder.half_height)) * m_energyScored.size() / (m_cylinder.half_height * 2);
            const auto ind = std::clamp(static_cast<std::size_t>(dose_ind_f), std::size_t { 0 }, m_energyScored.size() - 1);
            inter.value = m_dose[ind].dose();
        }
        return inter;
    }

    template <ParticleType P>
    void transport(P& p, RandomState& state) noexcept
    {
        if constexpr (std::is_same<P, ParticleTrack>::value) {
            m_tracker.registerParticle(p);
        }
        if constexpr (FORCEINTERACTIONS)
            transportForced(p, state);
        else
            transportRandom(p, state);
    }

    const std::vector<std::pair<double, EnergyScore>> depthEnergyScored() const
    {
        std::vector<std::pair<double, EnergyScore>> depth;
        depth.reserve(m_energyScored.size());

        const auto step = (2 * m_cylinder.half_height) / m_energyScored.size();
        const auto start = m_cylinder.center[2] - m_cylinder.half_height + step / 2;
        for (std::size_t i = 0; i < m_energyScored.size(); ++i) {
            depth.push_back(std::make_pair(start + step * i, m_energyScored[i]));
        }
        return depth;
    }

    const EnergyScore& energyScored(std::size_t index = 0) const
    {
        return m_energyScored[index];
    }

    void clearEnergyScored()
    {
        for (auto& d : m_energyScored) {
            d.clear();
        }
    }

    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        const auto totalVolume = m_cylinder.volume();
        const auto partVolume = totalVolume / m_dose.size();
        for (std::size_t i = 0; i < m_dose.size(); ++i) {
            m_dose[i].addScoredEnergy(m_energyScored[i], partVolume, m_materialDensity, calibration_factor);
        }
    }

    const DoseScore& doseScored(std::size_t index = 0) const
    {
        return m_dose[index];
    }

    void clearDoseScored()
    {
        for (auto& d : m_dose) {
            d.clear();
        }
    }

protected:
    void transportRandom(ParticleType auto& p, RandomState& state)
    {
        bool cont = basicshape::cylinder::pointInside(p.pos, m_cylinder);
        bool updateAtt = false;
        auto att = m_material.attenuationValues(p.energy);
        auto attSumInv = 1 / (att.sum() * m_materialDensity);
        while (cont) {
            if (updateAtt) {
                att = m_material.attenuationValues(p.energy);
                attSumInv = 1 / (att.sum() * m_materialDensity);
                updateAtt = false;
            }
            const auto stepLen = -std::log(state.randomUniform()) * attSumInv; // cm
            const auto intLen = intersect(p).intersection; // this can not be nullopt

            if (stepLen < intLen) {
                // interaction happends
                p.translate(stepLen);
                const auto intRes = interactions::template interact<NMaterialShells, LOWENERGYCORRECTION>(att, p, m_material, state);
                updateAtt = intRes.particleEnergyChanged;
                cont = intRes.particleAlive;

                const auto dose_ind_f = (p.pos[2] - (m_cylinder.center[2] - m_cylinder.half_height)) * m_energyScored.size() / (m_cylinder.half_height * 2);
                const auto ind = std::clamp(static_cast<std::size_t>(dose_ind_f), std::size_t { 0 }, m_energyScored.size() - 1);
                m_energyScored[ind].scoreEnergy(intRes.energyImparted);
            } else {
                // transport to border
                p.border_translate(intLen);
                cont = false;
            }
        }
    }

    void transportForced(ParticleType auto& p, RandomState& state) noexcept
    {
        bool cont = basicshape::cylinder::pointInside(p.pos, m_cylinder);
        while (cont) {
            const auto intLen = intersect(p).intersection; // this must be valid
            const auto intRes = interactions::template interactForced<NMaterialShells, LOWENERGYCORRECTION>(intLen, m_materialDensity, p, m_material, state);

            const auto dose_ind_f = (p.pos[2] - (m_cylinder.center[2] - m_cylinder.half_height)) * m_energyScored.size() / (m_cylinder.half_height * 2);
            const auto ind = std::clamp(static_cast<std::size_t>(dose_ind_f), std::size_t { 0 }, m_energyScored.size() - 1);
            m_energyScored[ind].scoreEnergy(intRes.energyImparted);
            cont = intRes.particleAlive && basicshape::cylinder::pointInside(p.pos, m_cylinder);
        }
    }

private:
    dxmc::basicshape::cylinder::Cylinder m_cylinder;
    double m_materialDensity = 1;
    Material<NMaterialShells> m_material;
    ParticleTracker m_tracker;
    std::vector<EnergyScore> m_energyScored;
    std::vector<DoseScore> m_dose;
};

}
