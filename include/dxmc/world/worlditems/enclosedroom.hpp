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
#include "dxmc/vectormath.hpp"
#include "dxmc/world/basicshapes/aabb.hpp"

#include <limits>
#include <optional>

namespace dxmc {

template <std::size_t NMaterialShells = 5, int Lowenergycorrection = 2>
class EnclosedRoom {
public:
    EnclosedRoom(double wallthickness = 10, const std::array<double, 6>& inner_aabb = { -1, -1, -1, 1, 1, 1 })
        : m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        m_wallThickness = std::max(std::abs(wallthickness), 0.001);
        setInnerRoomAABB(inner_aabb);
        m_density = NISTMaterials::density("Air, Dry (near sea level)");
    }

    void setWallThickness(double cm)
    {
        m_wallThickness = std::max(std::abs(cm), 0.001);
        m_outerAABB = outerAABfromInner(m_innerAABB, m_wallThickness);
    }

    void setMaterial(const Material<NMaterialShells>& material) { m_material = material; }

    void setMaterial(const Material<NMaterialShells>& material, double density)
    {
        m_material = material;
        setDensity(density);
    }

    void setDensity(double dens)
    {
        m_density = std::abs(dens);
    }

    void setInnerRoomAABB(const std::array<double, 6>& aabb)
    {
        m_innerAABB = aabb;
        m_outerAABB = outerAABfromInner(m_innerAABB, m_wallThickness);
    }

    void translate(const std::array<double, 3>& dist)
    {
        for (std::size_t i = 0; i < 3; ++i) {
            m_innerAABB[i] += dist[i];
            m_innerAABB[i + 3] += dist[i];
        }
    }

    std::array<double, 3> center() const
    {
        std::array center = {
            (m_innerAABB[0] + m_innerAABB[3]) / 2,
            (m_innerAABB[1] + m_innerAABB[4]) / 2,
            (m_innerAABB[2] + m_innerAABB[5]) / 2
        };
        return center;
    }

    const std::array<double, 6>& AABB() const
    {
        return m_outerAABB;
    }

    WorldIntersectionResult intersect(const ParticleType auto& p) const
    {
        const bool is_inside_inner = basicshape::AABB::pointInside(p.pos, m_innerAABB);
        WorldIntersectionResult intersect;
        if (is_inside_inner) {
            intersect = basicshape::AABB::intersect(p, m_innerAABB);
            intersect.rayOriginIsInsideItem = false;
        } else {
            intersect = basicshape::AABB::intersect(p, m_outerAABB);
            if (intersect.rayOriginIsInsideItem) {
                // test if we intersect inner aabb
                const auto intersect_inner = basicshape::AABB::intersect(p, m_innerAABB);
                if (intersect_inner.valid() && intersect_inner.intersection < intersect.intersection) {
                    intersect = intersect_inner;
                }
            }
        }
        return intersect;
    }

    template <typename U>
    VisualizationIntersectionResult<U> intersectVisualization(const ParticleType auto& p) const
    {
        VisualizationIntersectionResult<U> intersection = basicshape::AABB::template intersectVisualization<U>(p, m_innerAABB);
        if (intersection.valid()) {
            if (intersection.rayOriginIsInsideItem) {
                intersection.rayOriginIsInsideItem = false;
            } else {
                // we intersect inner box from outside and want to render closest walls invisible
                auto p_copy = p;
                p_copy.border_translate(intersection.intersection);
                const auto past_wall_intersection = basicshape::AABB::template intersectVisualization<U>(p_copy, m_innerAABB);
                intersection.intersection += past_wall_intersection.intersection;
                intersection.normal = past_wall_intersection.normal;
            }
        } else {
            intersection = basicshape::AABB::template intersectVisualization<U>(p, m_outerAABB);
        }
        return intersection;
    }

    const EnergyScore& energyScored(std::size_t index = 0) const
    {
        return m_energyScore;
    }

    void clearEnergyScored()
    {
        m_energyScore.clear();
    }

    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        std::array<double, 3> inner_sides;
        std::array<double, 3> outer_sides;
        for (std::size_t i = 0; i < 3; ++i) {
            inner_sides[i] = m_innerAABB[i + 3] - m_innerAABB[i];
            outer_sides[i] = m_outerAABB[i + 3] - m_outerAABB[i];
        }
        const auto inner_volume = std::reduce(inner_sides.cbegin(), inner_sides.cend(), 1.0, std::multiplies<>());
        const auto outer_volume = std::reduce(outer_sides.cbegin(), outer_sides.cend(), 1.0, std::multiplies<>());

        m_dose.addScoredEnergy(m_energyScore, outer_volume - inner_volume, m_density, calibration_factor);
    }

    const DoseScore& doseScored(std::size_t index = 0) const
    {
        return m_dose;
    }

    void clearDoseScored()
    {
        m_dose.clear();
    }

    void transport(ParticleType auto& p, RandomState& state) noexcept
    {
        bool cont = pointInside(p.pos);
        bool updateAtt = true;
        AttenuationValues att;
        double attSumInv;
        while (cont) {
            if (updateAtt) {
                att = m_material.attenuationValues(p.energy);
                attSumInv = 1 / (att.sum() * m_density);
                updateAtt = false;
            }
            const auto stepLen = -std::log(state.randomUniform()) * attSumInv; // cm
            const auto intLen = intersect(p).intersection; // this can not be nullopt

            if (stepLen < intLen) {
                // interaction happends
                p.translate(stepLen);
                const auto intRes = interactions::template interact<NMaterialShells, Lowenergycorrection>(att, p, m_material, state);
                m_energyScore.scoreEnergy(intRes.energyImparted);
                cont = intRes.particleAlive;
                updateAtt = intRes.particleEnergyChanged;
            } else {
                // transport to border
                p.border_translate(intLen);
                cont = false;
            }
        }
    }

protected:
    static std::array<double, 6> outerAABfromInner(const std::array<double, 6>& inner, double wall_thickness)
    {
        std::array<double, 6> aabb = {
            inner[0] - wall_thickness,
            inner[1] - wall_thickness,
            inner[2] - wall_thickness,
            inner[3] + wall_thickness,
            inner[4] + wall_thickness,
            inner[5] + wall_thickness
        };
        return aabb;
    }

    bool pointInside(const std::array<double, 3>& p) const noexcept
    {
        return basicshape::AABB::pointInside(p, m_outerAABB) && !basicshape::AABB::pointInside(p, m_innerAABB);
    }

private:
    double m_wallThickness = 10;
    std::array<double, 6> m_innerAABB = { -1, -1, -1, 1, 1, 1 };
    std::array<double, 6> m_outerAABB = { -1, -1, -1, 1, 1, 1 };
    double m_density = 1;
    Material<NMaterialShells> m_material;
    EnergyScore m_energyScore;
    DoseScore m_dose;
};
}