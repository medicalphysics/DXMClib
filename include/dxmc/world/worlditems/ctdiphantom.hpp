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
#include "dxmc/world/basicshapes/cylinder.hpp"
#include "dxmc/world/statickdtree.hpp"
#include "dxmc/world/worlditems/worlditembase.hpp"

#include <array>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

namespace dxmc {

template <Floating T, int NMaterialShells = 5>
class CTDIPhantom final : public WorldItemBase<T> {
public:
    CTDIPhantom(T radius = T { 16 }, T height = T { 15 }, const std::array<T, 3>& pos = { 0, 0, 0 })
        : WorldItemBase<T>()
        , m_radius(radius)
        , m_half_height(height * T { 0.5 })
        , m_center(pos)
        , m_pmma(Material2<T, NMaterialShells>::byNistName("Polymethyl Methacralate (Lucite, Perspex)").value())
        , m_air(Material2<T, NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        m_pmma_density = NISTMaterials<T>::density("Polymethyl Methacralate (Lucite, Perspex)");
        m_air_density = NISTMaterials<T>::density("Air, Dry (near sea level)");

        std::vector<CTDIAirHole> holes;
        holes.reserve(5);
        std::array<T, 10> positions = { 0, 0, 0, 1, 1, 0, 0, -1, -1, 0 };
        std::uint8_t index = 0;
        for (std::size_t i = 0; i < 10; i = i + 2) {
            const auto x = positions[i];
            const auto y = positions[i + 1];
            CTDIAirHole hole { .pos = m_center, .radii = holeRadii(), .half_height = m_half_height, .index = ++index };
            hole.pos[0] += x * (m_radius - holeEdgeDistance());
            hole.pos[1] += y * (m_radius - holeEdgeDistance());
            holes.push_back(hole);
        }
        m_kdtree = StaticKDTree<3, T, CTDIAirHole>(holes);
    }

    void translate(const std::array<T, 3>& dist) noexcept override
    {
        m_center = vectormath::add(m_center, dist);
        m_kdtree.translate(dist);
    }
    std::array<T, 3> center() const noexcept override
    {
        return m_center;
    }

    const DoseScore<T>& dose(std::size_t index = 0) const override
    {
        return m_dose[index];
    }

    std::array<T, 6> AABB() const noexcept override
    {
        std::array aabb {
            m_center[0] - m_radius,
            m_center[1] - m_radius,
            m_center[2] - m_half_height,
            m_center[0] + m_radius,
            m_center[1] + m_radius,
            m_center[2] + m_half_height
        };
        return aabb;
    }

    WorldIntersectionResult<T> intersect(const Particle<T>& p) const noexcept override
    {
        return basicshape::cylinder::intersect(p, m_center, m_radius, m_half_height);
    }

    void transport(Particle<T>& p, RandomState& state) noexcept override
    {
        bool updateAtt = true;
        AttenuationValues<T> att;
        T attSumInv;

        bool cont = basicshape::cylinder::pointInside(p.pos, m_center, m_radius, m_half_height);
        while (cont) {
            if (updateAtt) {
                att = m_pmma.attenuationValues(p.energy);
                attSumInv = 1 / (att.sum() * m_pmma_density);
                updateAtt = false;
            }
            const auto stepLen = -std::log(state.randomUniform<T>()) * attSumInv; // cm
            const auto intLen = intersect(p).intersection; // this can not be nullopt
            const std::array<T, 2> tbox { T { 0 }, intLen };
            const auto intHoles = m_kdtree.intersect(p, tbox);

            if (intHoles.valid()) {
                if (intHoles.rayOriginIsInsideItem) {
                    const auto holeAtt = m_air.attenuationValues(p.energy);
                    const auto interactionProb = 1 - std::exp(-intHoles.intersection * holeAtt.sum() * m_air_density);
                    const auto intRes = interactions::interactForced(interactionProb, att, p, m_air, state);
                    m_dose[intHoles.item->index].scoreEnergy(intRes.energyImparted);
                    updateAtt = intRes.particleEnergyChanged;
                    // transport particle across hole (particle is most likely alive)
                    p.border_translate(intHoles.intersection);
                    cont = intRes.particleAlive && basicshape::cylinder::pointInside(p.pos, m_center, m_radius, m_half_height);
                } else {
                    if (stepLen < intHoles.intersection) {
                        // interaction happends
                        p.translate(stepLen);
                        const auto intRes = interactions::interact(att, p, m_pmma, state);
                        m_dose[0].scoreEnergy(intRes.energyImparted);
                        updateAtt = intRes.particleEnergyChanged;
                        cont = intRes.particleAlive;
                    } else {

                        // transport particle to hole
                        p.border_translate(intHoles.intersection);
                        // find distance of hole crossing
                        const auto dist = intHoles.item->intersect(p).intersection;
                        const auto holeAtt = m_air.attenuationValues(p.energy);
                        const auto interactionProb = 1 - std::exp(-dist * holeAtt.sum() * m_air_density);
                        const auto intRes = interactions::interactForced(interactionProb, att, p, m_air, state);
                        m_dose[intHoles.item->index].scoreEnergy(intRes.energyImparted);
                        updateAtt = intRes.particleEnergyChanged;
                        // transport particle across hole (particle is most likely alive)
                        p.border_translate(dist);
                        cont = intRes.particleAlive && basicshape::cylinder::pointInside(p.pos, m_center, m_radius, m_half_height);
                    }
                }
            } else {
                if (stepLen < intLen) {
                    // interaction happends
                    p.translate(stepLen);
                    const auto intRes = interactions::interact(att, p, m_pmma, state);
                    m_dose[0].scoreEnergy(intRes.energyImparted);
                    updateAtt = intRes.particleEnergyChanged;
                    cont = intRes.particleAlive;
                } else {
                    // transport to border
                    p.border_translate(intLen);
                    cont = false;
                }
            }
        }
    }

protected:
    struct CTDIAirHole {
        std::array<T, 3> pos;
        T radii;
        T half_height;
        std::uint8_t index;

        void translate(const std::array<T, 3>& d) noexcept
        {
            for (std::size_t i = 0; i < 3; ++i)
                pos[i] += d[i];
        }
        const std::array<T, 3>& center() const noexcept { return pos; }
        std::array<T, 6> AABB() const noexcept
        {
            std::array<T, 6> aabb {
                pos[0] - radii,
                pos[1] - radii,
                pos[2] - half_height,
                pos[0] + radii,
                pos[1] + radii,
                pos[2] + half_height,
            };
            return aabb;
        }
        auto intersect(const Particle<T>& p) const noexcept
        {
            const std::array<T, 6> aabb {
                pos[0] - radii,
                pos[1] - radii,
                pos[2] - half_height,
                pos[0] + radii,
                pos[1] + radii,
                pos[2] + half_height,
            };
            const auto aabbintersection = basicshape::AABB::intersect(p, aabb);
            if (aabbintersection.valid())
                return basicshape::cylinder::intersect(p, pos, radii, half_height);
            else
                return aabbintersection;
        }
    };

    static constexpr T holeEdgeDistance() noexcept
    {
        return T { 1 };
    }
    static constexpr T holeRadii() noexcept
    {
        return T { 0.5 };
    }

private:
    std::array<T, 3> m_center;
    T m_radius = 0;
    T m_half_height = 0;
    T m_pmma_density = 0;
    T m_air_density = 0;
    std::array<DoseScore<T>, 6> m_dose;
    StaticKDTree<3, T, CTDIAirHole> m_kdtree;
    Material2<T, NMaterialShells> m_pmma;
    Material2<T, NMaterialShells> m_air;
};

}
