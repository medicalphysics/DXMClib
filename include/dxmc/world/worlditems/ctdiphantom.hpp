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

template <Floating T, int NMaterialShells = 5, int LOWENERGYCORRECTION = 2>
class CTDIPhantom final : public WorldItemBase<T> {
public:
    CTDIPhantom(T radius = T { 16 }, T height = T { 15 }, const std::array<T, 3>& pos = { 0, 0, 0 }, const std::array<T, 3>& direction = { 0, 0, 1 }, T holeHeight = T { 10 })
        : WorldItemBase<T>()
        , m_pmma(Material<T, NMaterialShells>::byNistName("Polymethyl Methacralate (Lucite, Perspex)").value())
        , m_air(Material<T, NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        radius = std::abs(radius);
        height = std::abs(height);
        m_cylinder.center = pos;
        m_cylinder.direction = direction;
        m_cylinder.radius = radius;
        m_cylinder.half_height = height / 2;

        m_pmma_density = NISTMaterials<T>::density("Polymethyl Methacralate (Lucite, Perspex)");
        m_air_density = NISTMaterials<T>::density("Air, Dry (near sea level)");

        holeHeight = std::min(height, std::abs(holeHeight));

        std::vector<CTDIAirHole> holes;
        holes.reserve(5);
        std::array<T, 10> positions = { 0, 0, 0, 1, 1, 0, 0, -1, -1, 0 };
        std::uint8_t index = 0;
        for (std::size_t i = 0; i < 10; i = i + 2) {
            const auto x = positions[i];
            const auto y = positions[i + 1];
            CTDIAirHole hole;
            hole.cylinder.center = m_cylinder.center;
            hole.cylinder.half_height = holeHeight / 2;
            hole.cylinder.direction = m_cylinder.direction;
            hole.cylinder.radius = holeRadii();
            hole.index = ++index;
            hole.cylinder.center[0] += x * (radius - holeEdgeDistance());
            hole.cylinder.center[1] += y * (radius - holeEdgeDistance());
            holes.push_back(hole);
        }
        m_kdtree = StaticKDTree<3, T, CTDIAirHole>(holes);
    }

    void translate(const std::array<T, 3>& dist) noexcept override
    {
        m_cylinder.center = vectormath::add(m_cylinder.center, dist);
        m_kdtree.translate(dist);
    }
    std::array<T, 3> center() const noexcept override
    {
        return m_cylinder.center;
    }
    void clearDose() override
    {
        for (auto& d : m_dose) {
            d.clear();
        }
    }

    void setHoleMaterial(const std::string& nistName, T density)
    {
        auto m = Material<T, NMaterialShells>::byNistName(nistName);
        if (m) {
            m_air = m.value();
            m_air_density = density;
        }
    }

    const EnergyScore<T>& dose(std::size_t index = 0) const override
    {
        return m_dose[index];
    }

    std::array<T, 6> AABB() const noexcept override
    {
        return basicshape::cylinder::cylinderAABB(m_cylinder);
    }

    WorldIntersectionResult<T> intersect(const Particle<T>& p) const noexcept override
    {
        return basicshape::cylinder::intersect(p, m_cylinder);
    }

    VisualizationIntersectionResult<T, WorldItemBase<T>> intersectVisualization(const Particle<T>& p) const noexcept override
    {
        return basicshape::cylinder::intersectVisualization<T, WorldItemBase<T>>(p, m_cylinder);
    }

    void transport(Particle<T>& p, RandomState& state) noexcept override
    {
        bool updateAtt = true;
        AttenuationValues<T> att;
        T attSumInv;

        bool cont = basicshape::cylinder::pointInside(p.pos, m_cylinder);
        while (cont) {
            if (updateAtt) {
                att = m_pmma.attenuationValues(p.energy);
                attSumInv = 1 / (att.sum() * m_pmma_density);
                updateAtt = false;
            }
            const auto stepLen = -std::log(state.randomUniform<T>()) * attSumInv; // cm
            const auto intCTDI = intersect(p); // this can not be nullopt
            const auto intLen = intCTDI.intersection;
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
                    cont = intRes.particleAlive && basicshape::cylinder::pointInside(p.pos, m_cylinder);
                } else {
                    if (stepLen < intHoles.intersection) {
                        // interaction happends
                        p.translate(stepLen);
                        const auto intRes = interactions::template interact<T, NMaterialShells, LOWENERGYCORRECTION>(att, p, m_pmma, state);
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
                        const auto intRes = interactions::template interactForced<T, NMaterialShells, LOWENERGYCORRECTION>(interactionProb, att, p, m_air, state);
                        m_dose[intHoles.item->index].scoreEnergy(intRes.energyImparted);
                        updateAtt = intRes.particleEnergyChanged;
                        // transport particle across hole (particle is most likely alive)
                        p.border_translate(dist);
                        cont = intRes.particleAlive && basicshape::cylinder::pointInside(p.pos, m_cylinder);
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
        basicshape::cylinder::Cylinder<T> cylinder;
        std::uint8_t index = 0;

        void translate(const std::array<T, 3>& d) noexcept
        {
            for (std::size_t i = 0; i < 3; ++i)
                cylinder.center[i] += d[i];
        }
        const std::array<T, 3>& center() const noexcept { return cylinder.center; }
        std::array<T, 6> AABB() const noexcept
        {
            return basicshape::cylinder::cylinderAABB(cylinder);
        }
        auto intersect(const Particle<T>& p) const noexcept
        {
            return basicshape::cylinder::intersect(p, cylinder);
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
    basicshape::cylinder::Cylinder<T> m_cylinder;
    T m_pmma_density = 0;
    T m_air_density = 0;
    std::array<EnergyScore<T>, 6> m_dose;
    StaticKDTree<3, T, CTDIAirHole> m_kdtree;
    Material<T, NMaterialShells> m_pmma;
    Material<T, NMaterialShells> m_air;
};

}
