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
#include "dxmc/world/basicshapes/cylinder.hpp"
#include "dxmc/world/statickdtree.hpp"
#include "dxmc/world/worlditems/worlditembase.hpp"

#include <array>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

namespace dxmc {

template <int NMaterialShells = 5, int LOWENERGYCORRECTION = 2>
class CTDIPhantom final : public WorldItemBase {
public:
    CTDIPhantom(double radius = 16, double height = 15, const std::array<double, 3>& pos = { 0, 0, 0 }, const std::array<double, 3>& direction = { 0, 0, 1 })
        : WorldItemBase()
        , m_pmma(Material<NMaterialShells>::byNistName("Polymethyl Methacralate (Lucite, Perspex)").value())
        , m_air(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        radius = std::abs(radius);
        height = std::max(std::abs(height), holeHeight());
        m_cylinder.center = pos;
        m_cylinder.direction = direction;
        m_cylinder.radius = radius;
        m_cylinder.half_height = height / 2;

        m_pmma_density = NISTMaterials::density("Polymethyl Methacralate (Lucite, Perspex)");
        m_air_density = NISTMaterials::density("Air, Dry (near sea level)");

        std::vector<CTDIAirHole> holes;
        holes.reserve(5);
        std::array<double, 10> positions = { 0, 0, 0, 1, 1, 0, 0, -1, -1, 0 };
        std::uint8_t index = 0;
        for (std::size_t i = 0; i < 10; i = i + 2) {
            const auto x = positions[i];
            const auto y = positions[i + 1];
            CTDIAirHole hole;
            hole.cylinder.center = m_cylinder.center;
            hole.cylinder.half_height = holeHeight() / 2;
            hole.cylinder.direction = m_cylinder.direction;
            hole.cylinder.radius = holeRadii();
            hole.index = ++index;
            hole.cylinder.center[0] += x * (radius - holeEdgeDistance());
            hole.cylinder.center[1] += y * (radius - holeEdgeDistance());
            holes.push_back(hole);
        }
        m_kdtree = StaticKDTree<3, CTDIAirHole>(holes);
    }

    void translate(const std::array<double, 3>& dist) noexcept override
    {
        m_cylinder.center = vectormath::add(m_cylinder.center, dist);
        m_kdtree.translate(dist);
    }

    std::array<double, 3> center() const noexcept override
    {
        return m_cylinder.center;
    }

    void setHoleMaterial(const std::string& nistName, double density)
    {
        auto m = Material<double, NMaterialShells>::byNistName(nistName);
        if (m) {
            m_air = m.value();
            m_air_density = density;
        }
    }

    const EnergyScore& energyScored(std::size_t index = 0) const override
    {
        return m_energyScore[index];
    }

    void clearEnergyScored() override
    {
        for (auto& d : m_energyScore) {
            d.clear();
        }
    }

    void addEnergyScoredToDoseScore(double calibration_factor = 1) override
    {
        const auto holeVolume = holeHeight() * std::numbers::pi_v<double> * holeRadii() * holeRadii();

        for (std::size_t i = 1; i < m_energyScore.size(); ++i) {
            m_dose[i].addScoredEnergy(m_energyScore[i], holeVolume, m_air_density, calibration_factor);
        }

        const auto totalVolume = m_cylinder.radius * m_cylinder.radius * 2 * m_cylinder.half_height * std::numbers::pi_v<double>;
        const auto pmmaVolume = totalVolume - 5 * holeVolume;
        m_dose[0].addScoredEnergy(m_energyScore[0], pmmaVolume, m_pmma_density, calibration_factor);
    }

    const DoseScore& doseScored(std::size_t index = 0) const override
    {
        return m_dose[index];
    }

    const double centerDoseScored() const
    {
        return m_dose[1].dose();
    }

    const double pheriferyDoseScored() const
    {
        return (m_dose[2].dose() + m_dose[3].dose() + m_dose[4].dose() + m_dose[5].dose()) / 4;
    }

    void clearDoseScored() override
    {
        for (auto& d : m_dose) {
            d.clear();
        }
    }

    std::array<double, 6> AABB() const noexcept override
    {
        return basicshape::cylinder::cylinderAABB(m_cylinder);
    }

    WorldIntersectionResult intersect(const Particle& p) const noexcept override
    {
        return basicshape::cylinder::intersect(p, m_cylinder);
    }

    VisualizationIntersectionResult<WorldItemBase> intersectVisualization(const Particle& p) const noexcept override
    {
        auto res = basicshape::cylinder::intersectVisualization<WorldItemBase>(p, m_cylinder);
        if (res.valid()) {
            const std::array tbox = { res.intersection, res.intersection + m_cylinder.radius };
            auto holes = m_kdtree.intersect(p, tbox);
            if (holes.valid()) {
                res.value = m_dose[holes.item->index].dose();
            } else {
                res.value = m_dose[0].dose();
            }
        }
        return res;
    }

    void transport(Particle& p, RandomState& state) noexcept override
    {
        bool updateAtt = true;
        AttenuationValues att;
        double attSumInv;
        bool cont = basicshape::cylinder::pointInside(p.pos, m_cylinder);
        while (cont) {
            if (updateAtt) {
                att = m_pmma.attenuationValues(p.energy);
                attSumInv = 1 / (att.sum() * m_pmma_density);
                updateAtt = false;
            }
            const auto stepLen = -std::log(state.randomUniform()) * attSumInv; // cm
            const auto intCTDI = intersect(p); // this can not be nullopt
            const auto intLen = intCTDI.intersection;
            const std::array<double, 2> tbox { 0.0, intLen };
            const auto intHoles = m_kdtree.intersect(p, tbox);

            if (intHoles.valid()) {
                if (intHoles.rayOriginIsInsideItem) {
                    const auto intRes = interactions::template interactForced<NMaterialShells, LOWENERGYCORRECTION>(intHoles.intersection, m_air_density, p, m_air, state);
                    m_energyScore[intHoles.item->index].scoreEnergy(intRes.energyImparted);
                    updateAtt = intRes.particleEnergyChanged;
                    cont = intRes.particleAlive && basicshape::cylinder::pointInside(p.pos, m_cylinder);
                } else {
                    if (stepLen < intHoles.intersection) {
                        // interaction happends before particle hit hole
                        p.translate(stepLen);
                        const auto intRes = interactions::template interact<NMaterialShells, LOWENERGYCORRECTION>(att, p, m_pmma, state);
                        m_energyScore[0].scoreEnergy(intRes.energyImparted);
                        updateAtt = intRes.particleEnergyChanged;
                        cont = intRes.particleAlive;
                    } else {
                        // transport particle to hole
                        p.border_translate(intHoles.intersection);
                        // find distance of hole crossing
                        const auto dist = intHoles.item->intersect(p).intersection;
                        const auto intRes = interactions::template interactForced<NMaterialShells, LOWENERGYCORRECTION>(dist, m_air_density, p, m_air, state);
                        m_energyScore[intHoles.item->index].scoreEnergy(intRes.energyImparted);
                        updateAtt = intRes.particleEnergyChanged;
                        cont = intRes.particleAlive && basicshape::cylinder::pointInside(p.pos, m_cylinder);
                    }
                }
            } else {
                if (stepLen < intLen) {
                    // interaction happends
                    p.translate(stepLen);
                    const auto intRes = interactions::interact(att, p, m_pmma, state);
                    m_energyScore[0].scoreEnergy(intRes.energyImparted);
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

    static constexpr double holeRadii() noexcept
    {
        return 0.5;
    }
    static constexpr double holeHeight() noexcept
    {
        return 10.0;
    }

protected:
    struct CTDIAirHole {
        basicshape::cylinder::Cylinder cylinder;
        std::uint8_t index = 0;

        void translate(const std::array<double, 3>& d) noexcept
        {
            for (std::size_t i = 0; i < 3; ++i)
                cylinder.center[i] += d[i];
        }
        const std::array<double, 3>& center() const noexcept { return cylinder.center; }
        std::array<double, 6> AABB() const noexcept
        {
            return basicshape::cylinder::cylinderAABB(cylinder);
        }
        auto intersect(const Particle& p) const noexcept
        {
            return basicshape::cylinder::intersect(p, cylinder);
        }
    };

    static constexpr double holeEdgeDistance() noexcept
    {
        return 1.0;
    }

private:
    basicshape::cylinder::Cylinder m_cylinder;
    double m_pmma_density = 0;
    double m_air_density = 0;
    std::array<EnergyScore, 6> m_energyScore;
    std::array<DoseScore, 6> m_dose;
    StaticKDTree<3, CTDIAirHole> m_kdtree;
    Material<NMaterialShells> m_pmma;
    Material<NMaterialShells> m_air;
};
}
