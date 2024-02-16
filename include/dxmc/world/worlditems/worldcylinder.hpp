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
#include "dxmc/world/basicshapes/cylinder.hpp"
#include "dxmc/world/dosescore.hpp"
#include "dxmc/world/worlditems/worlditembase.hpp"

#include <limits>
#include <optional>

namespace dxmc {

template <std::size_t NMaterialShells = 5, int LOWENERGYCORRECTION = 2>
class WorldCylinder final : public WorldItemBase {
public:
    WorldCylinder(double radius = 16, double height = 10, const std::array<double, 3>& center = { 0, 0, 0 }, const std::array<double, 3>& dir = { 0, 0, 1 })
        : WorldItemBase()
        , m_material(Material<double, NMaterialShells>::byNistName("Polymethyl Methacralate (Lucite, Perspex)").value())
    {
        m_cylinder.radius = std::abs(radius);
        m_cylinder.half_height = std::abs(height) / 2;
        m_cylinder.center = center;
        m_cylinder.direction = vectormath::normalized(dir);
        m_materialDensity = NISTMaterials<double>::density("Polymethyl Methacralate (Lucite, Perspex)");
    }

    void setMaterial(const Material<double, NMaterialShells>& material)
    {
        m_material = material;
    }

    void setMaterial(const Material<double, NMaterialShells>& material, double density)
    {
        m_material = material;
        setMaterialDensity(density);
    }

    void setMaterialDensity(double density)
    {
        m_materialDensity = std::abs(density);
    }

    bool setNistMaterial(const std::string& nist_name)
    {
        const auto mat = Material<double, NMaterialShells>::byNistName(nist_name);
        if (mat) {
            m_material = mat.value();
            m_materialDensity = NISTMaterials<double>::density(nist_name);
            return true;
        }
        return false;
    }

    void translate(const std::array<double, 3>& dist) override
    {
        m_cylinder.center = vectormath::add(m_cylinder.center, dist);
    }

    std::array<double, 3> center() const override
    {
        return m_cylinder.center;
    }

    std::array<double, 6> AABB() const override
    {
        return basicshape::cylinder::cylinderAABB(m_cylinder);
    }

    void setRadius(double r)
    {
        m_cylinder.radius = std::abs(r);
    }

    void setHeight(double h)
    {
        m_cylinder.half_height = std::abs(h / 2);
    }

    WorldIntersectionResult intersect(const Particle& p) const noexcept override
    {
        return basicshape::cylinder::intersect(p, m_cylinder);
    }

    VisualizationIntersectionResult<WorldItemBase> intersectVisualization(const Particle& p) const noexcept override
    {
        auto inter = basicshape::cylinder::template intersectVisualization<WorldItemBase>(p, m_cylinder);
        if (inter.valid())
            inter.value = m_dose.dose();
        return inter;
    }

    void transport(Particle& p, RandomState& state) noexcept override
    {
        bool cont = basicshape::cylinder::pointInside(p.pos, m_cylinder);
        bool updateAtt = true;
        AttenuationValues<double> att;
        double attSumInv;
        while (cont) {
            if (updateAtt) {
                att = m_material.attenuationValues(p.energy);
                attSumInv = 1 / (att.sum() * m_materialDensity);
                updateAtt = false;
            }
            const auto stepLen = -std::log(state.randomUniform()) * attSumInv; // cm
            const auto intLen = intersect(p);

            if (stepLen < intLen.intersection) {
                // interaction happends
                p.translate(stepLen);
                const auto intRes = interactions::template interact<NMaterialShells, LOWENERGYCORRECTION>(att, p, m_material, state);
                m_energyScored.scoreEnergy(intRes.energyImparted);
                updateAtt = intRes.particleEnergyChanged;
                cont = intRes.particleAlive;
            } else {
                // transport to border
                p.border_translate(intLen.intersection);
                cont = false;
            }
        }
    }

    const EnergyScore& energyScored(std::size_t index = 0) const override
    {
        return m_energyScored;
    }

    void clearEnergyScored() override
    {
        m_energyScored.clear();
    }

    void addEnergyScoredToDoseScore(double calibration_factor = 1) override
    {
        m_dose.addScoredEnergy(m_energyScored, m_cylinder.volume(), m_materialDensity, calibration_factor);
    }

    const DoseScore& doseScored(std::size_t index = 0) const override
    {
        return m_dose;
    }

    void clearDoseScored() override
    {
        m_dose.clear();
    }

protected:
private:
    basicshape::cylinder::Cylinder m_cylinder;
    double m_materialDensity = 1;
    Material<double, NMaterialShells> m_material;
    EnergyScore m_energyScored;
    DoseScore m_dose;
};

}
