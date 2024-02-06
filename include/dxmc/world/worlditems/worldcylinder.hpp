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
#include "dxmc/world/dosescore.hpp"
#include "dxmc/world/worlditems/worlditembase.hpp"

#include <limits>
#include <optional>

namespace dxmc {

template <Floating T, std::size_t NMaterialShells = 5, int LOWENERGYCORRECTION = 2>
class WorldCylinder final : public WorldItemBase<T> {
public:
    WorldCylinder(T radius = T { 16 }, T height = T { 10 }, const std::array<T, 3>& center = { 0, 0, 0 }, const std::array<T, 3>& dir = { 0, 0, 1 })
        : WorldItemBase<T>()
        , m_material(Material<T, NMaterialShells>::byNistName("Polymethyl Methacralate (Lucite, Perspex)").value())
    {
        m_cylinder.radius = std::abs(radius);
        m_cylinder.half_height = std::abs(height) / 2;
        m_cylinder.center = center;
        m_cylinder.direction = vectormath::normalized(dir);
        m_materialDensity = NISTMaterials<T>::density("Polymethyl Methacralate (Lucite, Perspex)");
    }

    void setMaterial(const Material<T, NMaterialShells>& material)
    {
        m_material = material;
    }

    void setMaterial(const Material<T, NMaterialShells>& material, T density)
    {
        m_material = material;
        setMaterialDensity(density);
    }

    void setMaterialDensity(T density)
    {
        m_materialDensity = std::abs(density);
    }

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
        m_cylinder.center = vectormath::add(m_cylinder.center, dist);
    }

    std::array<T, 3> center() const override
    {
        return m_cylinder.center;
    }

    std::array<T, 6> AABB() const override
    {
        return basicshape::cylinder::cylinderAABB(m_cylinder);
    }

    void setRadius(const T r)
    {
        m_cylinder.radius = std::abs(r);
    }

    void setHeight(T h)
    {
        m_cylinder.half_height = std::abs(h / 2);
    }

    WorldIntersectionResult<T> intersect(const Particle<T>& p) const noexcept override
    {
        return basicshape::cylinder::intersect(p, m_cylinder);
    }

    VisualizationIntersectionResult<T, WorldItemBase<T>> intersectVisualization(const Particle<T>& p) const noexcept override
    {
        auto inter = basicshape::cylinder::template intersectVisualization<T, WorldItemBase<T>>(p, m_cylinder);
        if (inter.valid())
            inter.value = m_dose.dose();
        return inter;
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
            const auto stepLen = -std::log(state.randomUniform<T>()) * attSumInv; // cm
            const auto intLen = intersect(p);

            if (stepLen < intLen.intersection) {
                // interaction happends
                p.translate(stepLen);
                const auto intRes = interactions::template interact<T, NMaterialShells, LOWENERGYCORRECTION>(att, p, m_material, state);
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
        m_dose.addScoredEnergy(m_energyScored, m_cylinder.volume(), m_materialDensity, calibration_factor);
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
private:
    basicshape::cylinder::Cylinder<T> m_cylinder;
    T m_materialDensity = 1;
    Material<T, NMaterialShells> m_material;
    EnergyScore<T> m_energyScored;
    DoseScore<T> m_dose;
};

}
