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
public:
    TG195World42(T radius = T { 16 }, T height = T { 10 }, const std::array<T, 3>& pos = { 0, 0, 0 })
        : WorldItemBase<T>()
        , m_radius(radius)
        , m_halfHeight(height * T { 0.5 })
        , m_center(pos)
        , m_material(Material2<T, NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        m_materialDensity = NISTMaterials<T>::density("Air, Dry (near sea level)");

        m_centerChild.center = m_center;
        m_centerChild.radii = T { 0.5 };
        m_centerChild.halfHeight = 5;
        m_periferyChild.center = { m_center[0] - (m_radius - 1), m_center[1], m_center[2] };
        m_periferyChild.radii = T { 0.5 };
        m_periferyChild.halfHeight = 5;
    }

    void setMaterial(const Material2<T, NMaterialShells>& material)
    {
        m_material = material;
    }

    void setMaterialDensity(T density) { m_materialDensity = density; }

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
            m_center[2] - m_halfHeight,
            m_center[0] + m_radius,
            m_center[1] + m_radius,
            m_center[2] + m_halfHeight
        };
        return aabb;
    }

    void clearDose() override
    {
        m_dose.clear();
        m_centerChild.dose.clear();
        m_periferyChild.dose.clear();
    }

    WorldIntersectionResult<T> intersect(const Particle<T>& p) const noexcept override
    {
        return basicshape::cylinder::intersect(p, m_center, m_radius, m_halfHeight);
    }

    void transport(Particle<T>& p, RandomState& state) noexcept override
    {
        bool cont = basicshape::cylinder::pointInside(p.pos, m_center, m_radius, m_halfHeight);
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

                if (basicshape::cylinder::pointInside(p.pos, m_centerChild.center, m_centerChild.radii, m_centerChild.halfHeight)) {
                    m_centerChild.dose.scoreEnergy(intRes.energyImparted);
                } else if (basicshape::cylinder::pointInside(p.pos, m_periferyChild.center, m_periferyChild.radii, m_periferyChild.halfHeight)) {
                    m_periferyChild.dose.scoreEnergy(intRes.energyImparted);
                } else {
                    m_dose.scoreEnergy(intRes.energyImparted);
                }

                updateAtt = intRes.particleEnergyChanged;
                cont = intRes.particleAlive;
            } else {
                // transport to border
                p.border_translate(intLen.intersection);
                cont = false;
            }
        }
    }

    const DoseScore<T>& doseCenterCylinder() const
    {
        return m_centerChild.dose;
    }
    const DoseScore<T>& dosePeriferyCylinder() const
    {
        return m_periferyChild.dose;
    }
    const DoseScore<T>& dose(std::size_t index = 0) const override
    {
        return m_dose;
    }

protected:
    struct CylinderChild {
        std::array<T, 3> center;
        T radii;
        T halfHeight;
        DoseScore<T> dose;
    };

private:
    T m_radius = 0;
    T m_halfHeight = 0;
    std::array<T, 3> m_center;
    T m_materialDensity = 1;
    Material2<T, NMaterialShells> m_material;
    DoseScore<T> m_dose;
    CylinderChild m_centerChild;
    CylinderChild m_periferyChild;
};

}
