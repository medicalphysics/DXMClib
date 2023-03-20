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
#include "dxmc/material/massenergytransfer.hpp"
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
    struct CylinderChild {
        std::array<T, 3> center = { 0, 0, 0 };
        T radii = 0;
        T halfHeight = 0;
        DoseScore<T> dose;
        std::array<T, 6> aabb = { 0, 0, 0, 0, 0, 0 };
    };

public:
    TG195World42(T radius = T { 16 }, T height = T { 10 }, const std::array<T, 3>& pos = { 0, 0, 0 })
        : WorldItemBase<T>()
        , m_radius(radius)
        , m_halfHeight(height * T { 0.5 })
        , m_center(pos)
        , m_material(Material2<T, NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
        , m_massEnergyTransfer(MassEnergyTransfer<T>("Air, Dry(near sea level)"))
    {
        m_materialDensity = NISTMaterials<T>::density("Air, Dry (near sea level)");

        m_centerChild.center = m_center;
        m_centerChild.radii = T { 0.5 };
        m_centerChild.halfHeight = 5;
        m_periferyChild.center = { m_center[0] - (m_radius - 1), m_center[1], m_center[2] };
        m_periferyChild.radii = T { 0.5 };
        m_periferyChild.halfHeight = 5;
        for (std::size_t i = 0; i < 3; ++i) {
            m_centerChild.aabb[i] = m_centerChild.center[i] - m_centerChild.radii;
            m_centerChild.aabb[i + 3] = m_centerChild.center[i] + m_centerChild.radii;
            m_periferyChild.aabb[i] = m_periferyChild.center[i] - m_periferyChild.radii;
            m_periferyChild.aabb[i + 3] = m_periferyChild.center[i] + m_periferyChild.radii;
            if (i == 2) {
                m_centerChild.aabb[i] = m_centerChild.center[i] - m_centerChild.halfHeight;
                m_centerChild.aabb[i + 3] = m_centerChild.center[i] + m_centerChild.halfHeight;
                m_periferyChild.aabb[i] = m_periferyChild.center[i] - m_periferyChild.halfHeight;
                m_periferyChild.aabb[i + 3] = m_periferyChild.center[i] + m_periferyChild.halfHeight;
            }
        }
    }

    void setMaterial(const Material2<T, NMaterialShells>& material, const MassEnergyTransfer<T>& massEn)
    {
        m_material = material;
        m_massEnergyTransfer = massEn;
    }

    void setMaterialDensity(T density) { m_materialDensity = density; }

    void translate(const std::array<T, 3>& dist) override
    {
        m_center = vectormath::add(m_center, dist);
        m_centerChild.center = vectormath::add(m_centerChild.center, dist);
        m_periferyChild.center = vectormath::add(m_periferyChild.center, dist);
        for (std::size_t i = 0; i < 3; ++i) {
            m_periferyChild.aabb[i] += dist[i];
            m_periferyChild.aabb[i + 3] += dist[i];
            m_centerChild.aabb[i] += dist[i];
            m_centerChild.aabb[i + 3] += dist[i];
        }
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

    std::pair<CylinderChild*, std::array<T, 2>> intersectChild(const Particle<T>& p, RandomState& state)
    {
        CylinderChild* child = nullptr;
        auto intersectChildC = basicshape::AABB::intersectForwardInterval(p, m_centerChild.aabb);
        if (intersectChildC) {
            intersectChildC = basicshape::cylinder::intersectForwardInterval(
                p, m_centerChild.center, m_centerChild.radii, m_centerChild.halfHeight);
        }
        auto intersectChildP = basicshape::AABB::intersectForwardInterval(p, m_periferyChild.aabb);
        if (intersectChildP) {
            intersectChildP = basicshape::cylinder::intersectForwardInterval(
                p, m_periferyChild.center, m_periferyChild.radii, m_periferyChild.halfHeight);
        }
        std::array<T, 2> t = { 0, 0 };
        if (intersectChildC && intersectChildP) {
            // select a random child
            if (state.randomUniform<T>() < T { 0.5 }) {
                child = &m_centerChild;
                t = intersectChildC.value();
            } else {
                child = &m_periferyChild;
                t = intersectChildP.value();
            }
        } else if (intersectChildC) {
            child = &m_centerChild;
            t = intersectChildC.value();
        } else if (intersectChildC) {
            child = &m_periferyChild;
            t = intersectChildP.value();
        }
        return std::make_pair(child, t);
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

            T stepLen;
            // auto stepLen = -std::log(state.randomUniform<T>()) * attSumInv; // cm
            constexpr int VARRED = 1;
            // forced interactions
            if constexpr (VARRED == 1) {
                // Forced interactions by modifying path lenght sampling to intersection interval of ROI
                auto [child_intersected, t_child] = intersectChild(p, state);
                if (child_intersected && state.randomUniform<T>() > T { 0.5 }) { // we pass an object of interest
                    const auto p0 = T { 1 };
                    // const auto p0 = std::exp(-t_child[0] * att.sum() * m_materialDensity);
                    const auto p1 = std::exp(-t_child[1] * att.sum() * m_materialDensity);
                    const auto pp = p0 - p1;

                    stepLen = -attSumInv * std::log(1 - (1 - p1) * state.randomUniform<T>());
                    p.weight *= pp;
                } else {
                    stepLen = -std::log(state.randomUniform<T>()) * attSumInv;
                }
            } else if constexpr (VARRED == 2) {
                stepLen = -std::log(state.randomUniform<T>()) * attSumInv;

                if (basicshape::AABB::intersectForwardInterval(p, m_centerChild.aabb)) {
                    auto intersectChild = basicshape::cylinder::intersectForwardInterval(
                        p, m_centerChild.center, m_centerChild.radii, m_centerChild.halfHeight);
                    if (intersectChild) {
                        auto& t = intersectChild.value();
                        t[1] = std::min(stepLen, t[1]);
                        if (t[0] < t[1]) {
                            const auto uen = m_massEnergyTransfer(p.energy);
                            const auto EI = p.energy * m_materialDensity * uen * (t[1] - t[0]) * p.weight;
                            m_centerChild.dose.scoreEnergy(EI);
                        }
                    }
                }
                if (basicshape::AABB::intersectForwardInterval(p, m_periferyChild.aabb)) {
                    auto intersectChild = basicshape::cylinder::intersectForwardInterval(
                        p, m_periferyChild.center, m_periferyChild.radii, m_periferyChild.halfHeight);
                    if (intersectChild) {
                        auto& t = intersectChild.value();
                        t[1] = std::min(stepLen, t[1]);
                        if (t[0] < t[1]) {
                            const auto uen = m_massEnergyTransfer(p.energy);
                            const auto EI = p.energy * m_materialDensity * uen * (t[1] - t[0]) * p.weight;
                            m_periferyChild.dose.scoreEnergy(EI);
                        }
                    }
                }

            } else {
                stepLen = -std::log(state.randomUniform<T>()) * attSumInv;
            }

            const auto intLen = intersect(p);
            if (stepLen < intLen.intersection) {
                // interaction happends
                p.translate(stepLen);
                const auto intRes = interactions::template interact<T, NMaterialShells, LOWENERGYCORRECTION>(att, p, m_material, state);
                if constexpr (VARRED != 2) {
                    if (basicshape::cylinder::pointInside(p.pos, m_centerChild.center, m_centerChild.radii, m_centerChild.halfHeight)) {
                        m_centerChild.dose.scoreEnergy(intRes.energyImparted);
                    } else if (basicshape::cylinder::pointInside(p.pos, m_periferyChild.center, m_periferyChild.radii, m_periferyChild.halfHeight)) {
                        m_periferyChild.dose.scoreEnergy(intRes.energyImparted);
                    } else {
                        // m_dose.scoreEnergy(intRes.energyImparted);
                    }
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

private:
    T m_radius = 0;
    T m_halfHeight = 0;
    std::array<T, 3> m_center;
    T m_materialDensity = 1;
    Material2<T, NMaterialShells> m_material;
    MassEnergyTransfer<T> m_massEnergyTransfer;
    DoseScore<T> m_dose;
    CylinderChild m_centerChild;
    CylinderChild m_periferyChild;
};
}
