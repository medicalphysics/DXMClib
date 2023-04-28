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

Copyright 2023 Erlend Andersen
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
class TG195World3Breast final : public WorldItemBase<T> {

public:
    TG195World3Breast()
        : WorldItemBase<T>()
        , m_skin_material(Material2<T, NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
        , m_tissue_material(Material2<T, NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        translateBox(m_dose_boxes[0], { T { 5 }, T { -5 }, T { 0 } });
        translateBox(m_dose_boxes[1], { T { 2 }, T { 0 }, T { 0 } });
        translateBox(m_dose_boxes[2], { T { 5 }, T { 0 }, T { 0 } });
        translateBox(m_dose_boxes[3], { T { 8 }, T { 0 }, T { 0 } });
        translateBox(m_dose_boxes[4], { T { 5 }, T { 5 }, T { 0 } });
        translateBox(m_dose_boxes[5], { T { 5 }, T { 0 }, T { -1.5f } });
        translateBox(m_dose_boxes[6], { T { 5 }, T { 0 }, T { 1.5f } });
    }

    void setTissueMaterial(const Material2<T, NMaterialShells>& material, T dens)
    {
        m_tissue_material = material;
        m_tissue_density = std::abs(dens);
    }

    void setSkinMaterial(const Material2<T, NMaterialShells>& material, T dens)
    {
        m_skin_material = material;
        m_skin_density = std::abs(dens);
    }

    void translate(const std::array<T, 3>& dist) override
    {
        m_center = vectormath::add(m_center, dist);
        for (auto& b : m_dose_boxes) {
            for (std::size_t i = 0; i < 3; ++i) {
                b.aabb[i] += dist[i];
                b.aabb[i + 3] += dist[i];
            }
        }
    }

    std::array<T, 3> center() const override
    {
        return m_center;
    }

    std::array<T, 6> AABB() const override
    {
        std::array<T, 6> aabb {
            m_center[0],
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
        for (auto& b : m_dose_boxes) {
            b.dose.clear();
        }
    }

    WorldIntersectionResult<T> intersect(const Particle<T>& p) const noexcept override
    {
        return intersectHalfCylindar(p, m_center, m_radius, m_halfHeight);
    }

    VisualizationIntersectionResult<T, WorldItemBase<T>> intersectVisualization(const Particle<T>& p) const noexcept override
    {
        const auto aabb = AABB();
        auto cyl = basicshape::cylinder::template intersectVisualization<T, WorldItemBase<T>>(p, m_center, m_radius, m_halfHeight);
        const auto box = basicshape::AABB::template intersectVisualization<T, WorldItemBase<T>>(p, aabb);
        if (cyl.valid() && box.valid()) {
            if (basicshape::AABB::pointInside(p.pos, aabb)) {
                cyl.intersection = std::min(cyl.intersection, box.intersection);
                cyl.rayOriginIsInsideItem = true;
            } else {
                if (cyl.intersection < box.intersection) {
                    cyl.normal = { -1, 0, 0 };
                    cyl.intersection = box.intersection;
                }
                cyl.rayOriginIsInsideItem = false;
            }
        } else {
            cyl.intersectionValid = false;
        }

        return cyl;
    }

    void transport(Particle<T>& p, RandomState& state) noexcept override
    {
        bool cont = basicshape::cylinder::pointInside(p.pos, m_center, m_radius, m_halfHeight) && basicshape::AABB::pointInside(p.pos, AABB());
        while (cont) {
            const auto intBreast = intersectHalfCylindar(p, m_center, m_radius, m_halfHeight);
            if (intBreast.valid()) {
                const auto intTissue = intersectHalfCylindar(p, m_center, m_radius - m_skin_thick, m_halfHeight - m_skin_thick);
                if (intTissue.valid()) {
                    if (intTissue.rayOriginIsInsideItem) {
                        // only inside tissue
                        const auto att = m_tissue_material.attenuationValues(p.energy);
                        const auto attSumInv = 1 / (att.sum() * m_tissue_density);
                        const auto stepLen = -std::log(state.randomUniform<T>()) * attSumInv;
                        if (stepLen < intTissue.intersection) {
                            p.translate(stepLen);
                            const auto intRes = interactions::template interact<T, NMaterialShells, LOWENERGYCORRECTION>(att, p, m_tissue_material, state);
                            cont = intRes.particleAlive;
                            scoreDose(p, intRes.energyImparted);
                        } else {
                            p.border_translate(intTissue.intersection);
                        }
                    } else {
                        // starts in skin and goes to tissue
                        const auto att = m_skin_material.attenuationValues(p.energy);
                        const auto attSumInv = 1 / (att.sum() * m_skin_density);
                        const auto stepLen = -std::log(state.randomUniform<T>()) * attSumInv;
                        if (stepLen < intTissue.intersection) {
                            p.translate(stepLen);
                            const auto intRes = interactions::template interact<T, NMaterialShells, LOWENERGYCORRECTION>(att, p, m_skin_material, state);
                            cont = intRes.particleAlive;
                            m_skin_dose.scoreEnergy(intRes.energyImparted);
                        } else {
                            p.border_translate(intTissue.intersection);
                        }
                    }
                } else {
                    // only intersects skin
                    const auto att = m_skin_material.attenuationValues(p.energy);
                    const auto attSumInv = 1 / (att.sum() * m_skin_density);
                    const auto stepLen = -std::log(state.randomUniform<T>()) * attSumInv;
                    if (stepLen < intBreast.intersection) {
                        p.translate(stepLen);
                        const auto intRes = interactions::template interact<T, NMaterialShells, LOWENERGYCORRECTION>(att, p, m_skin_material, state);
                        cont = intRes.particleAlive;
                        m_skin_dose.scoreEnergy(intRes.energyImparted);
                    } else {
                        p.border_translate(intBreast.intersection);
                        cont = false;
                    }
                }
            } else {
                cont = false;
            }
        }
    }

    const DoseScore<T>& dose(std::size_t index = 0) const override
    {
        if (index < m_dose_boxes.size())
            return m_dose_boxes[index].dose;
        else if (index == m_dose_boxes.size())
            return m_skin_dose;
        return m_dose;
    }

protected:
    struct DoseBoxChild {
        std::array<T, 6> aabb = { -1, -1, -.5f, 1, 1, .5f };
        DoseScore<T> dose;
    };
    static void translateBox(DoseBoxChild& box, const std::array<T, 3>& vec)
    {
        for (std::size_t i = 0; i < 3; ++i) {
            box.aabb[i] += vec[i];
            box.aabb[i + 3] += vec[i];
        }
    }

    void scoreDose(const Particle<T>& p, T energy)
    {
        bool hit = false;
        for (auto& b : m_dose_boxes) {
            if (!hit && basicshape::AABB::pointInside(p.pos, b.aabb)) {
                b.dose.scoreEnergy(energy);
                hit = true;
            }
        }
        m_dose.scoreEnergy(energy);
    }

    static WorldIntersectionResult<T> intersectHalfCylindar(const Particle<T>& p, const std::array<T, 3>& center, const T radius, const T halfHeight)
    {
        const std::array<T, 6> aabb = {
            center[0],
            center[1] - radius,
            center[2] - halfHeight,
            center[0] + radius,
            center[1] + radius,
            center[2] + halfHeight
        };
        auto cyl = basicshape::cylinder::template intersect(p, center, radius, halfHeight);
        const auto box = basicshape::AABB::template intersect(p, aabb);
        if (cyl.valid() && box.valid()) {
            if (basicshape::AABB::pointInside(p.pos, aabb)) {
                cyl.intersection = std::min(cyl.intersection, box.intersection);
            } else {
                cyl.intersection = std::max(cyl.intersection, box.intersection);
                cyl.rayOriginIsInsideItem = false;
            }
        } else {
            cyl.intersectionValid = false;
        }

        return cyl;
    }

private:
    T m_radius = 10;
    T m_skin_thick = 0.2f;
    T m_halfHeight = 2.5f;
    std::array<T, 3> m_center = { 0, 0, 0 };
    T m_skin_density = 1;
    T m_tissue_density = 1;
    Material2<T, NMaterialShells> m_skin_material;
    Material2<T, NMaterialShells> m_tissue_material;
    DoseScore<T> m_dose;
    DoseScore<T> m_skin_dose;
    std::array<DoseBoxChild, 7> m_dose_boxes;
};
}
