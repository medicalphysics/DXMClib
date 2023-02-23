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

template <Floating T, int NMaterialShells = 5>
class DepthDose final : public WorldItemBase<T> {
public:
    DepthDose(T radius = T { 16 }, T height = T { 10 }, std::size_t resolution = 100, const std::array<T, 3>& pos = { 0, 0, 0 })
        : WorldItemBase<T>()
        , m_radius(radius)
        , m_half_height(height * T { 0.5 })
        , m_center(pos)
        , m_material(Material2<T, NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        m_materialDensity = NISTMaterials<T>::density("Air, Dry (near sea level)");
        m_dose.resize(resolution);
    }

    void setMaterial(const Material2<T, NMaterialShells>& material)
    {
        m_material = material;
    }

    void setMaterialDensity(T density) { m_materialDensity = density; }

    bool setNistMaterial(const std::string& nist_name)
    {
        const auto mat = Material2<T, NMaterialShells>::byNistName(nist_name);
        if (mat) {
            m_material = mat.value();
            m_materialDensity = NISTMaterials<T>::density(nist_name);
            return true;
        }
        return false;
    }
    const Material2<T, NMaterialShells>& material() const
    {
        return m_material;
    }
    T density() const { return m_materialDensity; }

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
            m_center[2] - m_half_height,
            m_center[0] + m_radius,
            m_center[1] + m_radius,
            m_center[2] + m_half_height
        };
        return aabb;
    }

    std::optional<T> intersectForward(const Particle<T>& p) const noexcept override
    {
        return basicshape::cylinder::intersectForward(p, m_center, m_radius, m_half_height);
    }

    void transport(Particle<T>& p, RandomState& state) noexcept override
    {
        bool cont = basicshape::cylinder::pointInside(p.pos, m_center, m_radius, m_half_height);
        bool updateAtt = false;
        auto att = m_material.attenuationValues(p.energy);
        auto attSumInv = 1 / (att.sum() * m_materialDensity);
        while (cont) {
            if (updateAtt) {
                att = m_material.attenuationValues(p.energy);
                attSumInv = 1 / (att.sum() * m_materialDensity);
                updateAtt = false;
            }
            const auto stepLen = -std::log(state.randomUniform<T>()) * attSumInv; // cm
            const auto intLen = intersectForward(p).value(); // this can not be nullopt

            if (stepLen < intLen) {
                // interaction happends
                p.translate(stepLen);
                const auto intRes = interactions::interact(att, p, m_material, state);
                updateAtt = intRes.particleEnergyChanged;
                cont = intRes.particleAlive;

                const auto dose_ind_f = (p.pos[2] - (m_center[2] - m_half_height)) * m_dose.size() / (m_half_height * 2);
                const auto ind = std::clamp(static_cast<std::size_t>(dose_ind_f), std::size_t { 0 }, m_dose.size() - 1);
                m_dose[ind].scoreEnergy(intRes.energyImparted);
            } else {
                // transport to border
                p.border_translate(intLen);
                cont = false;
            }
        }
    }

    const std::vector<DoseScore<T>>& doseVector() const
    {
        return m_dose;
    }
    const std::vector<T> depthVector() const
    {
        std::vector<T> depth(m_dose.size());
        const auto step = (2 * m_half_height) / depth.size();
        const auto start = m_center[2] - m_half_height + step / 2;
        for (std::size_t i = 0; i < depth.size(); ++i) {
            depth[i] = start + step * i;
        }
        return depth;
    }

    const DoseScore<T>& dose(std::size_t index = 0) const override
    {
        return m_dose[index];
    }

protected:
private:
    T m_radius = 0;
    T m_half_height = 0;
    std::array<T, 3> m_center;
    T m_materialDensity = 1;
    Material2<T, NMaterialShells> m_material;
    std::vector<DoseScore<T>> m_dose;
};

}
