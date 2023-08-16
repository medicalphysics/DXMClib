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

template <Floating T, std::size_t NMaterialShells = 5, int Lowenergycorrection = 2>
class DepthDose final : public WorldItemBase<T> {
public:
    DepthDose(T radius = T { 16 }, T height = T { 10 }, std::size_t resolution = 100, const std::array<T, 3>& pos = { 0, 0, 0 })
        : WorldItemBase<T>()
        , m_material(Material<T, NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        m_cylinder.center = pos;
        m_cylinder.direction = { 0, 0, 1 };
        m_cylinder.radius = radius;
        m_cylinder.half_height = height / 2;

        m_materialDensity = NISTMaterials<T>::density("Air, Dry (near sea level)");
        m_energyScored.resize(resolution);
    }

    void setMaterial(const Material<T, NMaterialShells>& material)
    {
        m_material = material;
    }

    void setMaterialDensity(T density) { m_materialDensity = density; }

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
    const Material<T, NMaterialShells>& material() const
    {
        return m_material;
    }
    T density() const { return m_materialDensity; }

    void translate(const std::array<T, 3>& dist) override
    {
        m_cylinder.center = vectormath::add(m_cylinder.center, dist);
    }

    void clearEnergyScored() override
    {
        for (auto& d : m_energyScored) {
            d.clear();
        }
    }

    std::array<T, 3> center() const override
    {
        return m_cylinder.center;
    }

    std::array<T, 6> AABB() const override
    {
        std::array<T, 6> aabb {
            m_cylinder.center[0] - m_cylinder.radius,
            m_cylinder.center[1] - m_cylinder.radius,
            m_cylinder.center[2] - m_cylinder.half_height,
            m_cylinder.center[0] + m_cylinder.radius,
            m_cylinder.center[1] + m_cylinder.radius,
            m_cylinder.center[2] + m_cylinder.half_height
        };
        return aabb;
    }

    WorldIntersectionResult<T> intersect(const Particle<T>& p) const noexcept override
    {
        return basicshape::cylinder::intersect(p, m_cylinder);
    }

    VisualizationIntersectionResult<T, WorldItemBase<T>> intersectVisualization(const Particle<T>& p) const noexcept override
    {
        return basicshape::cylinder::template intersectVisualization<T, WorldItemBase<T>>(p, m_cylinder);
    }

    void transport(Particle<T>& p, RandomState& state) noexcept override
    {
        bool cont = basicshape::cylinder::pointInside(p.pos, m_cylinder);
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
            const auto intLen = intersect(p).intersection; // this can not be nullopt

            if (stepLen < intLen) {
                // interaction happends
                p.translate(stepLen);
                const auto intRes = interactions::template interact<T, NMaterialShells, Lowenergycorrection>(att, p, m_material, state);
                updateAtt = intRes.particleEnergyChanged;
                cont = intRes.particleAlive;

                const auto dose_ind_f = (p.pos[2] - (m_cylinder.center[2] - m_cylinder.half_height)) * m_energyScored.size() / (m_cylinder.half_height * 2);
                const auto ind = std::clamp(static_cast<std::size_t>(dose_ind_f), std::size_t { 0 }, m_energyScored.size() - 1);
                m_energyScored[ind].scoreEnergy(intRes.energyImparted);
            } else {
                // transport to border
                p.border_translate(intLen);
                cont = false;
            }
        }
    }

    const std::vector<std::pair<T, EnergyScore<T>>> depthEnergyScored() const
    {
        std::vector<std::pair<T, EnergyScore<T>>> depth;
        depth.reserve(m_energyScored.size());

        const auto step = (2 * m_cylinder.half_height) / m_energyScored.size();
        const auto start = m_cylinder.center[2] - m_cylinder.half_height + step / 2;
        for (std::size_t i = 0; i < m_energyScored.size(); ++i) {
            depth.push_back(std::make_pair(start + step * i, m_energyScored[i]));
        }
        return depth;
    }

    const EnergyScore<T>& energyScored(std::size_t index = 0) const override
    {
        return m_energyScored[index];
    }

protected:
private:
    dxmc::basicshape::cylinder::Cylinder<T> m_cylinder;
    T m_materialDensity = 1;
    Material<T, NMaterialShells> m_material;
    std::vector<EnergyScore<T>> m_energyScored;
};

}
