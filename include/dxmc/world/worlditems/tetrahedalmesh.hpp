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

#include "dxmc/constants.hpp"
#include "dxmc/dxmcrandom.hpp"
#include "dxmc/floating.hpp"
#include "dxmc/interactions.hpp"
#include "dxmc/material/material.hpp"
#include "dxmc/material/nistmaterials.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/worlditems/tetrahedalmesh/tetrahedalmeshkdtree.hpp"
#include "dxmc/world/worlditems/tetrahedalmesh/tetrahedron.hpp"
#include "dxmc/world/worlditems/worlditembase.hpp"

#include <array>
#include <string>
#include <string_view>
#include <vector>

namespace dxmc {

template <Floating T, int NMaterialShells = 5, int LOWENERGYCORRECTION = 2>
class TetrahedalMesh final : public WorldItemBase<T> {
public:
    TetrahedalMesh()
        : WorldItemBase<T>()
    {
    }

    TetrahedalMesh(std::vector<Tetrahedron<T>>&& tets, const std::vector<T>& collectionDensities, const std::vector<Material<T, NMaterialShells>>& materials, const std::vector<std::string>& collectionNames = {})
        : WorldItemBase<T>()
    {
        // finding max collectionIdx and Material index
        const auto maxCollectionIdx = std::transform_reduce(
            std::execution::par_unseq, tets.cbegin(), tets.cend(), std::uint16_t { 0 },
            [](const auto lh, const auto rh) { return std::max(lh, rh); },
            [](const auto& t) { return t.collection(); });
        const auto maxMaterialIdx = std::transform_reduce(
            std::execution::par_unseq, tets.cbegin(), tets.cend(), std::uint16_t { 0 },
            [](const auto lh, const auto rh) { return std::max(lh, rh); },
            [](const auto& t) { return t.materialIndex(); });
        if (collectionDensities.size() <= maxCollectionIdx || materials.size() <= maxMaterialIdx)
            return;
        m_materials = materials;
        m_collections.reserve(collectionDensities.size());

        std::vector<std::atomic<T>> volumes(maxCollectionIdx + 1);
        std::for_each(std::execution::par_unseq, volumes.begin(), volumes.end(), [](auto& v) { v.store(T { 0 }); });
        std::for_each(std::execution::par_unseq, tets.cbegin(), tets.cend(), [&volumes](const auto& tet) {
            const auto idx = tet.collection();
            volumes[idx].fetch_add(tet.volume());
        });

        for (std::size_t i = 0; i <= maxCollectionIdx; ++i) {
            const auto d = collectionDensities[i];
            const auto v = volumes[i].load();
            m_collections.emplace_back(d, v);
        }

        if (collectionNames.size() == m_collections.size())
            m_collectionNames = collectionNames;
        else
            m_collectionNames.resize(m_collections.size());
        m_dose.resize(m_collections.size());
        m_kdtree.setData(std::move(tets));
        m_aabb = m_kdtree.AABB();
        m_aabb = expandAABB(m_aabb);
    }

    void translate(const std::array<T, 3>& dist) override
    {
        m_kdtree.translate(dist);
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] += dist[i];
            m_aabb[i + 3] += dist[i];
        }
    }

    std::array<T, 3> center() const override
    {
        std::array<T, 3> c;
        for (std::size_t i = 0; i < 3; ++i) {
            c[i] = (m_aabb[i] + m_aabb[i + 3]) / 2;
        }
        return c;
    }

    std::array<T, 6> AABB() const override
    {
        return m_aabb;
    }
    WorldIntersectionResult<T> intersect(const Particle<T>& p) const override
    {
        const auto res = m_kdtree.intersect(p, m_aabb);
        WorldIntersectionResult<T> w;

        if (res.valid()) {
            w.rayOriginIsInsideItem = res.rayOriginIsInsideItem();
            w.intersection = res.intersection();
            w.intersectionValid = true;
        }
        return w;
    }

    VisualizationIntersectionResult<T, WorldItemBase<T>> intersectVisualization(const Particle<T>& p) const override
    {
        const auto res = m_kdtree.intersect(p, m_aabb);
        VisualizationIntersectionResult<T, WorldItemBase<T>> w;
        if (res.valid()) {
            w.rayOriginIsInsideItem = res.rayOriginIsInsideItem();
            w.intersection = res.intersection();
            w.intersectionValid = true;
            w.item = this;
            w.normal = w.rayOriginIsInsideItem ? res.normal_exit : res.normal_enter;
            vectormath::normalize(w.normal);
        }
        return w;
    }

    const EnergyScore<T>& energyScored(std::size_t index = 0) const override
    {
        return m_collections.at(index).energyScored;
    }

    void clearEnergyScored() override
    {
        for (auto& c : m_collections)
            c.energyScored.clear();
    }

    void addEnergyScoredToDoseScore(T calibration_factor = 1) override
    {
        for (std::size_t i = 0; i < m_collections.size(); ++i) {
            const auto& c = m_collections[i];
            m_dose[i].addScoredEnergy(c.energyScored, c.volume, c.density, calibration_factor);
        }
    }

    const DoseScore<T>& doseScored(std::size_t index = 0) const override
    {
        return m_dose.at(index);
    }

    void clearDoseScored() override
    {
        for (auto& d : m_dose)
            d.clear();
    }

    void transport(Particle<T>& p, RandomState& state) override
    {
        TetrahedalMeshIntersectionResult<T, Tetrahedron<T>> inter = m_kdtree.intersect(p, m_aabb);
        bool updateAtt = true;
        std::uint16_t currentCollection;
        std::uint16_t currentMaterialIdx;
        AttenuationValues<T> att;
        T attSumInv;

        while (inter.valid() && inter.rayOriginIsInsideItem()) {
            if (updateAtt) {
                currentCollection = inter.item->collection();
                currentMaterialIdx = inter.item->materialIndex();
                const auto& material = m_materials[currentMaterialIdx];
                att = material.attenuationValues(p.energy);
                attSumInv = 1 / (att.sum() * m_collections[currentCollection].density);
                updateAtt = false;
            }

            const auto stepLen = -std::log(state.randomUniform<T>()) * attSumInv; // cm

            if (stepLen < inter.t_exit) {
                // interaction happends
                p.translate(stepLen);
                const auto& material = m_materials[currentMaterialIdx];
                const auto intRes = interactions::template interact<T, NMaterialShells, LOWENERGYCORRECTION>(att, p, material, state);
                auto& energyScored = m_collections[currentCollection].energyScored;
                energyScored.scoreEnergy(intRes.energyImparted);
                updateAtt = true;
                if (intRes.particleAlive) {
                    inter = m_kdtree.intersect(p, m_aabb);
                } else {
                    inter.item = nullptr; // we exits
                }
            } else {
                // transport to border of tetrahedron
                p.border_translate(inter.t_exit);
                inter = m_kdtree.intersect(p, m_aabb);
                if (inter.valid()) {
                    updateAtt = currentCollection != inter.item->collection();
                }
            }
        }
    }

    std::size_t numberOfCollections() const { return m_collections.size(); }
    std::size_t numberOfMaterials() const { return m_materials.size(); }
    void setMaterial(const Material<T, NMaterialShells>& material, std::size_t index)
    {
        if (index < m_materials.size())
            m_materials[index] = material;
    }

protected:
    [[nodiscard]] static std::array<T, 6> expandAABB(std::array<T, 6> aabb)
    {
        for (std::size_t i = 0; i < 3; ++i) {
            aabb[i] -= GEOMETRIC_ERROR<T>();
            aabb[i + 3] += GEOMETRIC_ERROR<T>();
        }
        return aabb;
    }

private:
    struct Collection {
        EnergyScore<T> energyScored;
        const T density = 0;
        const T volume = 0;
        Collection(T dens, T volume)
            : density(dens)
            , volume(volume)
        {
        }
    };

    std::array<T, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };
    TetrahedalMeshKDTree<T> m_kdtree;
    std::vector<Collection> m_collections;
    std::vector<DoseScore<T>> m_dose;
    std::vector<Material<T, NMaterialShells>> m_materials;
    std::vector<std::string> m_collectionNames;
};
}
