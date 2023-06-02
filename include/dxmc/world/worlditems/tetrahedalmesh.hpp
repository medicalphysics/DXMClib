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
        m_collections.reserve(1);
        m_collections.emplace_back(T { 0 });
    }

    TetrahedalMesh(std::vector<Tetrahedron<T>>&& tets, const std::vector<T>& collectionDensities, const std::vector<Material2<T, NMaterialShells>>& materials, const std::vector<std::string>& collectionNames = {})
    {
        // finding mac collectionIdx and Material index
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
        for (auto d : collectionDensities)
            m_collections.emplace_back(d);

        if (collectionNames.size() == m_collections.size())
            m_collectionNames = collectionNames;
        else
            m_collectionNames.resize(m_collections.size());
        m_kdtree.setData(std::move(tets));
        m_aabb = m_kdtree.AABB();
        expandAABB();
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
            w.intersection = res.intersection;
            w.intersectionValid = true;
            w.rayOriginIsInsideItem = res.rayOriginIsInsideItem;
        }
        return w;
    }

    VisualizationIntersectionResult<T, WorldItemBase<T>> intersectVisualization(const Particle<T>& p) const override
    {
        const auto res = m_kdtree.intersect<6100>(p, m_aabb);
        VisualizationIntersectionResult<T, WorldItemBase<T>> w;
        if (res.valid()) {
            w.intersection = res.intersection;
            w.rayOriginIsInsideItem = res.rayOriginIsInsideItem;
            w.intersectionValid = true;
            w.item = this;
            w.normal = res.normal;
        }
        return w;
    }

    const DoseScore<T>& dose(std::size_t index = 0) const override
    {
        if (index < m_collections.size())
            return m_collections[index].dose;
        return m_collections[0].dose;
    }

    void clearDose() override
    {
        for (auto& c : m_collections)
            c.dose.clear();
    }

    void transport(Particle<T>& p, RandomState& state) override
    {
        TetrahedalMeshIntersectionResult<T, Tetrahedron<T>> inter = m_kdtree.intersect(p, m_aabb);
        bool updateAtt = true;

        std::uint16_t currentCollection;
        std::uint16_t currentMaterialIdx;
        AttenuationValues<T> att;
        T attSumInv;

        while (inter.valid() && inter.rayOriginIsInsideItem) {

            if (updateAtt) {
                currentCollection = inter.item->collection();
                currentMaterialIdx = inter.item->materialIndex();
                const auto& material = m_materials[currentMaterialIdx];
                att = material.attenuationValues(p.energy);
                attSumInv = 1 / (att.sum() * m_collections[currentCollection].density);
                updateAtt = false;
            }

            const auto stepLen = -std::log(state.randomUniform<T>()) * attSumInv; // cm

            if (stepLen < inter.intersection) {
                // interaction happends
                p.translate(stepLen);
                inter.intersection -= stepLen;
                const auto& material = m_materials[currentMaterialIdx];
                const auto intRes = interactions::template interact<T, NMaterialShells, LOWENERGYCORRECTION>(att, p, material, state);
                auto& dose = m_collections[currentCollection].dose;
                dose.scoreEnergy(intRes.energyImparted);
                updateAtt = intRes.particleEnergyChanged;
                if (!intRes.particleAlive)
                    inter.rayOriginIsInsideItem = false; // we exits
            } else {
                // transport to border of tetrahedron
                p.border_translate(inter.intersection);
                inter = m_kdtree.intersect(p, m_aabb);
                if (inter.valid())
                    updateAtt = currentCollection != inter.item->collection();
            }
        }
    }

protected:
    void expandAABB(T extra = 0.0001f)
    {
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] -= extra;
            m_aabb[i + 3] += extra;
        }
    }

private:
    struct Collection {
        DoseScore<T> dose;
        const T density = 0;
        Collection(T dens)
            : density(dens)
        {
        }
    };

    std::array<T, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };
    TetrahedalMeshKDTree<T> m_kdtree;
    std::vector<Collection> m_collections;
    std::vector<Material2<T, NMaterialShells>> m_materials;
    std::vector<std::string> m_collectionNames;
};
}
