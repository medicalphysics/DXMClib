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
#include "dxmc/interpolation.hpp"
#include "dxmc/material/material.hpp"
#include "dxmc/material/nistmaterials.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/worlditems/tetrahedalmesh/tetrahedalmeshgrid.hpp"
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

    TetrahedalMesh(const std::vector<Tetrahedron<T>>& tets, const std::vector<T>& collectionDensities, const std::vector<Material<T, NMaterialShells>>& materials, const std::vector<std::string>& collectionNames = {}, std::array<int, 3> depth = { 8, 8, 8 })
        : WorldItemBase<T>()
    {
        setData(tets, collectionDensities, materials, collectionNames, depth);
    }

    TetrahedalMesh(const std::vector<Tetrahedron<T>>& tets, const std::vector<T>& collectionDensities, const std::vector<Material<T, NMaterialShells>>& materials, const std::vector<std::string>& collectionNames = {}, int max_depth = 8)
        : WorldItemBase<T>()
    {
        setData(tets, collectionDensities, materials, collectionNames, { max_depth, max_depth, max_depth });
    }

    bool setData(const std::vector<Tetrahedron<T>>& tets, const std::vector<T>& collectionDensities, const std::vector<Material<T, NMaterialShells>>& materials, const std::vector<std::string>& collectionNames = {}, std::array<int, 3> depth = { 8, 8, 8 })
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
            return false;
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

        m_grid.setData(tets, depth);
        generateWoodcockStepTable();
        return true;
    }

    void translate(const std::array<T, 3>& dist) override
    {
        m_grid.translate(dist);
    }

    std::array<T, 3> center() const override
    {
        const auto& aabb = m_grid.AABB();
        const auto [low, high] = vectormath::splice(aabb);
        auto c = vectormath::add(low, high);
        return vectormath::scale(c, T { 0.5 });
    }

    std::array<T, 6> AABB() const override
    {
        return m_grid.AABB();
    }

    WorldIntersectionResult<T> intersect(const Particle<T>& p) const override
    {
        WorldIntersectionResult<T> res;
        if (const auto kres = m_grid.intersect(p); kres.valid()) {
            res.intersection = kres.intersection;
            res.rayOriginIsInsideItem = kres.rayOriginIsInsideItem;
            res.intersectionValid = kres.valid();
        }
        return res;
    }

    VisualizationIntersectionResult<T, WorldItemBase<T>> intersectVisualization(const Particle<T>& p) const override
    {
        VisualizationIntersectionResult<T, WorldItemBase<T>> res;
        if (const auto kres = m_grid.intersect(p); kres.valid()) {
            res.intersection = kres.intersection;
            res.rayOriginIsInsideItem = kres.rayOriginIsInsideItem;
            res.intersectionValid = kres.valid();
            res.item = this;
            const auto collection = kres.item->collection();
            res.value = kres.item->doseScored().dose();
            const auto hit_pos = vectormath::add(p.pos, vectormath::scale(p.dir, kres.intersection));
            res.normal = kres.item->normal(hit_pos);
        }
        return res;
    }

    const EnergyScore<T>& energyScored(std::size_t index = 0) const override
    {
        const auto tets = m_grid.tetrahedrons();
        return tets.at(index).energyImparted();
    }

    void clearEnergyScored() override
    {
        m_grid.clearEnergyScored();
    }

    void addEnergyScoredToDoseScore(T calibration_factor = 1) override
    {
        auto& tets = m_grid.tetrahedrons();
        std::for_each(std::execution::par_unseq, tets.begin(), tets.end(), [=](auto& tet) {
            const auto cidx = tet.collection();
            const auto dens = this->m_collections[cidx].density;
            tet.addEnergyScoredToDoseScore(dens, calibration_factor);
        });
    }

    const DoseScore<T>& doseScored(std::size_t index = 0) const override
    {
        const auto tets = m_grid.tetrahedrons();
        return tets.at(index).doseScored();
    }

    void clearDoseScored() override
    {
        // for (auto& d : m_dose)
        //     d.clear();
        m_grid.clearDoseScored();
    }

    void transport(Particle<T>& p, RandomState& state) override
    {
        bool still_inside = true;
        T attMaxInv;
        bool updateAtt = true;
        do {
            if (updateAtt) {
                attMaxInv = 1 / interpolate(m_woodcockStepTableLin, p.energy);
                bool updateAtt = false;
            }

            // making interaction step
            const auto steplen = -std::log(state.randomUniform<T>()) * attMaxInv;
            p.translate(steplen);

            // finding current tet
            const auto currentTet = m_grid.pointInside(p.pos);

            if (!currentTet) {
                // we are outside item, backtrack and return
                const Particle<T> pback = { .pos = p.pos, .dir = vectormath::scale(p.dir, T { -1 }) };
                const auto inter_back = intersect(pback);
                if (inter_back.valid()) {
                    p.border_translate(-inter_back.intersection);
                }
                return;
            }

            // is interaction virtual?
            const auto materialIdx = currentTet->materialIndex();
            const auto collectionIdx = currentTet->collection();
            const auto attenuation = m_materials[materialIdx].attenuationValues(p.energy);
            const auto attSum = attenuation.sum() * m_collections[collectionIdx].density;
            if (state.randomUniform<T>() < attSum * attMaxInv) {
                // we have a real interaction
                const auto intRes = interactions::template interact<T, NMaterialShells, LOWENERGYCORRECTION>(attenuation, p, m_materials[materialIdx], state);
                currentTet->scoreEnergy(intRes.energyImparted);
                still_inside = intRes.particleAlive;
                updateAtt = intRes.particleEnergyChanged;
            }
        } while (still_inside);
    }

    std::size_t numberOfCollections() const { return m_collections.size(); }
    std::size_t numberOfMaterials() const { return m_materials.size(); }
    std::size_t numberOfTetrahedra() const { return m_grid.tetrahedrons().size(); }
    const std::vector<Tetrahedron<T>>& tetrahedrons() const
    {
        return m_grid.tetrahedrons();
    }

    void setMaterial(const Material<T, NMaterialShells>& material, std::size_t index)
    {
        if (index < m_materials.size())
            m_materials[index] = material;
    }

protected:
    void generateWoodcockStepTable()
    {
        std::vector<T> energy;
        {
            T e = std::log(MIN_ENERGY<T>());
            const T emax = std::log(MAX_ENERGY<T>());
            const T estep = (emax - e) / 10;
            while (e <= emax) {
                energy.push_back(e);
                e += estep;
            }
        }
        // adding edges;
        for (const auto& mat : m_materials) {
            for (std::size_t i = 0; i < mat.numberOfShells(); ++i) {
                const auto& shell = mat.shell(i);
                const auto e = shell.bindingEnergy + T { 0.01 };
                if (e > MIN_ENERGY<T>()) {
                    energy.push_back(std::log(e));
                }
            }
        }
        std::sort(energy.begin(), energy.end());
        auto remove = std::unique(energy.begin(), energy.end());
        energy.erase(remove, energy.end());
        std::transform(std::execution::par_unseq, energy.cbegin(), energy.cend(), energy.begin(), [](const auto e) { return std::exp(e); });

        // finding max density for each material;
        std::vector<T> dens(m_materials.size(), T { 0 });
        for (const auto& tet : m_grid.tetrahedrons()) {
            const auto i = tet.materialIndex();
            const auto c = tet.collection();
            dens[i] = std::max(m_collections[c].density, dens[i]);
        }

        // finding max attenuation for each energy
        std::vector<T> att(energy.size(), T { 0 });
        for (std::size_t mIdx = 0; mIdx < m_materials.size(); ++mIdx) {
            const auto& mat = m_materials[mIdx];
            const auto d = dens[mIdx];
            for (std::size_t i = 0; i < energy.size(); ++i) {
                const auto aval = mat.attenuationValues(energy[i]);
                att[i] = std::max(aval.sum() * d, att[i]);
            }
        }
        std::vector<std::pair<T, T>> data(energy.size());
        std::transform(std::execution::par_unseq, energy.cbegin(), energy.cend(), att.cbegin(), data.begin(), [](const auto e, const auto a) {
            return std::make_pair(e, a);
        });

        m_woodcockStepTableLin = data;
    }

private:
    struct Collection {
        const T density = 0;
        const T volume = 0;
        Collection(T dens, T volume)
            : density(dens)
            , volume(volume)
        {
        }
    };

    TetrahedalMeshGrid<T> m_grid;
    std::vector<std::pair<T, T>> m_woodcockStepTableLin;
    std::vector<Collection> m_collections;
    std::vector<Material<T, NMaterialShells>> m_materials;
    std::vector<std::string> m_collectionNames;
};
}
