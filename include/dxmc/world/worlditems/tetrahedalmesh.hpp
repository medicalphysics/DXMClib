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
        m_dose.resize(m_collections.size());

        // sorting tets
        m_tets = tets;
        std::sort(std::execution::par_unseq, m_tets.begin(), m_tets.end(), [](const auto& lh, const auto& rh) {
            constexpr std::array<T, 3> n = { 1, 1, 1 };
            return vectormath::dot(lh.center(), n) < vectormath::dot(rh.center(), n);
        });
        calculateAABB();
        // setting up grid
        for (int i = 0; i < 3; ++i)
            m_gridDimensions[i] = std::clamp(depth[i], 1, 1000);
        assignGrid();
        generateWoodcockStepTable();
        return true;
    }

    void translate(const std::array<T, 3>& dist) override
    {

        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] += dist[i];
            m_aabb[i + 3] += dist[i];
        }
        std::for_each(std::execution::par_unseq, m_tets.begin(), m_tets.end(), [&](auto& tri) {
            tri.translate(dist);
        });
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

    WorldIntersectionResult<T> intersect(const Particle<T>& particle) const override
    {
        const auto inter = basicshape::AABB::intersectForwardInterval(particle, m_aabb);
        WorldIntersectionResult<T> res;
        if (inter) {
            const auto kres = intersect(particle, *inter);
            res.intersection = kres.intersection;
            res.rayOriginIsInsideItem = kres.rayOriginIsInsideItem;
            res.intersectionValid = kres.valid();
        }
        return res;
    }

    VisualizationIntersectionResult<T, WorldItemBase<T>> intersectVisualization(const Particle<T>& p) const override
    {
        VisualizationIntersectionResult<T, WorldItemBase<T>> w;
        const auto inter = basicshape::AABB::intersectForwardInterval(p, m_aabb);
        if (inter) {
            if (const auto res = intersect(p, *inter); res.valid()) {
                w.rayOriginIsInsideItem = res.rayOriginIsInsideItem;
                w.intersection = res.intersection;
                w.intersectionValid = true;
                w.item = this;
                const auto collection = res.item->collection();
                w.value = m_dose[collection].dose();
                const auto hit_pos = vectormath::add(p.pos, vectormath::scale(p.dir, res.intersection));
                w.normal = res.item->normal(hit_pos);
            }
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

    const std::vector<DoseScore<T>>& getDoseScores() const
    {
        return m_dose;
    }

    void clearDoseScored() override
    {
        for (auto& d : m_dose)
            d.clear();
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

            if (!basicshape::AABB::pointInside(p.pos, m_aabb)) {
                // we have left AABB, backtrack to right after mesh
                const Particle<T> pback = { .pos = p.pos, .dir = vectormath::scale(p.dir, T { -1 }) };
                const auto inter_back = intersect(pback);
                if (inter_back.valid())
                    p.border_translate(-inter_back.intersection);
                return;
            }

            // finding current tet
            const auto idx = getGridIndices<true>(p.pos);
            const auto idx_flat = getGridIndicesFlat(idx);
            const Tetrahedron<T>* tet = nullptr;
            {
                auto it = m_gridIndices[idx_flat].cbegin();
                auto end = m_gridIndices[idx_flat].cend();
                while (it != end) {
                    const auto& ctet = m_tets[*it];
                    if (ctet.pointInside(p.pos)) {
                        tet = &ctet;
                        it = end;
                    } else {
                        it++;
                    }
                }
            }
            if (!tet) {
                // we are outside item, backtrack and return
                const Particle<T> pback = { .pos = p.pos, .dir = vectormath::scale(p.dir, T { -1 }) };
                const auto inter_back = intersect(pback);
                if (inter_back.valid())
                    p.border_translate(-inter_back.intersection);
                return;
            }

            // is interaction virtual?
            const auto materialIdx = tet->materialIndex();
            const auto collectionIdx = tet->collection();
            const auto attenuation = m_materials[materialIdx].attenuationValues(p.energy);
            const auto attSum = attenuation.sum() * m_collections[collectionIdx].density;
            if (state.randomUniform<T>() < attSum * attMaxInv) {
                // we have a real interaction
                const auto intRes = interactions::template interact<T, NMaterialShells, LOWENERGYCORRECTION>(attenuation, p, m_materials[materialIdx], state);
                m_collections[collectionIdx].energyScored.scoreEnergy(intRes.energyImparted);
                still_inside = intRes.particleAlive;
                updateAtt = intRes.particleEnergyChanged;
            }

        } while (still_inside);
    }

    std::size_t numberOfCollections() const { return m_collections.size(); }
    std::size_t numberOfMaterials() const { return m_materials.size(); }
    void setMaterial(const Material<T, NMaterialShells>& material, std::size_t index)
    {
        if (index < m_materials.size())
            m_materials[index] = material;
    }

protected:
    void calculateAABB()
    {
        m_aabb = { std::numeric_limits<T>::max(),
            std::numeric_limits<T>::max(),
            std::numeric_limits<T>::max(),
            std::numeric_limits<T>::lowest(),
            std::numeric_limits<T>::lowest(),
            std::numeric_limits<T>::lowest() };
        for (const auto& tet : m_tets) {
            for (const auto& v : tet.vertices()) {
                for (std::size_t i = 0; i < 3; ++i) {
                    m_aabb[i] = std::min(m_aabb[i], v[i]);
                    m_aabb[i + 3] = std::max(m_aabb[i + 3], v[i]);
                }
            }
        }

        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] -= GEOMETRIC_ERROR<T>();
            m_aabb[i + 3] += GEOMETRIC_ERROR<T>();
        }

        for (std::size_t i = 0; i < 3; ++i)
            m_gridSpacing[i] = (m_aabb[i + 3] - m_aabb[i]) / m_gridDimensions[i];
    }

    void assignGrid()
    {
        const auto [start, stop] = vectormath::splice(m_aabb);
        const auto size = std::reduce(m_gridDimensions.cbegin(), m_gridDimensions.cend(), 1, std::multiplies {});
        m_gridIndices.resize(size);

        const std::array<int, 3> N = { m_gridDimensions[0] - 1, m_gridDimensions[1] - 1, m_gridDimensions[2] - 1 };
        auto caster = [N](const std::array<T, 3>& v) -> std::array<int, 3> {
            std::array<int, 3> vi = {
                std::clamp(static_cast<int>(v[0]), 0, N[0]),
                std::clamp(static_cast<int>(v[1]), 0, N[1]),
                std::clamp(static_cast<int>(v[2]), 0, N[2])
            };
            return vi;
        };
        const std::array<T, 3> inv_spacing = { 1 / m_gridSpacing[0], 1 / m_gridSpacing[1], 1 / m_gridSpacing[2] };
        for (std::size_t i = 0; i < m_tets.size(); ++i) {
            const auto tet_aabb = m_tets[i].AABB();
            const auto [tet_start, tet_stop] = vectormath::splice(tet_aabb);
            const auto start_ind = caster(vectormath::scale(vectormath::subtract(tet_start, start), inv_spacing));
            const auto stop_ind = caster(vectormath::scale(vectormath::subtract(tet_stop, start), inv_spacing));
            for (int z = start_ind[2]; z <= stop_ind[2]; ++z)
                for (int y = start_ind[1]; y <= stop_ind[1]; ++y)
                    for (int x = start_ind[0]; x <= stop_ind[0]; ++x) {
                        const int idx = x + m_gridDimensions[0] * y + m_gridDimensions[0] * m_gridDimensions[1] * z;
                        m_gridIndices[idx].push_back(i);
                    }
        }
        std::for_each(std::execution::par_unseq, m_gridIndices.begin(), m_gridIndices.end(), [](auto& v) { v.shrink_to_fit(); });
    }

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
        for (const auto& tet : m_tets) {
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

    template <bool CHECK_BOUNDS = false>
    std::array<int, 3> getGridIndices(const std::array<T, 3>& pos) const
    {
        if constexpr (CHECK_BOUNDS) {
            std::array<int, 3> res = {
                std::clamp(static_cast<int>((pos[0] - m_aabb[0]) / m_gridSpacing[0]), 0, m_gridDimensions[0] - 1),
                std::clamp(static_cast<int>((pos[1] - m_aabb[1]) / m_gridSpacing[1]), 0, m_gridDimensions[1] - 1),
                std::clamp(static_cast<int>((pos[2] - m_aabb[2]) / m_gridSpacing[2]), 0, m_gridDimensions[2] - 1)
            };
            return res;
        } else {
            std::array<int, 3> res = {
                static_cast<int>((pos[0] - m_aabb[0]) / m_gridSpacing[0]),
                static_cast<int>((pos[1] - m_aabb[1]) / m_gridSpacing[1]),
                static_cast<int>((pos[2] - m_aabb[2]) / m_gridSpacing[2])
            };
            return res;
        }
    }

    int getGridIndicesFlat(const std::array<int, 3>& idx) const
    {
        return idx[0] + m_gridDimensions[0] * (idx[1] + m_gridDimensions[1] * idx[2]);
    }

    static inline int argmin3(const std::array<T, 3>& a)
    {
        return a[0] < a[2] ? a[0] < a[1] ? 0 : 1 : a[2] < a[1] ? 2
                                                               : 1;
    }

    template <std::uint16_t COLLECTION = 65535>
    KDTreeIntersectionResult<T, const Tetrahedron<T>> intersect(const Particle<T>& p, const std::array<T, 2>& t) const
    {
        auto idx = getGridIndices<true>(vectormath::add(p.pos, vectormath::scale(p.dir, t[0])));
        const std::array<int, 3> step = {
            p.dir[0] < 0 ? -1 : 1,
            p.dir[1] < 0 ? -1 : 1,
            p.dir[2] < 0 ? -1 : 1
        };
        const std::array<T, 3> delta = {
            m_gridSpacing[0] / std::abs(p.dir[0]),
            m_gridSpacing[1] / std::abs(p.dir[1]),
            m_gridSpacing[2] / std::abs(p.dir[2])
        };

        std::array<T, 3> tmax;
        for (int i = 0; i < 3; ++i) {
            tmax[i] = step[i] > 0 ? (m_aabb[i] + (idx[i] + 1) * m_gridSpacing[i] - p.pos[i]) / p.dir[i] : (m_aabb[i] + idx[i] * m_gridSpacing[i] - p.pos[i]) / p.dir[i];
        };

        int dimension = argmin3(tmax);
        KDTreeIntersectionResult<T, const Tetrahedron<T>> res;
        res.intersection = std::numeric_limits<T>::max();
        bool cont = true;
        while (cont) {
            // we have a valid voxel, check intersections
            // const auto voxel_ind = idx[0] + (idx[1] + idx[2] * m_N[1]) * m_N[0];
            const auto voxel_ind = getGridIndicesFlat(idx);
            for (const auto tetIdx : m_gridIndices[voxel_ind]) {
                const auto& tet = m_tets[tetIdx];
                if constexpr (COLLECTION == 65535) {
                    const auto res_cand = tet.intersect(p);
                    if (res_cand.valid() && (res_cand.intersection <= tmax[dimension]) && (res_cand.intersection < res.intersection)) {
                        res.intersection = res_cand.intersection;
                        res.rayOriginIsInsideItem = res_cand.rayOriginIsInsideItem;
                        res.item = &tet;
                        if (res.rayOriginIsInsideItem) // early exit if we are inside tet
                            return res;
                    }
                } else {
                    if (tet.collection == COLLECTION) {
                        const auto res_cand = tet.intersect(p);
                        if (res_cand.valid() && (res_cand.intersection <= tmax[dimension]) && (res_cand.intersection < res.intersection)) {
                            res.intersection = res_cand.intersection;
                            res.rayOriginIsInsideItem = res_cand.rayOriginIsInsideItem;
                            res.item = &tet;
                            if (res.rayOriginIsInsideItem) // early exit if we are inside tet
                                return res;
                        }
                    }
                }
            }
            idx[dimension] += step[dimension];
            if (!res.valid() && (0 <= idx[dimension]) && (idx[dimension] < m_gridDimensions[dimension])) {
                tmax[dimension] += delta[dimension];
                dimension = argmin3(tmax);
            } else {
                cont = false;
            }
        }
        return res;
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
    std::array<int, 3> m_gridDimensions = { 8, 8, 8 };
    std::array<T, 3> m_gridSpacing = { 1, 1, 1 };
    std::vector<std::vector<std::size_t>> m_gridIndices;
    std::vector<Tetrahedron<T>> m_tets;
    std::vector<std::pair<T, T>> m_woodcockStepTableLin;
    std::vector<Collection> m_collections;
    std::vector<DoseScore<T>> m_dose;
    std::vector<Material<T, NMaterialShells>> m_materials;
    std::vector<std::string> m_collectionNames;
};
}
