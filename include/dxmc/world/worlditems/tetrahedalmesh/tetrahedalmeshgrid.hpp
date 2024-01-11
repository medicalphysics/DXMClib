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

#include "dxmc/floating.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/world/basicshapes/aabb.hpp"
#include "dxmc/world/kdtreeintersectionresult.hpp"
#include "dxmc/world/worlditems/tetrahedalmesh/tetrahedron.hpp"

#include <algorithm>
#include <execution>
#include <limits>
#include <ranges>
#include <vector>

namespace dxmc {

template <Floating T>
class TetrahedalMeshGrid {
public:
    TetrahedalMeshGrid() { }
    TetrahedalMeshGrid(const std::vector<Tetrahedron<T>>& tets, int depth = 8)
    {
        setData(tets, depth);
    }
    TetrahedalMeshGrid(const std::vector<Tetrahedron<T>>& tets, int dx, int dy, int dz)
    {
        setData(tets, dx, dy, dz);
    }
    TetrahedalMeshGrid(const std::vector<Tetrahedron<T>>& tets, const std::array<int, 3> depth)
    {
        setData(tets, depth);
    }

    void setData(const std::vector<Tetrahedron<T>>& tets, int depth = 8)
    {
        std::array<int, 3> d = { depth, depth, depth };
        setData(tets, d);
    }

    void setData(const std::vector<Tetrahedron<T>>& tets, int dx, int dy, int dz)
    {
        std::array<int, 3> d = { dx, dy, dz };
        setData(tets, d);
    }

    void setData(const std::vector<Tetrahedron<T>>& tets, const std::array<int, 3>& depth)
    {
        for (int i = 0; i < 3; ++i)
            m_N[i] = std::clamp(depth[i], 1, 1000);

        m_tets = tets;
        // sorting tets by diagonal
        std::sort(std::execution::par_unseq, m_tets.begin(), m_tets.end(), [](const auto& lh, const auto& rh) {
            constexpr std::array<T, 3> n = { 1, 1, 1 };
            return vectormath::dot(lh.center(), n) < vectormath::dot(rh.center(), n);
        });

        calculateAABB();
        assignGrid();
    }

    const std::vector<Tetrahedron<T>>& tetrahedrons() const
    {
        return m_tets;
    }

    void resample(const std::array<int, 3>& depth)
    {
        for (int i = 0; i < 3; ++i)
            m_N[i] = std::clamp(depth[i], 1, 1000);
        assignGrid();
    }

    const std::array<T, 6>& AABB() const
    {
        return m_aabb;
    }

    template <std::regular_invocable<Tetrahedron<T>> F>
    void apply(const F& func, auto policy = std::execution::par_unseq)
    {
        std::for_each(policy, m_tets.begin(), m_tets.end(), func);
    }

    void translate(const std::array<T, 3>& dist)
    {
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] += dist[i];
            m_aabb[i + 3] += dist[i];
        }
        std::for_each(std::execution::par_unseq, m_tets.begin(), m_tets.end(), [&](auto& tri) {
            tri.translate(dist);
        });
    }

    template <std::uint16_t COLLECTION = 65535>
    KDTreeIntersectionResult<T, const Tetrahedron<T>> intersect(const Particle<T>& particle, const std::array<T, 6>& aabb) const
    {
        const auto inter = basicshape::AABB::intersectForwardInterval<T, false>(particle, aabb);
        return inter ? intersect<COLLECTION>(particle, *inter) : KDTreeIntersectionResult<T, const Tetrahedron<T>> {};
    }

    template <std::uint16_t COLLECTION = 65535>
    KDTreeIntersectionResult<T, const Tetrahedron<T>> intersectExit(const Particle<T>& particle, const std::array<T, 6>& aabb) const
    {
        const auto inter = basicshape::AABB::intersectForwardInterval<T, false>(particle, aabb);
        if (inter) {
            auto p_reverse = particle;
            const auto& t = *inter;
            p_reverse.translate(t[1]);
            p_reverse.dir = vectormath::scale(p_reverse.dir, T { -1 });
            auto res = intersect<COLLECTION>(p_reverse, *inter);
            res.intersection = t[1] - res.intersection;
            return res;
        }
        return KDTreeIntersectionResult<T, const Tetrahedron<T>> {};
    }

    const Tetrahedron<T>* pointInside(const std::array<T, 3>& pos) const
    {
        const auto idx = getIndices<true>(pos);
        const auto idx_flat = idx[0] + idx[1] * m_N[0] + idx[2] * m_N[0] * m_N[1];
        for (const auto tet_idx : m_grid[idx_flat]) {
            const auto& tet = m_tets[tet_idx];
            if (tet.pointInside(pos)) {
                return &tet;
            }
        }
        return nullptr;
    }

protected:
    template <bool CHECK_BOUNDS = false>
    std::array<int, 3> getIndices(const std::array<T, 3>& pos) const
    {
        if constexpr (CHECK_BOUNDS) {
            std::array<int, 3> res = {
                std::clamp(static_cast<int>((pos[0] - m_aabb[0]) / m_spacing[0]), 0, m_N[0] - 1),
                std::clamp(static_cast<int>((pos[1] - m_aabb[1]) / m_spacing[1]), 0, m_N[1] - 1),
                std::clamp(static_cast<int>((pos[2] - m_aabb[2]) / m_spacing[2]), 0, m_N[2] - 1)
            };
            return res;
        } else {
            std::array<int, 3> res = {
                static_cast<int>((pos[0] - m_aabb[0]) / m_spacing[0]),
                static_cast<int>((pos[1] - m_aabb[1]) / m_spacing[1]),
                static_cast<int>((pos[2] - m_aabb[2]) / m_spacing[2])
            };
            return res;
        }
    }

    static inline int argmin3(const std::array<T, 3>& a)
    {
        return a[0] < a[2] ? a[0] < a[1] ? 0 : 1 : a[2] < a[1] ? 2
                                                               : 1;
    }

    template <std::uint16_t COLLECTION = 65535>
    KDTreeIntersectionResult<T, const Tetrahedron<T>> intersect(const Particle<T>& p, const std::array<T, 2>& t) const
    {
        auto idx = getIndices<true>(vectormath::add(p.pos, vectormath::scale(p.dir, t[0])));
        const std::array<int, 3> step = {
            p.dir[0] < 0 ? -1 : 1,
            p.dir[1] < 0 ? -1 : 1,
            p.dir[2] < 0 ? -1 : 1
        };
        const std::array<T, 3> delta = {
            m_spacing[0] / std::abs(p.dir[0]),
            m_spacing[1] / std::abs(p.dir[1]),
            m_spacing[2] / std::abs(p.dir[2])
        };

        std::array<T, 3> tmax;
        for (int i = 0; i < 3; ++i) {
            tmax[i] = step[i] > 0 ? (m_aabb[i] + (idx[i] + 1) * m_spacing[i] - p.pos[i]) / p.dir[i] : (m_aabb[i] + idx[i] * m_spacing[i] - p.pos[i]) / p.dir[i];
        };

        int dimension = argmin3(tmax);
        KDTreeIntersectionResult<T, const Tetrahedron<T>> res;
        res.intersection = std::numeric_limits<T>::max();
        bool cont = true;
        while (cont) {
            // we have a valid voxel, check intersections
            // const auto voxel_ind = idx[0] + (idx[1] + idx[2] * m_N[1]) * m_N[0];
            const auto voxel_ind = idx[0] + idx[1] * m_N[0] + idx[2] * m_N[0] * m_N[1];
            for (const auto& tetIdx : m_grid[voxel_ind]) {
                const auto& tet = m_tets[tetIdx];
                if constexpr (COLLECTION == 65535) {
                    const auto res_cand = tet.intersect(p);
                    if (res_cand.valid() && res_cand.intersection <= tmax[dimension] && res_cand.intersection < res.intersection) {
                        res.intersection = res_cand.intersection;
                        res.rayOriginIsInsideItem = res_cand.rayOriginIsInsideItem;
                        res.item = &tet;
                        if (res.rayOriginIsInsideItem) // early exit if we are inside tet
                            return res;
                    }
                } else {
                    if (tet.collection == COLLECTION) {
                        const auto res_cand = tet.intersect(p);
                        if (res_cand.valid() && res_cand.intersection <= tmax[dimension] && res_cand.intersection < res.intersection) {
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
            if (!res.valid() && 0 <= idx[dimension] && idx[dimension] < m_N[dimension]) {
                tmax[dimension] += delta[dimension];
                dimension = argmin3(tmax);
            } else {
                cont = false;
            }
        }
        return res;
    }

    void calculateAABB()
    {
        m_aabb = { std::numeric_limits<T>::max(),
            std::numeric_limits<T>::max(),
            std::numeric_limits<T>::max(),
            std::numeric_limits<T>::lowest(),
            std::numeric_limits<T>::lowest(),
            std::numeric_limits<T>::lowest() };
        for (const auto tet : m_tets) {
            for (const auto& v : tet.vertices()) {
                for (std::size_t i = 0; i < 3; ++i) {
                    m_aabb[i] = std::min(m_aabb[i], v[i]);
                    m_aabb[i + 3] = std::max(m_aabb[i + 3], v[i]);
                }
            }
        }
        for (std::size_t i = 0; i < 3; ++i)
            m_spacing[i] = (m_aabb[i + 3] - m_aabb[i]) / m_N[i];
    }

    void assignGrid()
    {
        const auto [start, stop] = vectormath::splice(m_aabb);
        const auto size = std::reduce(m_N.cbegin(), m_N.cend(), 1, std::multiplies {});
        m_grid.resize(size);

        const std::array<int, 3> N = { m_N[0] - 1, m_N[1] - 1, m_N[2] - 1 };
        auto caster = [N](const std::array<T, 3>& v) -> std::array<int, 3> {
            std::array<int, 3> vi = {
                std::clamp(static_cast<int>(v[0]), 0, N[0]),
                std::clamp(static_cast<int>(v[1]), 0, N[1]),
                std::clamp(static_cast<int>(v[2]), 0, N[2])
            };
            return vi;
        };
        const std::array<T, 3> inv_spacing = { 1 / m_spacing[0], 1 / m_spacing[1], 1 / m_spacing[2] };
        for (std::size_t i = 0; i < m_tets.size(); ++i) {
            const auto tet_aabb = m_tets[i].AABB();
            const auto [tet_start, tet_stop] = vectormath::splice(tet_aabb);
            const auto start_ind = caster(vectormath::scale(vectormath::subtract(tet_start, start), inv_spacing));
            const auto stop_ind = caster(vectormath::scale(vectormath::subtract(tet_stop, start), inv_spacing));
            for (int z = start_ind[2]; z <= stop_ind[2]; ++z)
                for (int y = start_ind[1]; y <= stop_ind[1]; ++y)
                    for (int x = start_ind[0]; x <= stop_ind[0]; ++x) {
                        const int idx = x + m_N[0] * y + m_N[0] * m_N[1] * z;
                        m_grid[idx].push_back(i);
                    }
        }
        std::for_each(std::execution::par_unseq, m_grid.begin(), m_grid.end(), [](auto& v) { v.shrink_to_fit(); });
    }

private:
    std::array<T, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };
    std::array<T, 6> m_spacing = { 1, 1, 1 };
    std::vector<std::vector<std::size_t>> m_grid;
    std::vector<Tetrahedron<T>> m_tets;
    std::array<int, 3> m_N = { 8, 8, 8 };
    // int m_N = 32;
};
}