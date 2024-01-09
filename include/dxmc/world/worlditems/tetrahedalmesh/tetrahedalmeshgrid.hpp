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
#include "dxmc/world/worlditems/tetrahedalmesh/tetrahedron.hpp"

#include <algorithm>
#include <execution>
#include <limits>
#include <ranges>
#include <vector>

namespace dxmc {

template <Floating T, int N = 16>
    requires N > 1 && N <= 1000;
class TetrahedalMeshGrid {
public:
    TetrahedalMeshGrid() { }
    TetrahedalMeshGrid(const std::vector<Tetrahedron<T>>& tets)
    {
    }
    void setData(const std::vector<Tetrahedron<T>>& tets)
    {
        m_tets = tets;
        // sorting tets by diagonal
        std::sort(std::execution::par_unseq, m_tets.begin(), m_tets.end(), [](const auto& lh, const auto& rh) {
            constexpr std::array<T, 3> n = { 1, 1, 1 };
            return vectormath::dot(lh.center(), n) < vectormath::dot(rh.center(), n);
        });
        calculateAABB();
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
    TetrahedalMeshIntersectionResult<T, Tetrahedron<T>> intersect(const Particle<T>& particle, const std::array<T, 6>& aabb) const
    {
        const auto inter = basicshape::AABB::intersectForwardInterval<T, false>(particle, aabb);
        return inter ? intersect<COLLECTION>(particle, *inter) : TetrahedalMeshIntersectionResult<T, Tetrahedron<T>> {};
    }

protected:
    template <bool CHECK_BOUNDS = false>
    std::array<int, 3> getIndices(const std::array<T, 3>& pos)
    {
        if constexpr (CHECK_BOUNDS) {
            std::array<int, 3> res = {
                static_cast<int>((pos[0] - m_aabb[0]) / m_spacing[0]),
                static_cast<int>((pos[1] - m_aabb[1]) / m_spacing[1]),
                static_cast<int>((pos[2] - m_aabb[2]) / m_spacing[2])
            };
            return res;
        } else {
            std::array<int, 3> res = {
                std::clamp(static_cast<int>((pos[0] - m_aabb[0]) / m_spacing[0]), 0, N - 1),
                std::clamp(static_cast<int>((pos[1] - m_aabb[1]) / m_spacing[1]), 0, N - 1),
                std::clamp(static_cast<int>((pos[2] - m_aabb[2]) / m_spacing[2]), 0, N - 1)
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
    TetrahedalMeshIntersectionResult<T, Tetrahedron<T>> intersect(const Particle<T>& particle, const std::array<T, 2>& t) const
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
            tmax[i] = step[i] > 0 ? (aabb[i] + (idx[i] + 1) * m_spacing[i] - p.pos[i]) / p.dir[i] : (aabb[i] + idx[i] * m_spacing[i] - p.pos[i]) / p.dir[i];
        };

        int dimension = argmin3<int, T>(tmax);
        TetrahedalMeshIntersectionResult<T, Tetrahedron<T>> res;
        res.intersection = std::numeric_limits<T>::max();
        bool cont = true;
        while (cont) {
            // we have a valid voxel, check intersections
            const auto voxel_ind = idx[0] + (idx[1] + idx[2] * N) * N;
            for (const auto& tetIdx : m_grid[voxel_ind]) {
                if constexpr (COLLECTION == 65535) {
                    const auto res_cand = m_tets[tetIdx].intersect(p);
                    if (res.valid() && res_cand.intersection < res.intersection && res_cand.intersection < tmax[dimension]) {
                        std::swap(res, res_cand);
                    }
                } else {
                    if (m_tets[tetIdx].collection == COLLECTION) {
                        const auto res_cand = m_tets[tetIdx].intersect(p);
                        if (res.valid() && res_cand.intersection < res.intersection && res_cand.intersection < tmax[dimension]) {
                            std::swap(res, res_cand);
                        }
                    }
                }
            }
            idx[dimension] += step[dimension];
            if (!res.valid() && 0 <= idx[dimension] && idx[dimension] < N) {
                tmax[dimension] += delta[dimension];
                dimension = argmin3<int, T>(tmax);
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
            m_spacing[i] = (m_aabb[i + 3] - m_aabb[i]) / N;
    }

    void assignGrid()
    {
        const auto [start, stop] = vectormath::splice(m_aabb);
        // Inefficient loop over grid, better to loop over tetrahedrons?
        for (int z = 0; z < N; ++z) {
            for (int y = 0; y < N; ++y) {
                for (int x = 0; x < N; ++x) {
                    const auto idx = x + y * N + z * N * N;
                    const std::array<T, 6> aabb = {
                        start[0] + x * m_spacing[0],
                        start[1] + y * m_spacing[1],
                        start[2] + z * m_spacing[2],
                        start[0] + (x + 1) * m_spacing[0],
                        start[1] + (y + 1) * m_spacing[1],
                        start[2] + (z + 1) * m_spacing[2],
                    };
                    for (std::size_t i = 0; i < m_tets.size(), ++i) {
                        const auto tet_aabb = m_tets[i].AABB();
                        if (basicshape::AABB::overlap(aabb, tet_aabb)) {
                            m_grid[idx].push_back(i);
                        }
                    }
                }
            }
        }
    }

private:
    std::array<T, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };
    std::array<T, 6> m_spacing = { 1, 1, 1 };
    std::array<std::vector<std::size_t>, N * N * N> m_grid;
    std::vector<Tetrahedron<T>> m_tets;
};
}