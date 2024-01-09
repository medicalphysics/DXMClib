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

template <Floating T, std::size_t N = 16>
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
        const std::array<T, 3> scale = { 1 / m_spacing[0], 1 / m_spacing[1], 1 / m_spacing[2] };
        constexpr auto max = std::numeric_limits<std::size_t>::max();

        for (std::size_t i = 0; i < m_tets.size(); ++i) {
            const auto& verts = m_tets[i].vertices();
            std::array<std::size_t, 4> flat_indices = { max, max, max, max };
            for (std::size_t k = 0; k < 4; ++k) {
                const auto d = vectormath::scale(vectormath::subtract(start, verts[k]), scale);
                std::array<std::size_t, 3> indices;
                for (std::size_t j = 0; j < 3; ++j) {
                    if (0 < d[j] && d[j] < N) {
                        indices[j] = static_cast<std::size_t>(d[j]);
                    }
                }
                flat_indices[k] = indices[0] + N * indices[1] + N * N * indices[2];
            }
            std::ranges::sort(flat_indices);
            auto end = std::unique(flat_indices.begin(), flat_indices.end());
            for (auto p = flat_indices.begin(); p != end; ++p) {
                m_grid[*p].push_back(i)
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