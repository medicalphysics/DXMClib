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

#include "dxmc/floating.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/basicshapes/tetrahedron.hpp"
#include "dxmc/world/worldintersectionresult.hpp"

#include <algorithm>
#include <array>
#include <execution>
#include <optional>

namespace dxmc {

template <Floating T>
class Tetrahedron {
public:
    Tetrahedron(const std::array<T, 3>& first, const std::array<T, 3>& second, const std::array<T, 3>& third, const std::array<T, 3>& fourth, std::uint8_t collectionIdx = 0)
        : m_collectionIdx(collectionIdx)
    {
        m_vertices[0] = first;
        m_vertices[1] = second;
        m_vertices[2] = third;
        m_vertices[3] = fourth;
    }

    Tetrahedron(const std::array<std::array<T, 3>, 4>& vertices, std::uint8_t collectionIdx = 0)
        : m_vertices(vertices)
        , m_collectionIdx(collectionIdx)
    {
    }

    Tetrahedron()
    {
    }

    std::uint8_t collection() const { return m_collectionIdx; }

    auto operator<=>(const Tetrahedron<T>& other) const = default;

    void translate(const std::array<T, 3>& dist)
    {
        std::for_each(std::execution::unseq, m_vertices.begin(), m_vertices.end(), [&](auto& vert) {
            for (std::size_t i = 0; i < 3; ++i) {
                vert[i] += dist[i];
            }
        });
    }

    void scale(T scale)
    {
        std::for_each(std::execution::unseq, m_vertices.begin(), m_vertices.end(), [&](auto& vert) {
            for (std::size_t i = 0; i < 3; ++i) {
                vert[i] *= scale;
            }
        });
    }

    const std::array<std::array<T, 3>, 4>& vertices() const
    {
        return m_vertices;
    }

    std::array<T, 3> center() const
    {
        std::array<T, 3> cent { 0, 0, 0 };
        for (const auto& vert : m_vertices) {
            for (std::size_t i = 0; i < 3; i++)
                cent[i] += vert[i];
        }
        constexpr T factor = 1 / T { 4 };
        for (std::size_t i = 0; i < 3; i++)
            cent[i] *= factor;
        return cent;
    }
    std::array<T, 6> AABB() const
    {
        std::array<T, 6> aabb {
            std::numeric_limits<T>::max(),
            std::numeric_limits<T>::max(),
            std::numeric_limits<T>::max(),
            std::numeric_limits<T>::lowest(),
            std::numeric_limits<T>::lowest(),
            std::numeric_limits<T>::lowest(),
        };
        for (std::size_t j = 0; j < 4; j++) {
            for (std::size_t i = 0; i < 3; i++) {
                aabb[i] = std::min(aabb[i], m_vertices[j][i]);
            }
            for (std::size_t i = 0; i < 3; i++) {
                const auto idx = i + 3;
                aabb[idx] = std::max(aabb[idx], m_vertices[j][i]);
            }
        }
        return aabb;
    }

    auto begin() { return m_vertices.begin(); }
    auto begin() const { return m_vertices.begin(); }
    auto cbegin() const { return m_vertices.cbegin(); }
    auto end() { return m_vertices.end(); }
    auto end() const { return m_vertices.end(); }
    auto cend() const { return m_vertices.cend(); }

    WorldIntersectionResult<T> intersect(const Particle<T>& particle) const
    {
        return basicshape::tetrahedron::intersect(particle, m_vertices[0], m_vertices[1], m_vertices[2], m_vertices[3]);
    }
    VisualizationIntersectionResult<T, WorldItemBase<T>> intersectVisualization(const Particle<T>& particle) const
    {
        return basicshape::tetrahedron::template intersectVisualization<T, WorldItemBase<T>>(particle, m_vertices[0], m_vertices[1], m_vertices[2], m_vertices[3]);
    }

private:
    std::array<std::array<T, 3>, 4> m_vertices;
    std::uint8_t m_collectionIdx = 0;
};
}