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

#include <algorithm>
#include <array>
#include <execution>
#include <optional>

namespace dxmc {

template <Floating T>
class Triangle {
public:
    Triangle(const std::array<T, 3>& first, const std::array<T, 3>& second, const std::array<T, 3>& third)
    {
        m_vertices[0] = first;
        m_vertices[1] = second;
        m_vertices[2] = third;
    }
    Triangle(const std::array<std::array<T, 3>, 3>& vertices)
        : m_vertices(vertices)
    {
    }
    Triangle(const T* first_element)
    {
        for (std::size_t i = 0; i < 3; ++i)
            for (std::size_t j = 0; j < 3; ++j) {
                const auto flatIdx = i * 3 + j;
                m_vertices[i][j] = *(first_element + flatIdx);
            }
    }
    auto operator<=>(const Triangle<T>& other) const = default;
    
    void translate(const std::array<T, 3>& dist)
    {
        std::for_each(std::execution::par_unseq, m_vertices.begin(), m_vertices.end(), [&](auto& vert) {
            for (std::size_t i = 0; i < 3; ++i) {
                vert[i] += dist[i];
            }
        });
    }
    void scale(T scale)
    {
        std::for_each(std::execution::par_unseq, m_vertices.begin(), m_vertices.end(), [&](auto& vert) {
            for (std::size_t i = 0; i < 3; ++i) {
                vert[i] *= scale;
            }
        });
    }
    
    const std::array<std::array<T, 3>, 3>& vertices() const
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
        constexpr T factor { 1 / 3.0 };
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
        for (std::size_t j = 0; j < 3; j++) {
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

    template <int FORWARD = 1>
    std::optional<T> intersect(const Particle<T>& particle) const
    {
        const auto& v1 = m_vertices[0];
        const auto& v2 = m_vertices[1];
        const auto& v3 = m_vertices[2];

        const std::array<T, 3> v1v2 { v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2] };
        const std::array<T, 3> v1v3 { v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2] };
        const auto& pvec = vectormath::cross(particle.dir, v1v3);
        const auto det = vectormath::dot(v1v2, pvec);

        if (std::abs(det) <= std::numeric_limits<T>::epsilon())
            return std::nullopt;

        const T invDet = T { 1 } / det;

        const std::array<T, 3> tvec { particle.pos[0] - v1[0], particle.pos[1] - v1[1], particle.pos[2] - v1[2] };
        const T u = vectormath::dot(tvec, pvec) * invDet;
        if (u < T { 0 } || u > T { 1 })
            return std::nullopt;

        const auto qvec = vectormath::cross(tvec, v1v2);
        const T v = vectormath::dot(particle.dir, qvec) * invDet;
        if (v < T { 0 } || u + v > T { 1 })
            return std::nullopt;

        const auto t = vectormath::dot(v1v3, qvec) * invDet;
        if constexpr (FORWARD == 1)
            return t > T { 0 } ? std::make_optional(t) : std::nullopt;
        else if constexpr (FORWARD == -1)
            return t < T { 0 } ? std::make_optional(t) : std::nullopt;
        else
            return std::make_optional(t);
    }

private:
    std::array<std::array<T, 3>, 3> m_vertices;
};
}