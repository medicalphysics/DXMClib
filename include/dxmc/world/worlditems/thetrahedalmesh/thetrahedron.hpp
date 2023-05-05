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
#include "dxmc/vectormath.hpp"

#include <algorithm>
#include <array>
#include <optional>

namespace dxmc {

template <Floating T>
class Thetrahedron {
public:
    Thetrahedron(const std::array<T, 3>& first, const std::array<T, 3>& second, const std::array<T, 3>& third, const std::array<T, 3>& fourth)
    {
        m_vertices[0] = first;
        m_vertices[1] = second;
        m_vertices[2] = third;
        m_vertices[3] = fourth;
    }

    Thetrahedron(const std::array<std::array<T, 3>, 4>& vertices)
        : m_vertices(vertices)
    {
    }

    Thetrahedron(const T* first_element)
    {
        for (std::size_t i = 0; i < 4; ++i)
            for (std::size_t j = 0; j < 3; ++j) {
                const auto flatIdx = i * 3 + j;
                m_vertices[i][j] = *(first_element + flatIdx);
            }
    }

    auto operator<=>(const Thetrahedron<T>& other) const = default;

    void translate(const std::array<T, 3>& dist)
    {
        std::for_each(m_vertices.begin(), m_vertices.end(), [&](auto& vert) {
            for (std::size_t i = 0; i < 3; ++i) {
                vert[i] += dist[i];
            }
        });
    }

    void scale(T scale)
    {
        std::for_each(m_vertices.begin(), m_vertices.end(), [&](auto& vert) {
            for (std::size_t i = 0; i < 3; ++i) {
                vert[i] *= scale;
            }
        });
    }

    const std::array<std::array<T, 4>, 3>& vertices() const
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
        constexpr T factor { 1 / 4.0 };
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

    std::optional<T> intersect(const Particle<T>& particle) const
    {
        // Plucker coordinates
        std::array<std::array<T, 3>, 2> p = { particle.dir, vectormath::cross(particle.dir, particle.pos) };

        auto s0_dir = vectormath::subtract(m_vertices[2], m_vertices[1]);
        std::array<std::array<T, 3>, 2> s0 = { s0_dir, vectormath::cross(s0_dir, m_vertices[1]) };
        auto s1_dir = vectormath::subtract(m_vertices[0], m_vertices[2]);
        std::array<std::array<T, 3>, 2> s1 = { s1_dir, vectormath::cross(s1_dir, m_vertices[2]) };
        auto s2_dir = vectormath::subtract(m_vertices[1], m_vertices[0]);
        std::array<std::array<T, 3>, 2> s2 = { s2_dir, vectormath::cross(s2_dir, m_vertices[0]) };

        auto r0 = inner_product(p, s0);
        auto r1 = inner_product(p, s1);
        auto r2 = inner_product(p, s2);

        constexpr T e = sizeof(T) == 4 ? 1e-5 : 1E-7;

        if ((r0 >= -e && r1 >= -e && r2 >= -e) || (r0 <= e && r1 <= e && r2 <= e)) {
            const auto idx = vectormath::argmax3<std::uint_fast32_t>(particle.dir);
            const auto pend = (r0 * m_vertices[0][idx] + r1 * m_vertices[1][idx] + r2 * m_vertices[2][idx]) / (r0 + r1 + r2);
            auto t = (pend - particle.pos[idx]) / particle.dir[idx];

            return std::make_optional(t);
        }
        return std::nullopt;
    }

protected:
    static T inner_product(const std::array<std::array<T, 3>, 2>& lh, const std::array<std::array<T, 3>, 2>& rh)
    {
        return vectormath::dot(lh[0], rh[1]) + vectormath::dot(lh[1], rh[0]);
    }

private:
    std::array<std::array<T, 3>, 3> m_vertices;
};
}