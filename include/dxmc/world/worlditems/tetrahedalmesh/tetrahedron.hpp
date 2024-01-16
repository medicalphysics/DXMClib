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

Copyright 2024 Erlend Andersen
*/

#pragma once

#include "dxmc/floating.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/dosescore.hpp"
#include "dxmc/world/energyscore.hpp"
#include "dxmc/world/worldintersectionresult.hpp"

#include <algorithm>
#include <array>
#include <execution>

namespace dxmc {

template <Floating T>
class Tetrahedron {
public:
    Tetrahedron(const std::array<T, 3>& first, const std::array<T, 3>& second, const std::array<T, 3>& third, const std::array<T, 3>& fourth, std::uint16_t collectionIdx = 0, std::uint16_t materialIdx = 0)
        : m_collectionIdx(collectionIdx)
        , m_materialIdx(materialIdx)
    {
        m_vertices[0] = first;
        m_vertices[1] = second;
        m_vertices[2] = third;
        m_vertices[3] = fourth;
    }

    Tetrahedron(const std::array<std::array<T, 3>, 4>& vertices, std::uint16_t collectionIdx = 0, std::uint16_t materialIdx = 0)
        : m_vertices(vertices)
        , m_collectionIdx(collectionIdx)
        , m_materialIdx(materialIdx)
    {
    }

    Tetrahedron()
    {
    }

    std::uint16_t collection() const { return m_collectionIdx; }
    void setCollection(std::uint16_t coll) { m_collectionIdx = coll; }
    std::uint16_t materialIndex() const { return m_materialIdx; }
    void setMaterialIndex(std::uint16_t idx) { m_materialIdx = idx; }

    auto operator<=>(const Tetrahedron<T>& other) const = default;

    std::array<T, 3> center() const
    {
        const auto c_sum = vectormath::add(m_vertices[0], m_vertices[1], m_vertices[2], m_vertices[3]);
        return vectormath::scale(c_sum, T { 0.25 });
    }

    void translate(const std::array<T, 3>& dist)
    {
        std::for_each(std::execution::unseq, m_vertices.begin(), m_vertices.end(), [&](auto& vert) {
            for (std::size_t i = 0; i < 3; ++i) {
                vert[i] += dist[i];
            }
        });
    }

    void rotate(const std::array<T, 3>& axis, T angle)
    {
        std::transform(std::execution::unseq, m_vertices.cbegin(), m_vertices.cend(), m_vertices.begin(), [&axis, angle](const auto& v) {
            return vectormath::rotate(v, axis, angle);
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

    T volume() const
    {
        const auto a = vectormath::subtract(m_vertices[1], m_vertices[0]);
        const auto b = vectormath::subtract(m_vertices[2], m_vertices[0]);
        const auto c = vectormath::subtract(m_vertices[3], m_vertices[0]);
        static constexpr T scale = 1 / T { 6 };
        return scale * std::abs(vectormath::tripleProduct(a, b, c));
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

        std::array<T, 2> t;
        int n_hits = 0;
        const auto h3 = intersectTriangle(m_vertices[0], m_vertices[1], m_vertices[2], particle);
        if (h3) {
            t[0] = *h3;
            n_hits++;
        }
        const auto h2 = intersectTriangle(m_vertices[1], m_vertices[0], m_vertices[3], particle);
        if (h2) {
            t[n_hits] = *h2;
            n_hits++;
        }
        if (n_hits < 2) {
            const auto h1 = intersectTriangle(m_vertices[2], m_vertices[3], m_vertices[0], particle);
            if (h1) {
                t[n_hits] = *h1;
                n_hits++;
            }
            if (n_hits == 1) {
                const auto h0 = intersectTriangle(m_vertices[3], m_vertices[2], m_vertices[1], particle);
                if (h0) {
                    t[n_hits] = *h0;
                    n_hits++;
                }
            }
        }
        WorldIntersectionResult<T> res;
        if (n_hits == 2) {
            const auto [min, max] = std::minmax(t[0], t[1]);
            res.rayOriginIsInsideItem = min <= 0 && max > 0;
            res.intersection = res.rayOriginIsInsideItem ? max : min;
            res.intersectionValid = res.intersection >= 0;
        }
        return res;
    }

    std::array<T, 3> normal(const std::array<T, 3>& point) const
    {
        auto distance = [](const std::array<T, 3>& pointf, const std::array<T, 3>& pointinplane, const std::array<T, 3>& normal_arr) {
            return std::abs(vectormath::dot(vectormath::subtract(pointf, pointinplane), normal_arr));
        };
        const std::array<std::array<T, 3>, 4> normals = {
            normalVector(m_vertices[0], m_vertices[1], m_vertices[2]),
            normalVector(m_vertices[1], m_vertices[0], m_vertices[3]),
            normalVector(m_vertices[2], m_vertices[3], m_vertices[0]),
            normalVector(m_vertices[3], m_vertices[2], m_vertices[1])
        };
        const std::array<T, 4> dist = {
            distance(point, m_vertices[0], normals[0]),
            distance(point, m_vertices[1], normals[1]),
            distance(point, m_vertices[2], normals[2]),
            distance(point, m_vertices[3], normals[3])
        };
        const auto idx = std::distance(dist.cbegin(), std::min_element(dist.cbegin(), dist.cend()));

        return normals[idx];
    }

    bool pointInside(const std::array<T, 3>& point) const
    {
        // F3 (V0V1V2), F2 (V1V0V3), F1 (V2V3V0), F0 (V3V2V1)
        /*const std::array<std::array<T, 3>, 4> normals = {
            normalVector<false>(m_vertices[0], m_vertices[1], m_vertices[2]),
            normalVector<false>(m_vertices[1], m_vertices[0], m_vertices[3]),
            normalVector<false>(m_vertices[2], m_vertices[3], m_vertices[0]),
            normalVector<false>(m_vertices[3], m_vertices[2], m_vertices[1])
        };*/

        const std::array<std::array<T, 3>, 4> normals = {
            normalVector<true>(m_vertices[0], m_vertices[1], m_vertices[2]),
            normalVector<true>(m_vertices[1], m_vertices[0], m_vertices[3]),
            normalVector<true>(m_vertices[2], m_vertices[3], m_vertices[0]),
            normalVector<true>(m_vertices[3], m_vertices[2], m_vertices[1])
        };

        bool inside = true;
        for (int i = 0; i < 4; ++i) {
            inside = inside && vectormath::dot(vectormath::subtract(m_vertices[i], point), normals[i]) >= 0;
        }
        return inside;
    }

    bool validVerticeOrientation() const
    {
        const auto c = center();
        return pointInside(c);
    }

    const DoseScore<T>& doseScored() const
    {
        return m_dose;
    }

    const EnergyScore<T>& energyImparted() const
    {
        return m_energy_imparted;
    }

    void clearEnergyScored()
    {
        m_energy_imparted.clear();
    }
    void clearDoseScored()
    {
        m_dose.clear();
    }

    void scoreEnergy(T energy)
    {
        m_energy_imparted.scoreEnergy(energy);
    }

    void addEnergyScoredToDoseScore(T density, T calibration_factor = 1)
    {
        m_dose.addScoredEnergy(m_energy_imparted, volume(), density, calibration_factor);
    }

protected:
    template <bool NORMALIZE = true>
    static std::array<T, 3> normalVector(const std::array<T, 3>& p0, const std::array<T, 3>& p1, const std::array<T, 3>& p2)
    {
        const auto s1 = vectormath::subtract(p1, p0);
        const auto s2 = vectormath::subtract(p2, p0);
        if constexpr (NORMALIZE)
            return vectormath::normalized(vectormath::cross(s2, s1));
        else
            return vectormath::cross(s2, s1);
    }

    static std::optional<T> intersectTriangle(const std::array<T, 3>& v0, const std::array<T, 3>& v1, const std::array<T, 3>& v2, const Particle<T>& p)
    {
        const auto E1 = vectormath::subtract(v1, v0);
        const auto TT = vectormath::subtract(p.pos, v0);
        const auto Q = vectormath::cross(TT, E1);

        const auto E2 = vectormath::subtract(v2, v0);
        const auto P = vectormath::cross(p.dir, E2);

        const auto det_inv = 1 / vectormath::dot(P, E1);

        const auto v = vectormath::dot(Q, p.dir) * det_inv;
        const auto u = vectormath::dot(P, TT) * det_inv;

        return v >= 0 && u >= 0 && (u + v) <= 1 ? std::make_optional(vectormath::dot(Q, E2) * det_inv)
                                                : std::nullopt;
    }

private:
    std::array<std::array<T, 3>, 4> m_vertices;
    EnergyScore<T> m_energy_imparted;
    DoseScore<T> m_dose;
    std::uint16_t m_collectionIdx = 0;
    std::uint16_t m_materialIdx = 0;
};
}