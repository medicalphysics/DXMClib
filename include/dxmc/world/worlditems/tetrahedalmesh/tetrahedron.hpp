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

namespace dxmc {

template <Floating T, typename U>
struct TetrahedalMeshIntersectionResult {
    const U* item = nullptr;
    T intersection;
    std::array<T, 3> normal;
    bool rayOriginIsInsideItem = false;

    inline bool valid() const
    {
        return item != nullptr;
    }
};

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

    TetrahedalMeshIntersectionResult<T, Tetrahedron<T>> intersect(const Particle<T>& particle) const
    {
        const auto res = forwardIntersect(particle);
        TetrahedalMeshIntersectionResult<T, Tetrahedron<T>> w;
        if (res.valid && res.t_exit > 0) {
            w.rayOriginIsInsideItem = res.t_enter < 0;
            w.intersection = w.rayOriginIsInsideItem ? res.t_exit : res.t_enter;
            w.normal = w.rayOriginIsInsideItem ? res.normal_exit : res.normal_enter;
            w.item = this;
        }
        return w;
    }

protected:
    static std::array<T, 3> normalVector(const std::array<T, 3>& p0, const std::array<T, 3>& p1, const std::array<T, 3>& p2)
    {
        const auto s1 = vectormath::subtract(p1, p0);
        const auto s2 = vectormath::subtract(p2, p0);
        auto normal = vectormath::cross(s1, s2);
        vectormath::normalize(normal);
        return normal;
    }

    static T planeIntersect(const Particle<T>& p, const std::array<T, 3>& point, const std::array<T, 3>& normal)
    {
        // we assume ray start is in origo
        // point is point on plane

        // if not ray in origo
        // return vectormath::dot(vectormath::subtract(point, p.pos), normal) / vectormath::dot(p.dir, normal);
        return vectormath::dot(point, normal) / vectormath::dot(p.dir, normal);
    }

    static T planeIntersect(const Particle<T>& p, const std::array<T, 3>& v0, const std::array<T, 3>& v1, const std::array<T, 3>& v2)
    {
        // we assume ray start is in origo
        // v0 v1 v2 are point on plane
        const auto normal = normalVector(v0, v1, v2);
        // if not ray in origo
        // const auto t = vectormath::dot(vectormath::subtract(v0, p.pos), normal) / vectormath::dot(p.dir, normal);
        const auto t = vectormath::dot(v0, normal) / vectormath::dot(p.dir, normal);
        return t;
    }

    struct IntersectResult {
        T t_enter, t_exit;
        std::array<T, 3> normal_enter, normal_exit;
        bool valid = false;
    };

    IntersectResult forwardIntersect(const Particle<T>& p) const
    {
        // translate such as p.pos lies at origo

        const auto A = vectormath::subtract(m_vertices[0], p.pos);
        const auto B = vectormath::subtract(m_vertices[1], p.pos);
        const auto C = vectormath::subtract(m_vertices[2], p.pos);
        const auto D = vectormath::subtract(m_vertices[3], p.pos);

        // Let A B C D define tetrahedron such that faces are clockwise when viewed from front (outside)
        // i.e F3 = ABC, F2 = BAD, F1 = CDA, F0 = DCB

        // Write triple scalar product = (p * (q x r)) as [p q r] and Q = p.dir

        constexpr T e = sizeof(T) > 4 ? 1e-7 : 1e-5f;
        constexpr T ne = -e;
        IntersectResult res;

        bool noEnter = true;
        bool noExit = true;

        // test for F3
        // [Q A B] >= 0 && [Q B C] >= 0 && [Q C A] >= 0
        const auto QAB = vectormath::tripleProduct(p.dir, A, B);
        const auto QBC = vectormath::tripleProduct(p.dir, B, C);
        const auto QCA = vectormath::tripleProduct(p.dir, C, A);
        if (noEnter && QAB >= ne && QBC >= ne && QCA >= ne) {
            // enter
            noEnter = false;
            res.normal_enter = normalVector(A, B, C);
            res.t_enter = planeIntersect(p, A, res.normal_enter);
        } else if (noExit && QAB <= e && QBC <= e && QCA <= e) {
            // exit
            noExit = false;
            // t[1] = planeIntersect(p, A, B, C);
            res.normal_exit = normalVector(A, B, C);
            res.t_exit = planeIntersect(p, A, res.normal_exit);
        }

        // test for F2
        //  [Q B A] >= 0 && [Q A D] >= 0 && [Q D B] >= 0
        const auto QBA = -QAB; // tp(p.dir, B, A);
        const auto QAD = vectormath::tripleProduct(p.dir, A, D);
        const auto QDB = vectormath::tripleProduct(p.dir, D, B);
        if (noEnter && QBA >= ne && QAD >= ne && QDB >= ne) {
            // enter
            noEnter = false;
            // t[0] = planeIntersect(p, B, A, D);
            res.normal_enter = normalVector(B, A, D);
            res.t_enter = planeIntersect(p, B, res.normal_enter);
        } else if (noExit && QBA <= e && QAD <= e && QDB <= e) {
            // exit
            noExit = false;
            // t[1] = planeIntersect(p, B, A, D);
            res.normal_exit = normalVector(B, A, D);
            res.t_exit = planeIntersect(p, B, res.normal_exit);
        }

        // test for F1
        // [Q C D] >= 0 && [Q D A] >= 0 && [Q A C]
        const auto QCD = vectormath::tripleProduct(p.dir, C, D);
        const auto QDA = -QAD; // tp(p.dir, D, A);
        const auto QAC = -QCA; // tp(p.dir, A, C);
        if (noEnter && QCD >= ne && QDA >= ne && QAC >= ne) {
            // enter
            noEnter = false;
            // t[0] = planeIntersect(p, C, D, A);
            res.normal_enter = normalVector(C, D, A);
            res.t_enter = planeIntersect(p, C, res.normal_enter);
        } else if (noExit && QCD <= e && QDA <= e && QAC <= e) {
            // exit
            noExit = false;
            // t[1] = planeIntersect(p, C, D, A);
            res.normal_exit = normalVector(C, D, A);
            res.t_exit = planeIntersect(p, C, res.normal_exit);
        }

        const auto oneIntersection = (noEnter || noExit) || (!noEnter && !noExit);
        if (oneIntersection) {
            // test for F0
            // [Q D C] >= 0 && [Q C B] >= 0 && [Q B D] >= 0
            const auto QDC = -QCD; // tp(p.dir, D, C);
            const auto QCB = -QBC; // tp(p.dir, C, B);
            const auto QBD = -QDB; // tp(p.dir, B, D);
            if (noEnter && QDC >= ne && QCB >= ne && QBD >= ne) {
                // enter
                noEnter = false;
                // t[0] = planeIntersect(p, D, C, B);
                res.normal_enter = normalVector(D, C, B);
                res.t_enter = planeIntersect(p, D, res.normal_enter);
            } else if (noExit && QDC <= e && QCB <= e && QBD <= e) {
                // exit
                noExit = false;
                // t[1] = planeIntersect(p, D, C, B);
                res.normal_exit = normalVector(D, C, B);
                res.t_exit = planeIntersect(p, D, res.normal_exit);
            }
        }
        if (!noEnter && !noExit) { // two intersections
            res.valid = true;
        }
        return res;
    }

private:
    std::array<std::array<T, 3>, 4> m_vertices;
    std::uint16_t m_collectionIdx = 0;
    std::uint16_t m_materialIdx = 0;
};
}