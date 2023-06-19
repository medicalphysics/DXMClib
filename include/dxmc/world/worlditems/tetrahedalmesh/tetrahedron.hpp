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
    T t_enter, t_exit;
    std::array<T, 3> normal_enter, normal_exit;

    inline bool valid() const
    {
        return item != nullptr;
    }
    inline T intersection() const
    {
        return t_enter > T { 0 } ? t_enter : t_exit;
    }
    inline T rayOriginIsInsideItem() const
    {
        return t_enter <= T { 0 };
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

    inline TetrahedalMeshIntersectionResult<T, Tetrahedron<T>> intersect(const Particle<T>& particle) const
    {
        return forwardIntersect<false>(particle);
    }

    bool validVerticeOrientation() const
    {
        auto center = std::reduce(m_vertices.cbegin(), m_vertices.cend(), std::array<T, 3> { 0, 0, 0 }, [](const auto& lh, const auto& rh) { return vectormath::add(lh, rh); });
        for (auto& c : center)
            c /= T { 4 };

        std::array<T, 3> pv0, pv1;
        std::array<T, 4> proj;

        // F3 (V0V1V2), F2 (V1V0V3), F1 (V2V3V0), F0 (V3V2V1)
        {
            const auto d = vectormath::subtract(center, m_vertices[3]);
            const auto pv0 = vectormath::subtract(m_vertices[1], m_vertices[0]);
            const auto pv1 = vectormath::subtract(m_vertices[2], m_vertices[0]);
            const auto n = vectormath::cross(pv0, pv1);
            proj[3] = vectormath::dot(d, n);
        }
        {
            const auto d = vectormath::subtract(center, m_vertices[2]);
            const auto pv0 = vectormath::subtract(m_vertices[0], m_vertices[1]);
            const auto pv1 = vectormath::subtract(m_vertices[3], m_vertices[1]);
            const auto n = vectormath::cross(pv0, pv1);
            proj[2] = vectormath::dot(d, n);
        }
        {
            const auto d = vectormath::subtract(center, m_vertices[1]);
            const auto pv0 = vectormath::subtract(m_vertices[3], m_vertices[2]);
            const auto pv1 = vectormath::subtract(m_vertices[0], m_vertices[2]);
            const auto n = vectormath::cross(pv0, pv1);
            proj[1] = vectormath::dot(d, n);
        }
        {
            const auto d = vectormath::subtract(center, m_vertices[0]);
            const auto pv0 = vectormath::subtract(m_vertices[2], m_vertices[3]);
            const auto pv1 = vectormath::subtract(m_vertices[1], m_vertices[3]);
            const auto n = vectormath::cross(pv0, pv1);
            proj[0] = vectormath::dot(d, n);
        }

        bool valid = true;
        for (auto p : proj)
            valid = valid && p <= T { 0 };

        return valid;
    }

protected:
    template <bool NORMALIZE = true>
    static std::array<T, 3> normalVector(const std::array<T, 3>& p0, const std::array<T, 3>& p1, const std::array<T, 3>& p2)
    {
        const auto s1 = vectormath::subtract(p1, p0);
        const auto s2 = vectormath::subtract(p2, p0);
        auto normal = vectormath::cross(s1, s2);
        if constexpr (NORMALIZE)
            vectormath::normalize(normal);
        return normal;
    }

    template <bool PROPERNORMAL = true>
    TetrahedalMeshIntersectionResult<T, Tetrahedron<T>> forwardIntersect(const Particle<T>& p) const
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
        TetrahedalMeshIntersectionResult<T, Tetrahedron<T>> res;

        bool noEnter = true;
        bool noExit = true;

        // test for F3
        // [Q A B] >= 0 && [Q B C] >= 0 && [Q C A] >= 0
        const auto QAB = vectormath::tripleProduct(p.dir, A, B);
        const auto QBC = vectormath::tripleProduct(p.dir, B, C);
        const auto QCA = vectormath::tripleProduct(p.dir, C, A);
        if (noEnter && QAB >= ne && QBC >= ne && QCA >= ne) {
            // enter
            res.normal_enter = normalVector<PROPERNORMAL>(A, B, C);
            const auto den = vectormath::dot(p.dir, res.normal_enter);
            if (std::abs(den) > e) {
                noEnter = false;
                const auto t = vectormath::dot(A, res.normal_enter) / den;
                res.t_enter = t;
            }
        } else if (noExit && QAB <= e && QBC <= e && QCA <= e) {
            // exit
            res.normal_exit = normalVector<PROPERNORMAL>(A, B, C);
            const auto den = vectormath::dot(p.dir, res.normal_exit);
            if (std::abs(den) > e) {
                noExit = false;
                const auto t = vectormath::dot(A, res.normal_exit) / den;
                res.t_exit = t;
            }
        }

        // test for F2
        //  [Q B A] >= 0 && [Q A D] >= 0 && [Q D B] >= 0
        const auto QBA = -QAB; // tp(p.dir, B, A);
        const auto QAD = vectormath::tripleProduct(p.dir, A, D);
        const auto QDB = vectormath::tripleProduct(p.dir, D, B);
        if (noEnter && QBA >= ne && QAD >= ne && QDB >= ne) {
            // enter
            res.normal_enter = normalVector<PROPERNORMAL>(B, A, D);
            const auto den = vectormath::dot(p.dir, res.normal_enter);
            if (std::abs(den) > e) {
                noEnter = false;
                const auto t = vectormath::dot(B, res.normal_enter) / den;
                res.t_enter = t;
            }
        } else if (noExit && QBA <= e && QAD <= e && QDB <= e) {
            // exit
            res.normal_exit = normalVector<PROPERNORMAL>(B, A, D);
            const auto den = vectormath::dot(p.dir, res.normal_exit);
            if (std::abs(den) > e) {
                noExit = false;
                const auto t = vectormath::dot(B, res.normal_exit) / den;
                res.t_exit = t;
            }
        }

        // test for F1
        // [Q C D] >= 0 && [Q D A] >= 0 && [Q A C]
        const auto QCD = vectormath::tripleProduct(p.dir, C, D);
        const auto QDA = -QAD; // tp(p.dir, D, A);
        const auto QAC = -QCA; // tp(p.dir, A, C);
        if (noEnter && QCD >= ne && QDA >= ne && QAC >= ne) {
            // enter
            res.normal_enter = normalVector<PROPERNORMAL>(C, D, A);
            const auto den = vectormath::dot(p.dir, res.normal_enter);
            if (std::abs(den) > e) {
                noEnter = false;
                const auto t = vectormath::dot(C, res.normal_enter) / den;
                res.t_enter = t;
            }
        } else if (noExit && QCD <= e && QDA <= e && QAC <= e) {
            // exit
            res.normal_exit = normalVector<PROPERNORMAL>(C, D, A);
            const auto den = vectormath::dot(p.dir, res.normal_exit);
            if (std::abs(den) > e) {
                noExit = false;
                const auto t = vectormath::dot(C, res.normal_exit) / den;
                res.t_exit = t;
            }
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
                res.normal_enter = normalVector<PROPERNORMAL>(D, C, B);
                const auto den = vectormath::dot(p.dir, res.normal_enter);
                if (std::abs(den) > e) {
                    noEnter = false;
                    const auto t = vectormath::dot(D, res.normal_enter) / den;
                    res.t_enter = t;
                }
            } else if (noExit && QDC <= e && QCB <= e && QBD <= e) {
                // exit
                res.normal_exit = normalVector<PROPERNORMAL>(D, C, B);
                const auto den = vectormath::dot(p.dir, res.normal_exit);
                if (std::abs(den) > e) {
                    noExit = false;
                    const auto t = vectormath::dot(D, res.normal_exit) / den;
                    res.t_exit = t;
                }
            }
        }
        if (!noEnter && !noExit) { // two intersections
            if (res.t_exit > res.t_enter && res.t_exit > T { 0 }) // we ignore glancing hits
                res.item = this;
        }
        return res;
    }

private:
    std::array<std::array<T, 3>, 4> m_vertices;
    std::uint16_t m_collectionIdx = 0;
    std::uint16_t m_materialIdx = 0;
};
}