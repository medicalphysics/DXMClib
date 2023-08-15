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

#include <algorithm>
#include <array>
#include <execution>
#include <memory>
#include <optional>
#include <vector>

namespace dxmc {

template <typename U, typename T>
concept MeshKDTreeType = requires(U u, Particle<T> p, std::array<T, 3> vec, T scale) {
                             Floating<T>;
                             u <=> u;
                             u.translate(vec);
                             u.scale(scale);
                             {
                                 u.intersect(p)
                                 } -> std::same_as<std::optional<T>>;
                             {
                                 u.center()
                                 } -> std::same_as<std::array<T, 3>>;
                             {
                                 u.AABB()
                                 } -> std::same_as<std::array<T, 6>>;
                             {
                                 u.planeVector()
                                 } -> std::same_as<std::array<T, 3>>;
                         };

template <Floating T, MeshKDTreeType<T> U>
class MeshKDTree {
public:
    MeshKDTree() {};
    MeshKDTree(MeshKDTree<T, U>& other) = delete;
    MeshKDTree(const MeshKDTree<T, U>& other) = delete;

    MeshKDTree(MeshKDTree<T, U>&& other) noexcept
    {
        m_D = other.m_D;
        m_plane = other.m_plane;
        m_triangles = std::move(other.m_triangles);
        if (other.m_left) {
            m_left = std::unique_ptr(std::move(other.m_left));
        }
        if (other.m_right) {
            m_right = std::unique_ptr(std::move(other.m_right));
        }
    }

    MeshKDTree(const std::vector<U>& triangles, const std::size_t max_depth = 8)
    {
        setData(triangles, max_depth);
    }

    void setData(const std::vector<U>& triangles, const std::size_t max_depth = 8)
    {
        if (triangles.size() == 0)
            return;
        // finding aabb
        std::array<T, 6> aabb {
            std::numeric_limits<T>::max(),
            std::numeric_limits<T>::max(),
            std::numeric_limits<T>::max(),
            std::numeric_limits<T>::lowest(),
            std::numeric_limits<T>::lowest(),
            std::numeric_limits<T>::lowest(),
        };

        for (const auto& tri : triangles) {
            std::array<T, 6> aabb_tri;

            aabb_tri = tri.AABB();

            for (std::size_t i = 0; i < 3; ++i) {
                aabb[i] = std::min(aabb[i], aabb_tri[i]);
            }
            for (std::size_t i = 3; i < 6; ++i) {
                aabb[i] = std::max(aabb[i], aabb_tri[i]);
            }
        }
        const std::array<T, 3> extent { aabb[3] - aabb[0], aabb[4] - aabb[1], aabb[5] - aabb[2] };

        m_D = vectormath::argmax3<std::uint_fast32_t, T>(extent);

        const auto split = planeSplit(triangles);

        const auto fom = figureOfMerit(triangles, split);

        if (fom == triangles.size() || max_depth <= 1 || triangles.size() <= 1) {
            m_triangles = triangles;
            m_left = nullptr;
            m_right = nullptr;
        } else {
            m_plane = split;
            std::vector<U> left;
            std::vector<U> right;
            for (const auto& triangle : triangles) {
                const auto side = planeSide(triangle, m_plane, m_D);
                if (side <= 0)
                    left.push_back(triangle);
                if (side >= 0)
                    right.push_back(triangle);
            }
            m_left = std::make_unique<MeshKDTree<T, U>>(left, max_depth - 1);
            m_right = std::make_unique<MeshKDTree<T, U>>(right, max_depth - 1);
        }
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
        AABB_iterator(aabb);
        return aabb;
    }

    std::size_t depth() const
    {
        std::size_t teller = 0;
        depth_iterator(teller);
        return teller;
    }

    std::vector<U> items() const
    {
        std::vector<U> all;
        item_iterator(all);
        std::sort(all.begin(), all.end());
        auto last = std::unique(all.begin(), all.end());
        all.erase(last, all.end());
        return all;
    }

    void translate(const std::array<T, 3>& dist)
    {
        m_plane += dist[m_D];
        if (m_left) {
            m_left->translate(dist);
            m_right->translate(dist);
        } else {
            std::for_each(std::execution::par_unseq, m_triangles.begin(), m_triangles.end(), [&](auto& tri) {
                tri.translate(dist);
            });
        }
    }

    KDTreeIntersectionResult<T, const U> intersect(const Particle<T>& particle, const std::array<T, 6>& aabb) const
    {
        const auto inter = basicshape::AABB::intersectForwardInterval(particle, aabb);
        return inter ? intersect(particle, *inter) : KDTreeIntersectionResult<T, const U> {};
    }

protected:
    KDTreeIntersectionResult<T, const U> intersect(const Particle<T>& particle, const std::array<T, 2>& tbox) const
    {

        if (!m_left) { // this is a leaf
            // intersect triangles between tbox and return;
            KDTreeIntersectionResult<T, const U> res;
            res.intersection = std::numeric_limits<T>::max();

            for (const U& triangle : m_triangles) {
                const auto t_cand = triangle.intersect(particle);
                if (t_cand) {
                    if (*t_cand < res.intersection) {
                        if (tbox[0] <= *t_cand && *t_cand <= tbox[1]) {
                            res.intersection = *t_cand;
                            res.item = &triangle;
                            res.rayOriginIsInsideItem = vectormath::dot(particle.dir, triangle.planeVector()) > 0;
                        }
                    }
                }
            }
            return res;
        }

        // test for parallell beam
        if (std::abs(particle.dir[m_D]) <= std::numeric_limits<T>::epsilon()) {
            auto hit_left = m_left->intersect(particle, tbox);
            auto hit_right = m_right->intersect(particle, tbox);
            if (hit_left.valid() && hit_right.valid())
                return hit_left.intersection > hit_right.intersection ? hit_right : hit_left;
            if (hit_right.valid())
                return hit_right;
            return hit_left;
        }

        const MeshKDTree<T, U>* const front = particle.dir[m_D] > T { 0 } ? m_left.get() : m_right.get();
        const MeshKDTree<T, U>* const back = particle.dir[m_D] > T { 0 } ? m_right.get() : m_left.get();

        const auto t = (m_plane - particle.pos[m_D]) / particle.dir[m_D];

        if (t <= tbox[0]) {
            // back only
            return back->intersect(particle, tbox);
        } else if (t >= tbox[1]) {
            // front only
            return front->intersect(particle, tbox);
        }

        // both directions (start with front)
        const std::array<T, 2> t_front { tbox[0], t };
        auto hit = front->intersect(particle, t_front);
        if (hit.valid()) {
            if (hit.intersection <= t) {
                return hit;
            }
        }
        const std::array<T, 2> t_back { t, tbox[1] };
        return back->intersect(particle, t_back);
    }

    T planeSplit(const std::vector<U>& triangles) const
    {
        const auto N = triangles.size();
        std::vector<T> vals;
        vals.reserve(N);

        for (const auto& triangle : triangles) {
            const auto v = triangle.center()[m_D];
            vals.push_back(v);
        }
        std::sort(vals.begin(), vals.end());

        if (N % 2 == 1) {
            return vals[N / 2];
        } else {
            return (vals[N / 2] + vals[N / 2 - 1]) * T { 0.5 };
        }
    }

    int figureOfMerit(const std::vector<U>& triangles, const T planesep) const
    {
        int fom = 0;
        int shared = 0;
        for (const auto& triangle : triangles) {
            const auto side = planeSide(triangle, planesep, m_D);
            fom += side;
            if (side == 0) {
                shared++;
            }
        }
        return std::abs(fom) + shared;
    }
    static std::int_fast8_t planeSide(const U& triangle, const T plane, const std::uint_fast32_t D)
    {
        T max = std::numeric_limits<T>::lowest();
        T min = std::numeric_limits<T>::max();

        const auto& aabb = triangle.AABB();
        min = aabb[D];
        max = aabb[D + 3];

        if (lessOrEqual(max, plane))
            // if (max - plane <= epsilon())
            return -1;
        if (greaterOrEqual(min, plane))
            // if (plane - min <= epsilon())
            return 1;
        return 0;
    }
    void depth_iterator(std::size_t& teller) const
    {
        teller++;
        if (m_left)
            m_left->depth_iterator(teller);
    }
    void item_iterator(std::vector<U>& all) const
    {
        if (m_left) {
            m_left->item_iterator(all);
            m_right->item_iterator(all);
        } else {
            std::copy(m_triangles.cbegin(), m_triangles.cend(), std::back_inserter(all));
        }
    }
    void AABB_iterator(std::array<T, 6>& aabb) const
    {
        if (!m_left) {
            for (const auto& tri : m_triangles) {
                const auto aabb_tri = tri.AABB();
                for (std::size_t i = 0; i < 3; ++i) {
                    aabb[i] = std::min(aabb[i], aabb_tri[i]);
                }
                for (std::size_t i = 3; i < 6; ++i) {
                    aabb[i] = std::max(aabb[i], aabb_tri[i]);
                }
            }
        } else {
            m_left->AABB_iterator(aabb);
            m_right->AABB_iterator(aabb);
        }
    }
    constexpr static T epsilon()
    {
        // Huristic epsilon for triangle intersections
        return T { 11 } * std::numeric_limits<T>::epsilon();
    }
    constexpr static bool lessOrEqual(T a, T b)
    {
        return a - b <= epsilon() * a;
    }
    constexpr static bool greaterOrEqual(T a, T b)
    {
        return b - a <= epsilon() * a;
    }

private:
    std::uint_fast32_t m_D = 0;
    T m_plane = 0;
    std::vector<U> m_triangles;
    std::unique_ptr<MeshKDTree<T, U>> m_left = nullptr;
    std::unique_ptr<MeshKDTree<T, U>> m_right = nullptr;
    friend class MeshKDTree<T, U>;
};

}