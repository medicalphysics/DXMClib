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
#include <array>
#include <execution>
#include <memory>
#include <optional>
#include <vector>

namespace dxmc {

template <Floating T>
class TetrahedalMeshKDTree {
public:
    TetrahedalMeshKDTree() {};
    TetrahedalMeshKDTree(TetrahedalMeshKDTree<T>& other) = delete;
    TetrahedalMeshKDTree(const TetrahedalMeshKDTree<T>& other) = delete;

    TetrahedalMeshKDTree(TetrahedalMeshKDTree<T>&& other)
    {
        m_D = other.m_D;
        m_plane = other.m_plane;
        m_tets = std::move(other.m_tets);
        if (other.m_left) {
            m_left = std::unique_ptr(std::move(other.m_left));
        }
        if (other.m_right) {
            m_right = std::unique_ptr(std::move(other.m_right));
        }
    }

    TetrahedalMeshKDTree(std::vector<Tetrahedron<T>>&& tets, int max_depth = 8)
    {
        setData(std::move(tets), max_depth);
    }

    void setData(std::vector<Tetrahedron<T>>&& tets, int max_depth = 8)
    {
        if (tets.size() == 0)
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
        for (const auto& tet : tets) {
            for (const auto& v : tet) {
                for (std::size_t i = 0; i < 3; ++i) {
                    aabb[i] = std::min(aabb[i], v[i]);
                    aabb[i + 3] = std::max(aabb[i + 3], v[i]);
                }
            }
        }
        const std::array<T, 3> extent { aabb[3] - aabb[0], aabb[4] - aabb[1], aabb[5] - aabb[2] };

        m_D = vectormath::argmax3<std::uint_fast32_t, T>(extent);

        const auto split = planeSplit(tets);

        const auto fom = figureOfMerit(tets, split);

        if (fom == tets.size() || max_depth <= 1 || tets.size() <= 1) {
            m_tets = std::move(tets);
            m_left = nullptr;
            m_right = nullptr;
        } else {
            m_plane = split;
            std::vector<Tetrahedron<T>> left;
            left.reserve(tets.size() / 2);
            std::vector<Tetrahedron<T>> right;
            right.reserve(tets.size() / 2);
            for (const auto& tet : tets) {
                const auto side = planeSide(tet, m_plane, m_D);
                if (side <= 0)
                    left.push_back(tet);
                if (side >= 0)
                    right.push_back(tet);
            }
            tets.clear();
            tets.shrink_to_fit(); // clear does not free memory, we shrink to avoid memory explosion
            m_left = std::make_unique<TetrahedalMeshKDTree<T>>(std::move(left), max_depth - 1);
            m_right = std::make_unique<TetrahedalMeshKDTree<T>>(std::move(right), max_depth - 1);
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

    template <std::regular_invocable<Tetrahedron<T>> F>
    void apply(const F& func, auto policy = std::execution::par_unseq)
    {
        if (!m_left) {
            std::for_each(policy, m_tets.begin(), m_tets.end(), func);
        } else {
            m_left->apply(func);
            m_right->apply(func);
        }
    }

    std::size_t depth() const
    {
        std::size_t teller = 0;
        depth_iterator(teller);
        return teller;
    }

    void translate(const std::array<T, 3>& dist)
    {
        m_plane += dist[m_D];
        if (m_left) {
            m_left->translate(dist);
            m_right->translate(dist);
        } else {
            std::for_each(std::execution::par_unseq, m_tets.begin(), m_tets.end(), [&](auto& tri) {
                tri.translate(dist);
            });
        }
    }

    template <std::uint16_t COLLECTION = 65535>
    KDTreeIntersectionResult<T, const Tetrahedron<T>> intersect(const Particle<T>& particle, const std::array<T, 6>& aabb) const
    {
        const auto inter = basicshape::AABB::intersectForwardInterval<T, false>(particle, aabb);
        return inter ? intersect<COLLECTION>(particle, *inter) : KDTreeIntersectionResult<T, const Tetrahedron<T>> {};
    }

protected:
    template <std::uint16_t COLLECTION = 65535>
    KDTreeIntersectionResult<T, const Tetrahedron<T>> intersect(const Particle<T>& particle, const std::array<T, 2>& tbox) const
    {
        if (!m_left) { // this is a leaf
            // intersect tetrahedrons between tbox and return;
            KDTreeIntersectionResult<T, const Tetrahedron<T>> res;

            res.intersection = std::numeric_limits<T>::max();

            for (const auto& tet : m_tets) {
                if constexpr (COLLECTION == 65535) {
                    const auto res_cand = tet.intersect(particle);
                    if (res_cand.valid() && res_cand.intersection < res.intersection) {
                        res.intersection = res_cand.intersection;
                        res.rayOriginIsInsideItem = res_cand.rayOriginIsInsideItem;
                        res.item = &tet;
                        if (res.rayOriginIsInsideItem) // early exit if we are inside tet
                            return res;
                    }
                } else {
                    if (tet.collection == COLLECTION) {
                        const auto res_cand = tet.intersect(particle);
                        if (res_cand.valid() && res_cand.intersection < res.intersection) {
                            res.intersection = res_cand.intersection;
                            res.rayOriginIsInsideItem = res_cand.rayOriginIsInsideItem;
                            res.item = &tet;
                            if (res.rayOriginIsInsideItem) // early exit if we are inside tet
                                return res;
                        }
                    }
                }
            }
            return res;
        }

        // test for parallell beam
        if (std::abs(particle.dir[m_D]) <= std::numeric_limits<T>::epsilon()) {
            auto hit_left = m_left->template intersect<COLLECTION>(particle, tbox);
            auto hit_right = m_right->template intersect<COLLECTION>(particle, tbox);
            if (hit_left.valid() && hit_right.valid())
                return hit_left.intersection > hit_right.intersection ? hit_right : hit_left;
            if (hit_right.valid())
                return hit_right;
            return hit_left;
        }

        const TetrahedalMeshKDTree<T>* const front = particle.dir[m_D] > T { 0 } ? m_left.get() : m_right.get();
        const TetrahedalMeshKDTree<T>* const back = particle.dir[m_D] > T { 0 } ? m_right.get() : m_left.get();

        const auto t = (m_plane - particle.pos[m_D]) / particle.dir[m_D];

        if (t <= tbox[0]) {
            // back only
            return back->intersect<COLLECTION>(particle, tbox);
        } else if (t >= tbox[1]) {
            // front only
            return front->intersect<COLLECTION>(particle, tbox);
        }

        // both directions (start with front)
        const std::array<T, 2> t_front { tbox[0], t };
        auto hit = front->intersect<COLLECTION>(particle, t_front);
        if (hit.valid()) {
            if (hit.intersection <= t) {
                return hit;
            }
        }
        const std::array<T, 2> t_back { t, tbox[1] };
        return back->intersect<COLLECTION>(particle, t_back);
    }

    T planeSplit(std::vector<Tetrahedron<T>>& tets) const
    {
        const auto D = m_D;
        std::sort(std::execution::par_unseq, tets.begin(), tets.end(), [D](const auto& lh, const auto& rh) {
            const auto lc = lh.vertices()[0][D];
            const auto rc = rh.vertices()[0][D];
            return lc < rc;
        });
        const auto N = tets.size();
        T split;
        if (N % 2 == 1) {
            split = tets[N / 2].center()[D];
        } else {
            split = (tets[N / 2].center()[D] + tets[N / 2 - 1].center()[D]) / 2;
        }

        return split;
    }

    int figureOfMerit(const std::vector<Tetrahedron<T>>& tets, const T planesep) const
    {
        int fom = 0;
        int shared = 0;
        for (const auto& tet : tets) {
            const auto side = planeSide(tet, planesep, m_D);
            fom += side;
            if (side == 0) {
                shared++;
            }
        }
        return std::abs(fom) + shared;
    }
    static std::int_fast8_t planeSide(const Tetrahedron<T>& tet, const T plane, const std::uint_fast32_t D)
    {
        auto min = std::numeric_limits<T>::max();
        auto max = std::numeric_limits<T>::lowest();
        for (const auto& vert : tet) {
            min = std::min(min, vert[D]);
            max = std::max(max, vert[D]);
        }
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

    void item_iterator(std::vector<Tetrahedron<T>>& all) const
    {
        if (m_left) {
            m_left->item_iterator(all);
            m_right->item_iterator(all);
        } else {
            std::copy(m_tets.cbegin(), m_tets.cend(), std::back_inserter(all));
        }
    }
    void AABB_iterator(std::array<T, 6>& aabb) const
    {
        if (!m_left) {
            for (const auto& tet : m_tets) {
                for (const auto& v : tet) {
                    for (std::size_t i = 0; i < 3; ++i) {
                        aabb[i] = std::min(aabb[i], v[i]);
                        aabb[i + 3] = std::max(aabb[i + 3], v[i]);
                    }
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
    std::vector<Tetrahedron<T>> m_tets;
    std::unique_ptr<TetrahedalMeshKDTree<T>> m_left = nullptr;
    std::unique_ptr<TetrahedalMeshKDTree<T>> m_right = nullptr;
    // friend class TetrahedalMeshKDTree<T>;
};
}