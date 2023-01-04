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
#include "dxmc/world/worlditembase.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <execution>
#include <iterator>
#include <memory>
#include <numeric>
#include <optional>
#include <vector>

namespace dxmc {

template <Floating T>
class KDTree {
    friend class KDTree<T>;

public:
    KDTree() {};
    KDTree(const std::vector<const WorldItemBase<T>*>& items, const std::size_t max_depth = 8)
    {
        if (items.size() < 2) {
            for (const auto& item : items)
                m_items.push_back(item);
            return;
        }
        // finding aabb
        std::array<T, 6> aabb {
            std::numeric_limits<T>::max(),
            std::numeric_limits<T>::max(),
            std::numeric_limits<T>::max(),
            std::numeric_limits<T>::lowest(),
            std::numeric_limits<T>::lowest(),
            std::numeric_limits<T>::lowest(),
        };
        for (const auto& item : items) {
            std::array<T, 6> aabb_tri;
            aabb_tri = item->AABB();

            for (std::size_t i = 0; i < 3; ++i) {
                aabb[i] = std::min(aabb[i], aabb_tri[i]);
            }
            for (std::size_t i = 3; i < 6; ++i) {
                aabb[i] = std::max(aabb[i], aabb_tri[i]);
            }
        }
        const std::array<T, 3> extent { aabb[3] - aabb[0], aabb[4] - aabb[1], aabb[5] - aabb[2] };

        m_D = vectormath::argmax3<std::uint_fast32_t, T>(extent);

        const auto split = planeSplit(items);

        const auto fom = figureOfMerit(items, split);

        if (fom == items.size() || max_depth <= 1 || items.size() <= 1) {
            m_items = items;
        } else {
            m_plane = split;
            std::vector<const WorldItemBase<T>*> left;
            std::vector<const WorldItemBase<T>*> right;
            for (const auto& item : items) {
                const auto side = planeSide(item, m_plane, m_D);
                if (side <= 0)
                    left.push_back(item);
                if (side >= 0)
                    right.push_back(item);
            }
            m_left = std::make_unique<KDTree<T>>(left, max_depth - 1);
            m_right = std::make_unique<KDTree<T>>(right, max_depth - 1);
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
    std::vector<const WorldItemBase<T>*> items() const
    {
        std::vector<const WorldItemBase<T>*> all;
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
            std::for_each(std::execution::par_unseq, m_items.begin(), m_items.end(), [&](auto& tri) {
                tri->translate(dist);
            });
        }
    }

    IntersectionResult<T> intersect(const Particle<T>& particle, const std::array<T, 6>& aabb) const
    {
        const auto& inter = intersectAABB(particle, aabb);
        return inter ? intersect(particle, *inter) : IntersectionResult<T> {};
    }

protected:
    // delete this?
    std::optional<T> intersect(const Particle<T>& particle, const std::array<T, 6>& aabb, const std::array<T, 2>& tbox) const
    {
        const auto& inter = intersectAABB(particle, aabb);
        if (inter) {
            const std::array<T, 2> tbox_min { std::min((*inter)[0], tbox[0]), std::min((*inter)[1], tbox[1]) };
            return intersect(particle, tbox_min);
        } else
            return std::nullopt;
    }

    IntersectionResult<T> intersect(const Particle<T>& particle, const std::array<T, 2>& tbox) const
    {
        if (!m_left) { // this is a leaf
            // intersect triangles between tbox and return;

            IntersectionResult<T> res { .item = nullptr, .intersection = std::numeric_limits<T>::max() };
            for (const auto& item : m_items) {
                const auto t_cand = item->intersect(particle);
                if (t_cand.item) {
                    if (greaterOrEqual(t_cand.intersection, tbox[0]) && lessOrEqual(t_cand.intersection, tbox[1])) {
                        if (t_cand.intersection < res.intersection) {
                            res = t_cand;
                        }
                    }
                }
            }
            return res;
            // return greaterOrEqual(t, tbox[0]) && lessOrEqual(t, tbox[1]) ? std::make_optional(t) : std::nullopt;
        }

        // test for parallell beam
        if (std::abs(particle.dir[m_D]) <= std::numeric_limits<T>::epsilon()) {
            const auto hit_left = m_left->intersect(particle, tbox);
            const auto hit_right = m_right->intersect(particle, tbox);
            if (hit_left.item && hit_right.item)
                return hit_left.intersection < hit_right.intersection ? hit_left : hit_right;
            if (!hit_left.item)
                return hit_right;
            return hit_left;
        }

        const KDTree<T>* const front = particle.dir[m_D] > T { 0 } ? m_left.get() : m_right.get();
        const KDTree<T>* const back = particle.dir[m_D] > T { 0 } ? m_right.get() : m_left.get();

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
        const auto hit = front->intersect(particle, t_front);
        if (hit.item) {
            if (hit.intersection <= t) {
                return hit;
            }
        }
        const std::array<T, 2> t_back { t, tbox[1] };
        return back->intersect(particle, t_back);
    }
    static bool inPlaneBox(const std::array<T, 3>& pos, const std::array<T, 6>& aabb, const std::uint_fast8_t axis)
    {
        bool inside = true;
        for (std::uint_fast8_t i = 0; i < 3; ++i) {
            if (i != axis) {
                inside = inside && aabb[i] >= pos[i] && pos[i] <= aabb[i + 3];
            }
        }
        return inside;
    }
    static std::optional<std::array<T, 2>> intersectAABB(const Particle<T>& p, const std::array<T, 6>& aabb)
    {
        std::array<T, 2> t {
            T { 0 },
            std::numeric_limits<T>::max()
        };
        for (std::uint_fast8_t i = 0; i < 3; i++) {
            if (std::abs(p.dir[i]) > std::numeric_limits<T>::epsilon()) {
                const auto d_inv = T { 1 } / p.dir[i];

                if (p.pos[i] < aabb[i]) {
                    const auto t_cand = (aabb[i] - p.pos[i]) * d_inv;
                    const auto npos = vectormath::add(p.pos, vectormath::scale(p.dir, t_cand));
                    if (inPlaneBox(npos, aabb, i)) {
                        t[0] = std::max(t[0], t_cand);
                        t[1] = std::min(t[1], t_cand);
                    }
                } else if (p.pos[i] > aabb[i + 3]) {
                    const auto t_cand = (aabb[i + 3] - p.pos[i]) * d_inv;
                    const auto npos = vectormath::add(p.pos, vectormath::scale(p.dir, t_cand));
                    if (inPlaneBox(npos, aabb, i)) {
                        t[0] = std::max(t[0], t_cand);
                        t[1] = std::min(t[1], t_cand);
                    }
                } else {
                    if (p.dir[i] > T { 0 }) {
                        const auto t_cand = (aabb[i + 3] - p.pos[i]) * d_inv;
                        const auto npos = vectormath::add(p.pos, vectormath::scale(p.dir, t_cand));
                        if (inPlaneBox(npos, aabb, i)) {
                            t[0] = std::max(t[0], t_cand);
                            t[1] = std::min(t[1], t_cand);
                        }
                    } else {
                        const auto t_cand = (aabb[i] - p.pos[i]) * d_inv;
                        const auto npos = vectormath::add(p.pos, vectormath::scale(p.dir, t_cand));
                        if (inPlaneBox(npos, aabb, i)) {
                            t[0] = std::max(t[0], t_cand);
                            t[1] = std::min(t[1], t_cand);
                        }
                    }
                }
            }
        }
        return t[0] > t[1] ? std::nullopt : std::make_optional(t);
    }

    T planeSplit(const std::vector<const WorldItemBase<T>*>& items) const
    {
        const auto N = items.size();
        std::vector<T> vals;
        vals.reserve(N);

        for (const auto& item : items) {
            const auto v = item->center()[m_D];
            vals.push_back(v);
        }

        std::sort(vals.begin(), vals.end());

        if (N % 2 == 1) {
            return vals[N / 2];
        } else {
            return (vals[N / 2] + vals[N / 2 - 1]) * T { 0.5 };
        }
    }

    int figureOfMerit(const std::vector<const WorldItemBase<T>*>& items, const T planesep) const
    {
        int fom = 0;
        int shared = 0;
        for (const auto& item : items) {
            const auto side = planeSide(item, planesep, m_D);
            fom += side;
            if (side == 0) {
                shared++;
            }
        }
        return std::abs(fom) + shared;
    }
    static int planeSide(const WorldItemBase<T>* triangle, const T plane, const unsigned int D)
    {
        T max = std::numeric_limits<T>::lowest();
        T min = std::numeric_limits<T>::max();

        const auto& aabb = triangle->AABB();
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
    void item_iterator(std::vector<const WorldItemBase<T>*>& all) const
    {
        if (m_left) {
            m_left->item_iterator(all);
            m_right->item_iterator(all);
        } else {
            std::copy(m_items.cbegin(), m_items.cend(), std::back_inserter(all));
        }
    }
    void AABB_iterator(std::array<T, 6>& aabb) const
    {
        if (!m_left) {
            for (const auto& item : m_items) {
                const auto aabb_tri = item->AABB();
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
    std::vector<const WorldItemBase<T>*> m_items;

    std::unique_ptr<KDTree<T>> m_left = nullptr;
    std::unique_ptr<KDTree<T>> m_right = nullptr;
};

/*
template <Floating T, KDTreeType<T> U>
class KDTree {
public:
    // KDTree() {};
    KDTree(std::vector<U>& triangles, const std::size_t max_depth = 8)
    {
        m_node = KDTreeNode<T, U>(triangles, max_depth);
        m_aabb = m_node.AABB();
    }
    std::size_t depth() const
    {
        return m_node.depth();
    }

    void translate(const std::array<T, 3>& vec)
    {
        m_node.translate(vec);
    }

    std::array<T, 3> center() const
    {
        std::array<T, 3> center { (m_aabb[0] + m_aabb[3]) / 2, (m_aabb[1] + m_aabb[4]) / 2, (m_aabb[2] + m_aabb[5]) / 2 };
        return center;
    }
    std::array<T, 6> AABB() const
    {
        return m_node.AABB();
    }

    std::optional<T> intersect(const Particle<T>& particle) const
    {
        return m_node.intersect(particle, m_aabb);
    }

private:
    std::array<T, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };
    KDTreeNode<T, U> m_node;
};
*/
}