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
#include "dxmc/world/basicshapes/aabb.hpp"
#include "dxmc/world/kdtreeintersectionresult.hpp"
#include "dxmc/world/visualizationintersectionresult.hpp"
#include "dxmc/world/worlditems/worlditembase.hpp"

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

public:
    KDTree() {};
    KDTree(std::vector<WorldItemBase<T>*>& items, const std::size_t max_depth = 8)
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
            std::vector<WorldItemBase<T>*> left;
            std::vector<WorldItemBase<T>*> right;
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
    std::vector<WorldItemBase<T>*> items()
    {
        std::vector<WorldItemBase<T>*> all;
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

    KDTreeIntersectionResult<T, WorldItemBase<T>> intersect(const Particle<T>& particle, const std::array<T, 6>& aabb)
    {

        const auto& inter = basicshape::AABB::intersectForwardInterval(particle, aabb);
        return inter ? intersect(particle, *inter) : KDTreeIntersectionResult<T, WorldItemBase<T>> {};
    }

    KDTreeIntersectionResult<T, WorldItemBase<T>> intersect(const Particle<T>& particle, const std::array<T, 2>& tbox)
    {
        if (!m_left) { // this is a leaf
            // intersect triangles between tbox and return;

            KDTreeIntersectionResult<T, WorldItemBase<T>> res = { .item = nullptr, .intersection = std::numeric_limits<T>::max() };
            for (auto& item : m_items) {
                const auto t_cand = item->intersect(particle);
                if (t_cand.valid()) {
                    const T border_delta = t_cand.rayOriginIsInsideItem ? T { -1E-6 } : T { 1E-6 };
                    if (t_cand.intersection + border_delta < res.intersection) {
                        if (tbox[0] <= t_cand.intersection && t_cand.intersection <= tbox[1]) {
                            res.intersection = t_cand.intersection;
                            res.item = item;
                            res.rayOriginIsInsideItem = t_cand.rayOriginIsInsideItem;
                        }
                    }
                }
            }
            return res;
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

        KDTree<T>* const front = particle.dir[m_D] > T { 0 } ? m_left.get() : m_right.get();
        KDTree<T>* const back = particle.dir[m_D] > T { 0 } ? m_right.get() : m_left.get();

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

    VisualizationIntersectionResult<T, WorldItemBase<T>> intersectVisualization(const Particle<T>& particle, const std::array<T, 6>& aabb)
    {
        const auto& inter = basicshape::AABB::intersectForwardInterval(particle, aabb);
        return inter ? intersectVisualization(particle, *inter) : VisualizationIntersectionResult<T, WorldItemBase<T>> {};
    }

    VisualizationIntersectionResult<T, WorldItemBase<T>> intersectVisualization(const Particle<T>& particle, const std::array<T, 2>& tbox)
    {
        if (!m_left) { // this is a leaf
            // intersect triangles between tbox and return;

            VisualizationIntersectionResult<T, WorldItemBase<T>> res;
            res.intersection = std::numeric_limits<T>::max();
            for (auto& item : m_items) {
                const auto t_cand = item->intersectVisualization(particle);
                if (t_cand.valid()) {
                    if (t_cand.intersection < res.intersection) {
                        if (tbox[0] <= t_cand.intersection && t_cand.intersection <= tbox[1]) {
                            res = t_cand;
                            res.item = item;
                        }
                    }
                }
            }
            return res;
        }

        // test for parallell beam
        if (std::abs(particle.dir[m_D]) <= std::numeric_limits<T>::epsilon()) {
            const auto hit_left = m_left->intersectVisualization(particle, tbox);
            const auto hit_right = m_right->intersectVisualization(particle, tbox);
            if (hit_left.item && hit_right.item)
                return hit_left.intersection < hit_right.intersection ? hit_left : hit_right;
            if (!hit_left.item)
                return hit_right;
            return hit_left;
        }

        KDTree<T>* const front = particle.dir[m_D] > T { 0 } ? m_left.get() : m_right.get();
        KDTree<T>* const back = particle.dir[m_D] > T { 0 } ? m_right.get() : m_left.get();

        const auto t = (m_plane - particle.pos[m_D]) / particle.dir[m_D];

        if (t <= tbox[0]) {
            // back only
            return back->intersectVisualization(particle, tbox);
        } else if (t >= tbox[1]) {
            // front only
            return front->intersectVisualization(particle, tbox);
        }

        // both directions (start with front)
        const std::array<T, 2> t_front { tbox[0], t };
        const auto hit = front->intersectVisualization(particle, t_front);
        if (hit.item) {
            if (hit.intersection <= t) {
                return hit;
            }
        }
        const std::array<T, 2> t_back { t, tbox[1] };
        return back->intersectVisualization(particle, t_back);
    }

protected:
    T planeSplit(const std::vector<WorldItemBase<T>*>& items) const
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

    int figureOfMerit(const std::vector<WorldItemBase<T>*>& items, const T planesep) const
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
    static int planeSide(WorldItemBase<T>* triangle, const T plane, const unsigned int D)
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
    std::vector<WorldItemBase<T>*> m_items;

    std::unique_ptr<KDTree<T>> m_left = nullptr;
    std::unique_ptr<KDTree<T>> m_right = nullptr;
};

}