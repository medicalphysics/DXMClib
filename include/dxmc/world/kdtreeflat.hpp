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

#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/basicshapes/aabb.hpp"
#include "dxmc/world/kdtreeintersectionresult.hpp"
#include "dxmc/world/visualizationintersectionresult.hpp"
#include "dxmc/world/worlditems/worlditemtype.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <concepts>
#include <execution>
#include <iterator>
#include <memory>
#include <numeric>
#include <optional>
#include <variant>
#include <vector>

namespace dxmc {

template <WorldItemType... Us>
class KDTreeFlat {
public:
    KDTreeFlat() {};
    KDTreeFlat(std::vector<std::variant<Us...>*>& items, const std::size_t max_depth = 8)
    {
        setData(items, max_depth);
    }

    void setData(std::vector<std::variant<Us...>*>& items, const std::size_t max_depth = 8)
    {
        m_nodes.reserve(items.size()); // left nodes
        auto split_dim = splitAxis(items);
        auto split_val = splitPlane(items, split_dim);
        build(items, m_nodes, max_depth, split_val, split_dim);
    };

    std::array<double, 6> calculateAABB() const
    {
        std::array<double, 6> aabb {
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
        };
        for (const auto& node : m_nodes) {
            if (node.dimension == 4) {
                auto item = node.data.element;
                const auto aabb_tri = std::visit([](const auto& it) { return it.AABB(); }, *item);
                for (int i = 0; i < 3; ++i) {
                    aabb[i] = std::min(aabb[i], aabb_tri[i]);
                }
                for (int i = 3; i < 6; ++i) {
                    aabb[i] = std::max(aabb[i], aabb_tri[i]);
                }
            }
        }
        return aabb;
    }

    KDTreeIntersectionResult<std::variant<Us...>> intersect(const ParticleType auto& particle, const std::array<double, 6>& aabb)
    {
        const auto& inter = basicshape::AABB::intersectForwardInterval(particle, aabb);
        return inter ? intersect(particle, *inter) : KDTreeIntersectionResult<std::variant<Us...>> {};
    }

    KDTreeIntersectionResult<std::variant<Us...>> intersect(const ParticleType auto& particle, const std::array<double, 2>& tbox)
    {
        if (!m_left) { // this is a leaf
            // intersect triangles between tbox and return;

            KDTreeIntersectionResult<std::variant<Us...>> res = { .item = nullptr, .intersection = std::numeric_limits<double>::max() };
            for (auto& item : m_items) {
                // const auto t_cand = item->intersect(particle);
                const auto t_cand = std::visit([&particle](const auto& it) { return it.intersect(particle); }, *item);
                if (t_cand.valid()) {
                    const auto border_delta = t_cand.rayOriginIsInsideItem ? -GEOMETRIC_ERROR() : GEOMETRIC_ERROR();
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
        if (std::abs(particle.dir[m_D]) <= std::numeric_limits<double>::epsilon()) {
            const auto hit_left = m_left->intersect(particle, tbox);
            const auto hit_right = m_right->intersect(particle, tbox);
            if (hit_left.item && hit_right.item)
                return hit_left.intersection < hit_right.intersection ? hit_left : hit_right;
            if (!hit_left.item)
                return hit_right;
            return hit_left;
        }

        auto front = particle.dir[m_D] > 0 ? m_left.get() : m_right.get();
        auto back = particle.dir[m_D] > 0 ? m_right.get() : m_left.get();

        const auto t = (m_plane - particle.pos[m_D]) / particle.dir[m_D];

        if (t <= tbox[0]) {
            // back only
            return back->intersect(particle, tbox);
        } else if (t >= tbox[1]) {
            // front only
            return front->intersect(particle, tbox);
        }

        // both directions (start with front)
        const std::array<double, 2> t_front { tbox[0], t };
        const auto hit = front->intersect(particle, t_front);
        if (hit.item) {
            if (hit.intersection <= t) {
                return hit;
            }
        }
        const std::array<double, 2> t_back { t, tbox[1] };
        return back->intersect(particle, t_back);
    }

protected:
    static std::uint8_t splitAxis(const std::vector<std::variant<Us...>*>& items)
    {
        // finding aabb
        std::array<double, 6> aabb {
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
        };
        for (const auto& item : items) {
            const auto aabb_tri = std::visit([](const auto& it) { return it.AABB(); }, *item);
            for (std::size_t i = 0; i < 3; ++i) {
                aabb[i] = std::min(aabb[i], aabb_tri[i]);
            }
            for (std::size_t i = 3; i < 6; ++i) {
                aabb[i] = std::max(aabb[i], aabb_tri[i]);
            }
        }
        const std::array<double, 3> extent { aabb[3] - aabb[0], aabb[4] - aabb[1], aabb[5] - aabb[2] };
        return vectormath::argmax3<std::uint8_t>(extent);
    }

    static double splitPlane(const std::vector<std::variant<Us...>*>& items, const std::uint8_t dim)
    {
        const auto N = items.size();
        std::vector<double> vals;
        vals.reserve(N);

        for (const auto& item : items) {
            const auto v = std::visit([](const auto& it) { return it.center(); }, *item);
            vals.push_back(v[dim]);
        }

        std::sort(vals.begin(), vals.end());

        if (N % 2 == 1) {
            return vals[N / 2];
        } else {
            return (vals[N / 2] + vals[N / 2 - 1]) * 0.5;
        }
    }

    static int figureOfMerit(const std::vector<std::variant<Us...>*>& items, const double planesep, const std::uint8_t dim)
    {
        int fom = 0;
        int shared = 0;
        for (const auto& item : items) {
            const auto side = planeSide(item, planesep, dim);
            fom += side;
            if (side == 0) {
                shared++;
            }
        }
        return std::abs(fom) + shared;
    }

    static int planeSide(std::variant<Us...>* item, const double plane, const std::uint8_t D)
    {
        auto max = std::numeric_limits<double>::lowest();
        auto min = std::numeric_limits<double>::max();

        const auto aabb = std::visit([](const auto& it) { return it.AABB(); }, *item);

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
    constexpr static double epsilon()
    {
        // Huristic epsilon
        return 11 * std::numeric_limits<double>::epsilon();
    }

    constexpr static bool lessOrEqual(double a, double b)
    {
        return a - b <= epsilon() * a;
    }
    constexpr static bool greaterOrEqual(double a, double b)
    {
        return b - a <= epsilon() * a;
    }

private:
    struct Node {
        union data {
            double split;
            std::variant<Us...>* element;
        };
        // if dim == 4 this is a leaf and index n remaining elements, else this is a node and index is offset to right branch, left branch is next node
        std::uint32_t index = 0;
        std::uint8_t dimension = 4;
    };

    static void build(const std::vector<std::variant<Us...>*>& items, std::vector<Node>& nodes, const std::size_t depth, double split_val, std::uint8_t split_dim)
    {

        const int fom = figureOfMerit(items, split_val, split_dim);

        if (depth == 0 || fom == items.size() || items.size() < 2) {
            for (std::uint32_t i = items.size(); i > 0; --i) {
                Node node;
                node.data.element = items[i];
                node.dimension = 4;
                node.index = i - 1;
                nodes.push_back(node);
            }
        } else {
            // construct node
            const std::uint32_t current_Idx = static_cast<std::uint32_t>(nodes.size());
            Node node;
            node.dimension = split_dim;
            node.data.split = split_value;
            nodes.push_back(node);

            std::vector<std::variant<Us...>*> left, right;
            for (auto item : items) {
                const auto side = planeSide(item, split_val, split_dim);
                if (side <= 0)
                    left.push_back(item);
                if (side >= 0)
                    right.push_back(item);
            }

            // construct left side
            split_dim = splitAxis(left);
            split_val = splitPlane(left, split_dim);
            build(left, nodes, depth - 1, split_val, split_dim);
            // adding offset index to parent
            nodes[current_Idx].index = static_cast<std::uint32_t>(nodes.size()) - current_Idx;

            // construct right side
            split_dim = splitAxis(right);
            split_val = splitPlane(right, split_dim);
            build(right, nodes, depth - 1, split_val, split_dim);
        }
    };

    std::vector<Node> m_nodes;
};

}