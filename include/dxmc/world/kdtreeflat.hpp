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

    void setData(std::vector<std::variant<Us...>*>& items, const std::size_t max_depth = 8) {

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

        return aabb;
    }

    KDTreeIntersectionResult<std::variant<Us...>> intersect(const ParticleType auto& particle, const std::array<double, 6>& aabb)
    {
        const auto& inter = basicshape::AABB::intersectForwardInterval(particle, aabb);
        return inter ? intersect(particle, *inter, 0) : KDTreeIntersectionResult<std::variant<Us...>> {};
    }

    KDTreeIntersectionResult<std::variant<Us...>> intersect(const ParticleType auto& particle, const std::array<double, 2>& tbox, std::uint32_t node_idx)
    {
        // if leaf
        if (m_nodes[node_idx].dim == 4) {
            KDTreeIntersectionResult<std::variant<Us...>> res = { .item = nullptr, .intersection = std::numeric_limits<double>::max() };
            while (m_nodes[node_idx].index >= 0) {
                const auto item = nm_nodes[node_idx].data.element;
                const auto t_cand = std::visit([&particle](const auto& it) { return it.intersect(particle); }, *item);
                if (t_cand.intersection + border_delta < res.intersection) {
                    if (tbox[0] <= t_cand.intersection && t_cand.intersection <= tbox[1]) {
                        res.intersection = t_cand.intersection;
                        res.item = item;
                        res.rayOriginIsInsideItem = t_cand.rayOriginIsInsideItem;
                    }
                }
            }
            return res;
        }

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
    void build(std::vector<std::uint32_t>& indices, int max_depth = 8)
    {
        struct NodeTemplate {
            Node node;
            std::vector<std::uint32_t> indices;
        };

        std::vector<NodeTemplate> nodes(1);
        nodes.reserve(indices.size());

        nodes[0].indices = indices;
        std::size_t currentNodeIdx = 0;

        while (currentNodeIdx < m_nodes.size()) {
            auto& cnode = nodes[currentNodeIdx].node;
            auto& cind = nodes[currentNodeIdx].indices;
            auto split_dim = splitAxis(cind);
            auto split_val = splitPlane(cind, split_dim);
            auto fom = figureOfMerit(cind, split_dim, split_val);
            if (fom != cind.size()) {
                // leaf
                cnode.setLeaf();
                cnode.split_nelements.nelements = static_cast<std::uint32_t>(cind.size());
                cnode.setOffset(static_cast<std::uint32_t>(m_indices.size());
                for (auto i :cind)
                    m_indices.push_back(i);
                cind.clear();
                cind.shrink_to_fit();
            } else {
                // branch
                cnode.setDim(split_dim);
                cnode.split_nelements.split = split_val;
                cnode.setOffset(nodes.size());
                NodeTemplate left, right;
                populate left right;
                nodes.push_back(left);
                nodes.push_back(right);
            }
            ++currentNodeIdx;
        }

        // copy nodes;
    }

    std::uint32_t splitAxis(const std::vector<std::uint32_t>& indices)
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

        for (auto idx : indices) {
            auto* item = m_items[idx];
            const auto aabb_tri = std::visit([](const auto& it) { return it.AABB(); }, *item);
            for (std::size_t i = 0; i < 3; ++i) {
                aabb[i] = std::min(aabb[i], aabb_tri[i]);
            }
            for (std::size_t i = 3; i < 6; ++i) {
                aabb[i] = std::max(aabb[i], aabb_tri[i]);
            }
        }
        const std::array<double, 3> extent { aabb[3] - aabb[0], aabb[4] - aabb[1], aabb[5] - aabb[2] };
        return vectormath::argmax3<std::uint32_t>(extent);
    }

    float splitPlane(const std::vector<std::uint32_t>& indices, const std::uint32_t dim)
    {
        const auto N = items.size();
        std::vector<float> vals;
        vals.reserve(N);

        for (auto idx : indices) {
            auto* item = m_items[idx];
            const auto v = std::visit([](const auto& it) { return it.center(); }, *item);
            vals.push_back(static_cast<float>(v[dim]));
        }

        std::sort(vals.begin(), vals.end());

        if (N % 2 == 1) {
            return vals[N / 2];
        } else {
            return (vals[N / 2] + vals[N / 2 - 1]) / 2;
        }
    }

    int planeSide(std::variant<Us...>* item, const std::uint32_t D, const float plane)
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

    int figureOfMerit(const std::vector<std::uint32_t>& indices, const std::uint32_t dim, const float planesep)
    {
        int fom = 0;
        int shared = 0;
        for (auto idx : indices) {
            auto* item = m_items[idx];
            const auto side = planeSide(item, planesep, dim);
            fom += side;
            if (side == 0) {
                shared++;
            }
        }
        return std::abs(fom) + shared;
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
        struct {
            std::uint32_t dim : 2; // dimensjon for branch
            std::uint32_t offset : 29; // offset to first child (branch) or to first element (leaf)
            std::uint32_t flag : 1; // Node is leaf (1) or branch (0)
        } dim_offset_flag;

        union {
            float split = 0; // Union split is for branches
            std::uint32_t nelements; // Union nelements is for number of items
        } split_nelements;

        Node()
        {
            dim_offset_flag.dim = 0;
            dim_offset_flag.offset = 0;
            dim_offset_flag.flag = 0;
        }

        std::uint32_t dim()
        {
            return dim_offset_flag.dim;
        }
        void setDim(std::uint32_t dim)
        {
            dim_offset_flag.dim = dim;
        }

        void setLeaf()
        {
            dim_offset_flag.flag = std::uint32_t { 1 };
        }
        std::uint32_t isLeaf()
        {
            return dim_offset_flag.flag;
        }

        void setOffset(std::uint32_t offset)
        {
            dim_offset_flag.offset = offset;
        }
        std::uint32_t offset()
        {
            return dim_offset_flag.offset;
        }
    };

    std::vector<std::uint32_t> m_indices;
    std::vector<std::variant<Us...>*> m_items;
    std::vector<Node> m_nodes;
};
}