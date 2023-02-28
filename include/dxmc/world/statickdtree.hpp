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

#include <array>
#include <concepts>
#include <execution>
#include <numeric>
#include <optional>
#include <vector>

namespace dxmc {

template <typename U, typename T>
concept StaticKDTreeType = requires(U u, Particle<T> p, std::array<T, 3> vec) {
                               Floating<T>;

                               u.translate(vec);
                               {
                                   u.intersect(p)
                                   } -> std::same_as<WorldIntersectionResult<T>>;
                               {
                                   u.center()
                                   } -> std::convertible_to<std::array<T, 3>>;
                               {
                                   u.AABB()
                                   } -> std::convertible_to<std::array<T, 6>>;
                           };

template <int Depth, Floating T, StaticKDTreeType<T> U>
class StaticKDTree {
public:
    StaticKDTree() { }
    StaticKDTree(std::vector<U>& data)
    {
        if (data.size() == 0)
            return;

        const auto [D, plane] = planeSplit(data);

        m_D = D;
        m_plane = plane;

        // moving object depending on plane splits
        std::vector<U> left, right;
        for (std::size_t i = 0; i < data.size(); ++i) {
            const auto side = planeSide(data[i], m_D, m_plane);
            if (side == -1) {
                left.push_back(data[i]);
            } else if (side == 1) {
                right.push_back(data[i]);
            } else {
                left.push_back(data[i]);
                right.push_back(data[i]);
            }
        }

        m_left = StaticKDTree<Depth - 1, T, U>(left);
        m_right = StaticKDTree<Depth - 1, T, U>(right);
    }
    void translate(const std::array<T, 3>& dir)
    {
        m_plane += dir[m_D];
        m_left.translate(dir);
        m_right.translate(dir);
    }
    std::array<T, 6> AABB() const
    {

        auto left = m_left.AABB();
        const auto right = m_right.AABB();
        for (std::size_t i = 0; i < 3; ++i) {
            left[i] = std::min(left[i], right[i]);
        }
        for (std::size_t i = 3; i < 6; ++i) {
            left[i] = std::max(left[i], right[i]);
        }
        return left;
    }
    KDTreeIntersectionResult<T, const U> intersect(const Particle<T>& particle, const std::array<T, 6>& aabb) const
    {
        const auto tbox = basicshape::AABB::intersectForwardInterval(particle, aabb);
        return tbox ? intersect(particle, *tbox) : KDTreeIntersectionResult<T, const U> {};
    }
    KDTreeIntersectionResult<T, const U> intersect(const Particle<T>& particle, const std::array<T, 2>& tbox) const
    {
        // test for parallell beam, if parallell we must test both sides.
        if (std::abs(particle.dir[m_D]) <= std::numeric_limits<T>::epsilon()) {
            const auto hit_left = m_left.intersect(particle, tbox);
            const auto hit_right = m_right.intersect(particle, tbox);
            if (hit_left.valid() && hit_right.valid())
                return hit_left.intersection < hit_right.intersection ? hit_left : hit_right;
            return hit_left.valid() ? hit_left : hit_right;
        }

        const auto [front, back] = particle.dir[m_D] > T { 0 } ? std::make_pair(&m_left, &m_right) : std::make_pair(&m_right, &m_left);

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
        if (hit.valid()) {
            if (hit.intersection <= t) {
                return hit;
            }
        }
        const std::array<T, 2> t_back { t, tbox[1] };
        return back->intersect(particle, t_back);
    }

protected:
    static std::pair<std::uint_fast32_t, T> planeSplit(const std::vector<U>& data)
    {
        // split where we find best separation between objects

        std::array<std::vector<std::pair<T, T>>, 3> segs;
        for (const auto& u : data) {
            const auto aabb = u.AABB();
            for (std::size_t i = 0; i < 3; ++i) {
                const auto seg = std::make_pair(aabb[i], aabb[i + 3]);
                bool added = false;
                std::size_t teller = 0;
                while (!added && teller != segs[i].size()) {
                    auto& r_seg = segs[i][teller];
                    if (r_seg.second > seg.first || r_seg.first > seg.second) {
                        r_seg.first = std::min(r_seg.first, seg.first);
                        r_seg.second = std::min(r_seg.second, seg.second);
                        added = true;
                    }
                    ++teller;
                }
                if (!added) {
                    segs[i].push_back(seg);
                }
            }
        }

        const std::array<std::size_t, 3> n_splits = { segs[0].size(), segs[1].size(), segs[2].size() };
        const auto n_splits_max = std::max(n_splits[0], std::max(n_splits[1], n_splits[2]));

        if (n_splits_max == 1) {
            const std::uint_fast32_t D = 0;
            const auto plane = std::nextafter(segs[D][0].second, std::numeric_limits<T>::max());
            return std::make_pair(D, plane);
        } else {
            // finding dim with min 2 splits and max extent
            std::uint_fast32_t D = 0;
            T extent = 0;
            for (std::uint_fast32_t i = 0; i < 3; ++i) {
                if (n_splits[i] == n_splits_max) {
                    T min = std::numeric_limits<T>::max();
                    T max = std::numeric_limits<T>::lowest();
                    for (const auto& [fi, se] : segs[i]) {
                        min = std::min(min, fi);
                        max = std::max(max, se);
                    }
                    const auto ex_cand = max - min;
                    if (ex_cand > extent) {
                        D = i;
                        extent = ex_cand;
                    }
                }
            }

            const auto ind_split = n_splits_max / 2;
            const auto plane = (segs[D][ind_split - 1].second + segs[D][ind_split].first) * T { 0.5 };
            return std::make_pair(D, plane);
        }
    }

    static int planeSide(const U& obj, std::uint_fast32_t D, const T planesep)
    {
        const auto& aabb = obj.AABB();
        const auto min = aabb[D];
        const auto max = aabb[D + 3];

        if (max < planesep)
            return -1;
        if (min > planesep)
            return 1;
        return 0;
    }

private:
    std::uint_fast32_t m_D = 0;
    T m_plane = 0;
    StaticKDTree<Depth - 1, T, U> m_left;
    StaticKDTree<Depth - 1, T, U> m_right;
};

template <Floating T, StaticKDTreeType<T> U>
class StaticKDTree<0, T, U> {
public:
    StaticKDTree() { }
    StaticKDTree(std::vector<U>& data)
        : m_data(data)
    {
    }
    void translate(const std::array<T, 3>& dir)
    {
        for (auto& u : m_data)
            u.translate(dir);
    }
    std::array<T, 6> AABB() const noexcept
    {

        std::array<T, 6> aabb {
            std::numeric_limits<T>::max(),
            std::numeric_limits<T>::max(),
            std::numeric_limits<T>::max(),
            std::numeric_limits<T>::lowest(),
            std::numeric_limits<T>::lowest(),
            std::numeric_limits<T>::lowest(),
        };
        for (const auto& u : m_data) {
            const auto& aabb_obj = u.AABB();
            for (std::size_t i = 0; i < 3; ++i) {
                aabb[i] = std::min(aabb[i], aabb_obj[i]);
            }
            for (std::size_t i = 3; i < 6; ++i) {
                aabb[i] = std::max(aabb[i], aabb_obj[i]);
            }
        }

        return aabb;
    }

    static constexpr std::size_t depth()
    {
        return 0;
    }
    KDTreeIntersectionResult<T, const U> intersect(const Particle<T>& particle, const std::array<T, 2>& tbox) const noexcept
    {
        KDTreeIntersectionResult<T, const U> res { .item = nullptr, .intersection = std::numeric_limits<T>::max() };

        for (const auto& u : m_data) {
            const auto t_cand = u.intersect(particle);
            if (t_cand.valid())
                if (t_cand.intersection < res.intersection) {
                    res.intersection = t_cand.intersection;
                    res.rayOriginIsInsideItem = t_cand.rayOriginIsInsideItem;
                    res.item = &u;
                }
        }
        return res;
    }

private:
    std::vector<U> m_data;
};
}