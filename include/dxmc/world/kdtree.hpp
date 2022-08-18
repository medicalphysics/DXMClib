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
#include <cmath>
#include <concepts>
#include <execution>
#include <iterator>
#include <memory>
#include <numeric>
#include <optional>
#include <type_traits>
#include <vector>

namespace dxmc {

template <typename U>
concept KDTreeType = requires(U u, Particle<typename U::Type> p, std::array<typename U::Type, 3> vec)
{
    typename U::Type;
    Floating<typename U::Type>;
    u <=> u;
    u.translate(vec);
    {
        u.intersect(p)
        } -> std::same_as<std::optional<typename U::Type>>;

    {
        u.center()
        } -> std::same_as<std::array<typename U::Type, 3>>;

    {
        u.AABB()
        } -> std::same_as<std::array<typename U::Type, 6>>;
}
|| requires(U u, Particle<typename std::remove_pointer<U>::type::Type> p, std::array<typename std::remove_pointer<U>::type::Type, 3> vec)
{
    std::remove_pointer<U>::type::Type;
    Floating<typename std::remove_pointer<U>::type::Type>;
    *u <=> *u;
    u->translate(vec);
    {
        u->intersect(p)
        } -> std::same_as<std::optional<typename std::remove_pointer<U>::type::Type>>;

    {
        u->center()
        } -> std::same_as<std::array<typename std::remove_pointer<U>::type::Type, 3>>;

    {
        u->AABB()
        } -> std::same_as<std::array<typename std::remove_pointer<U>::type::Type, 6>>;
};

template <KDTreeType U>
class KDTree {
public:
    using T = typename std::remove_pointer<U>::type::Type;
    using Type = T;
    KDTree() {};
    KDTree(std::vector<U>& triangles, const std::size_t max_depth = 8)
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
            if constexpr (std::is_pointer<U>::value) {
                aabb_tri = tri->AABB();
            } else {
                aabb_tri = tri.AABB();
            }
            for (std::size_t i = 0; i < 3; ++i) {
                aabb[i] = std::min(aabb[i], aabb_tri[i]);
            }
            for (std::size_t i = 3; i < 6; ++i) {
                aabb[i] = std::max(aabb[i], aabb_tri[i]);
            }
        }
        const std::array<T, 3> extent { aabb[3] - aabb[0], aabb[4] - aabb[1], aabb[5] - aabb[2] };

        m_D = vectormath::argmax3<unsigned int, T>(extent);

        std::sort(triangles.begin(), triangles.end(), [&](const auto& lh, const auto& rh) {
            if constexpr (std::is_pointer<U>::value) {
                const auto lf_val = lh->center();
                const auto rh_val = rh->center();
                return lf_val[m_D] < rh_val[m_D];
            } else {
                const auto lf_val = lh.center();
                const auto rh_val = rh.center();
                return lf_val[m_D] < rh_val[m_D];
            }
        });

        const auto mean = planeSplit(triangles);

        const auto fom = figureOfMerit(triangles, mean);

        if (fom == triangles.size() || fom == 0 || max_depth <= 1) {
            m_triangles = triangles;
        } else {
            m_plane = mean;
            std::vector<U> left;
            std::vector<U> right;
            for (const auto& triangle : triangles) {
                const auto side = planeSide(triangle, m_plane, m_D);
                if (side <= 0)
                    left.push_back(triangle);
                if (side >= 0)
                    right.push_back(triangle);
            }
            m_left = std::make_unique<KDTree<U>>(left, max_depth - 1);
            m_right = std::make_unique<KDTree<U>>(right, max_depth - 1);
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

    std::optional<T> intersect(const Particle<T>& particle, const std::array<T, 6>& aabb, const std::array<T, 2>& tbox) const
    {
        const auto& inter = intersectAABB(particle, aabb);
        if (inter) {
            const std::array<T, 2> tbox_min { std::min((*inter)[0], tbox[0]), std::min((*inter)[1], tbox[1]) };
            return intersect(particle, tbox_min);
        } else
            return std::nullopt;
    }

    std::optional<T> intersect(const Particle<T>& particle, const std::array<T, 6>& aabb) const
    {
        const auto& inter = intersectAABB(particle, aabb);
        return inter ? intersect(particle, *inter) : std::nullopt;
    }
    std::optional<T> intersect(const Particle<T>& particle, const std::array<T, 2>& tbox) const
    {
        if (!m_left) { // this is a leaf
            // intersect triangles between tbox and return;
            T t = std::numeric_limits<T>::max();
            for (const U& triangle : m_triangles) {
                if constexpr (std::is_pointer<U>::value) {
                    const auto t_cand = triangle->intersect(particle);
                    if (t_cand)
                        t = std::min(t, *t_cand);
                } else {
                    const auto t_cand = triangle.intersect<0>(particle);
                    if (t_cand)
                        t = std::min(t, *t_cand);
                }
            }
            return t <= tbox[1] && t >= tbox[0] ? std::make_optional(t) : std::nullopt;
        }

        // test for parallell beam
        if (std::abs(particle.dir[m_D]) <= std::numeric_limits<T>::epsilon()) {
            const auto hit_left = m_left->intersect(particle, tbox);
            const auto hit_right = m_right->intersect(particle, tbox);
            if (hit_left && hit_right)
                return std::min(hit_left, hit_right);
            if (!hit_left)
                return hit_right;
            return hit_left;
        }

        const KDTree<U>* const front = particle.dir[m_D] > T { 0 } ? m_left.get() : m_right.get();
        const KDTree<U>* const back = particle.dir[m_D] > T { 0 } ? m_right.get() : m_left.get();

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
        if (hit) {
            if (*hit <= t) {
                return hit;
            }
        }
        const std::array<T, 2> t_back { t, tbox[1] };
        return back->intersect(particle, t_back);
    }

protected:
    static std::optional<std::array<T, 2>> intersectAABB(const Particle<T>& p, const std::array<T, 6>& aabb)
    {
        std::array<T, 2> t {
            std::numeric_limits<T>::lowest(),
            std::numeric_limits<T>::max()
        };
        for (std::size_t i = 0; i < 3; i++) {
            if (std::abs(p.dir[i]) > std::numeric_limits<T>::epsilon()) {
                const auto d_inv = T { 1 } / p.dir[i];
                const auto t1 = (aabb[i] - p.pos[i]) * d_inv;
                const auto t2 = (aabb[i + 3] - p.pos[i]) * d_inv;
                const auto t_min_cand = std::min(t1, t2);
                const auto t_max_cand = std::max(t1, t2);
                t[0] = std::max(t[0], t_min_cand);
                t[1] = std::min(t[1], t_max_cand);
            }
        }
        return t[0] > t[1] ? std::nullopt : std::make_optional(t);
    }

    T planeSplit(const std::vector<U>& triangles) const
    {
        const auto N = triangles.size();
        if constexpr (std::is_pointer<U>::value) {
            if (N % 2 == 1) {
                return triangles[N / 2]->center()[m_D];
            } else {
                return (triangles[N / 2]->center()[m_D] + triangles[N / 2 - 1]->center()[m_D]) / 2;
            }
        } else {
            if (N % 2 == 1) {
                return triangles[N / 2].center()[m_D];
            } else {
                return (triangles[N / 2].center()[m_D] + triangles[N / 2 - 1].center()[m_D]) / 2;
            }
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
    static int planeSide(const U& triangle, const T plane, const unsigned int D)
    {
        T max = std::numeric_limits<T>::lowest();
        T min = std::numeric_limits<T>::max();
        constexpr T epsilon = std::numeric_limits<T>::epsilon();

        if constexpr (std::is_pointer<U>::value) {
            const auto& aabb = triangle->AABB();
            min = aabb[D];
            max = aabb[D + 3];
        } else {
            const auto& aabb = triangle.AABB();
            min = aabb[D];
            max = aabb[D + 3];
        }
        if (max - plane <= -epsilon)
            return -1;
        if (min - plane >= epsilon)
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
                if constexpr (std::is_pointer<U>::value) {
                    const auto aabb_tri = tri->AABB();
                    for (std::size_t i = 0; i < 3; ++i) {
                        aabb[i] = std::min(aabb[i], aabb_tri[i]);
                    }
                    for (std::size_t i = 3; i < 6; ++i) {
                        aabb[i] = std::max(aabb[i], aabb_tri[i]);
                    }
                } else {
                    const auto aabb_tri = tri.AABB();
                    for (std::size_t i = 0; i < 3; ++i) {
                        aabb[i] = std::min(aabb[i], aabb_tri[i]);
                    }
                    for (std::size_t i = 3; i < 6; ++i) {
                        aabb[i] = std::max(aabb[i], aabb_tri[i]);
                    }
                }
            }
        } else {
            m_left->AABB_iterator(aabb);
            m_right->AABB_iterator(aabb);
        }
    }

private:
    unsigned int m_D = 0;
    T m_plane = 0;
    std::vector<U> m_triangles;

    std::unique_ptr<KDTree<U>> m_left = nullptr;
    std::unique_ptr<KDTree<U>> m_right = nullptr;
    friend class KDTree<U>;
};
}