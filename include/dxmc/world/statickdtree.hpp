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

#include <array>
#include <concepts>
#include <memory>
#include <numeric>
#include <optional>
#include <tuple>
#include <variant>
#include <vector>

namespace dxmc {

template <typename U, typename T>
concept KDTreeType = requires(U u, Particle<T> p, std::array<T, 3> vec) {
                         Floating<T>;
                         u.translate(vec);
                         {
                             u.intersect(p)
                             } -> std::same_as<std::optional<T>>;
                         {
                             u.center()
                             } -> std::convertible_to<std::array<T, 3>>;
                         {
                             u.AABB()
                             } -> std::convertible_to<std::array<T, 6>>;
                     };
template <typename U, typename... Us>
concept AnyKDTreeType = (... or std::same_as<U, Us>);

template <Floating T, KDTreeType<T>... Us>
class KDTreeNode {
public:
    KDTreeNode() { }
    KDTreeNode(const std::vector<std::variant<Us...>>& data, std::size_t max_depth = 8)
    {
        buildTree(data, max_depth);
    }
    KDTreeNode(const std::vector<std::variant<Us...>>& data, const std::vector<std::size_t>& idx, std::size_t max_depth = 8)
    {
        buildTree(data, idx, max_depth);
    }

    void clear()
    {
        if (m_left)
            m_left->clear();
        if (m_right)
            m_right->clear();
        m_idx.clear();
        m_left = nullptr;
        m_right = nullptr;
    }
    std::size_t depth() const
    {
        if (m_left)
            return std::max(m_left->depth(), m_right->depth()) + 1;
        return 0;
    }
    void buildTree(const std::vector<std::variant<Us...>>& data, std::size_t max_depth = 8)
    {
        std::vector<std::size_t> idx(data.size());
        std::iota(idx.begin(), idx.end(), std::size_t { 0 });
        buildTree(data, idx, max_depth);
    }
    void buildTree(const std::vector<std::variant<Us...>>& data, const std::vector<std::size_t>& idx, std::size_t max_depth = 8)
    {
        if (data.size() == 0)
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
        for (const auto& v : data) {
            std::array<T, 6> aabb_obj = std::visit([&aabb](const auto& obj) -> std::array<T, 6> { return obj.AABB(); }, v);
            for (std::size_t i = 0; i < 3; ++i) {
                aabb[i] = std::min(aabb[i], aabb_obj[i]);
            }
            for (std::size_t i = 3; i < 6; ++i) {
                aabb[i] = std::max(aabb[i], aabb_obj[i]);
            }
        }

        const std::array<T, 3> extent { aabb[3] - aabb[0], aabb[4] - aabb[1], aabb[5] - aabb[2] };

        const auto D = vectormath::argmax3<std::uint_fast32_t, T>(extent);
        m_D = D;

        // finding split
        const auto median_split = planeSplit(data, idx, m_D);
        const auto fom = figureOfMerit(data, idx, m_D, median_split);

        if (fom == idx.size() || max_depth <= 1 || idx.size() <= 1) {
            m_idx = idx;
        } else {
            m_plane = median_split;
            std::vector<std::size_t> left, right;

            auto func = [=](const auto& o) -> int {
                return KDTreeNode<T, Us...>::planeSide(o, m_plane, m_D);
            };
            for (const auto i : idx) {
                const auto side = std::visit(func, data[i]);
                if (side <= 0)
                    left.push_back(i);
                if (side >= 0)
                    right.push_back(i);
            }
            m_left = std::make_unique<KDTreeNode<T, Us...>>(data, left, max_depth - 1);
            m_right = std::make_unique<KDTreeNode<T, Us...>>(data, right, max_depth - 1);
        }
    }

    std::optional<T> intersect(const std::vector<std::variant<Us...>>& data, const Particle<T>& particle, const std::array<T, 6>& aabb, const std::array<T, 2>& tbox) const
    {
        const auto& inter = intersectAABB(particle, aabb);
        if (inter) {
            const std::array<T, 2> tbox_min { std::min((*inter)[0], tbox[0]), std::min((*inter)[1], tbox[1]) };
            return intersect(data, particle, tbox_min);
        } else
            return std::nullopt;
    }
    std::optional<T> intersect(const std::vector<std::variant<Us...>>& data, const Particle<T>& particle, const std::array<T, 6>& aabb) const
    {
        const auto& inter = intersectAABB(particle, aabb);
        return inter ? intersect(data, particle, *inter) : std::nullopt;
    }
    std::optional<T> intersect(const std::vector<std::variant<Us...>>& data, const Particle<T>& particle, const std::array<T, 2>& tbox) const
    {
        if (!m_left) { // this is a leaf
            // intersect triangles between tbox and return;
            T t = std::numeric_limits<T>::max();
            auto func = [&particle](const auto& obj) -> std::optional<T> { return obj.intersect(particle); };
            for (const auto i : m_idx) {
                const auto t_cand = std::visit(func, data[i]);
                if (t_cand) {
                    t = std::min(t, *t_cand);
                }
            }
            return greaterOrEqual(t, tbox[0]) && lessOrEqual(t, tbox[1]) ? std::make_optional(t) : std::nullopt;
        }

        // test for parallell beam
        if (std::abs(particle.dir[m_D]) <= std::numeric_limits<T>::epsilon()) {
            const auto hit_left = m_left->intersect(data, particle, tbox);
            const auto hit_right = m_right->intersect(data, particle, tbox);
            if (hit_left && hit_right)
                return std::min(hit_left, hit_right);
            if (!hit_left)
                return hit_right;
            return hit_left;
        }

        const auto* const front = particle.dir[m_D] > T { 0 } ? m_left.get() : m_right.get();
        const auto* const back = particle.dir[m_D] > T { 0 } ? m_right.get() : m_left.get();

        const auto t = (m_plane - particle.pos[m_D]) / particle.dir[m_D];

        if (t <= tbox[0]) {
            // back only
            return back->intersect(data, particle, tbox);
        } else if (t >= tbox[1]) {
            // front only
            return front->intersect(data, particle, tbox);
        }

        // both directions (start with front)
        const std::array<T, 2> t_front { tbox[0], t };
        const auto hit = front->intersect(data, particle, t_front);
        if (hit) {
            if (*hit <= t) {
                return hit;
            }
        }
        const std::array<T, 2> t_back { t, tbox[1] };
        return back->intersect(data, particle, t_back);
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

    static T planeSplit(const std::vector<std::variant<Us...>>& data, const std::vector<std::size_t>& idx, std::uint_fast32_t D)
    {
        // mean split
        auto func = [&](const auto& obj) -> T { return obj.center()[D]; };

        const auto plane_mean = std::transform_reduce(idx.cbegin(), idx.cend(), T { 0 }, std::plus<>(), [&](const auto& i) -> T {
            return std::visit(func, data[i]);
        }) / idx.size();

        return plane_mean;

        // median split
        /*
        std::sort(idx.begin(), idx.end(), [&](auto lh, auto rh) {
            const auto lh_value = std::visit(func, data[lh]);
            const auto rh_value = std::visit(func, data[rh]);
            return lh_value < rh_value;
        });

        const auto N = idx.size();
        if (N % 2 == 1) {
            return std::visit(func, data[N / 2]);
        }
        return (std::visit(func, data[N / 2]) + std::visit(func, data[N / 2 - 1])) / 2;
        */
    }

    static int figureOfMerit(const std::vector<std::variant<Us...>>& objects, const std::vector<std::size_t>& idx, std::uint_fast32_t D, T planesep)
    {
        int fom = 0;
        int shared = 0;
        auto func = [planesep, D](const auto& o) -> int {
            return KDTreeNode<T, Us...>::planeSide(o, planesep, D);
        };
        for (const auto& i : idx) {
            const auto side = std::visit(func, objects[i]);
            fom += side;
            if (side == 0)
                shared++;
        }
        return std::abs(fom) + shared;
    }

    template <AnyKDTreeType<Us...> U>
    static int planeSide(const U& obj, const T plane, std::uint_fast32_t D)
    {
        T max = std::numeric_limits<T>::lowest();
        T min = std::numeric_limits<T>::max();
        const auto aabb = obj.AABB();
        min = aabb[D];
        max = aabb[D + 3];

        if (lessOrEqual(max, plane))
            return -1;
        if (greaterOrEqual(min, plane))
            return 1;
        return 0;
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
    std::vector<std::size_t> m_idx;
    std::unique_ptr<KDTreeNode<T, Us...>> m_left = nullptr;
    std::unique_ptr<KDTreeNode<T, Us...>> m_right = nullptr;
};

template <Floating T, KDTreeType<T>... Us>
class KDTree {
public:
    void build()
    {
        m_node.buildTree(m_data);
        m_aabb = AABB();
    }
    void clear()
    {
        m_node.clear();
    }
    std::size_t depth() const
    {
        return m_node.depth();
    }
    template <AnyKDTreeType<Us...> U>
    void insert(U item)
    {
        m_data.emplace_back(std::variant<Us...>(item));
        clear();
    }
    template <AnyKDTreeType<Us...> U>
    void insert(const std::vector<U>& items)
    {
        for (const auto& item : items) {
            m_data.emplace_back(std::variant<Us...>(item));
        }
        clear();
    }

    void translate(const std::array<T, 3>& vec)
    {
        for (auto& v : m_data) {
            std::visit(
                [&vec](auto&& obj) {
                    obj.translate(vec);
                },
                v);
        }
        clear();
    }

    std::array<T, 3> center() const
    {
        std::array<T, 3> cnt { 0, 0, 0 };
        std::size_t count = 0;

        for (auto& v : m_data) {
            auto c = std::visit([&cnt, &count](auto&& obj) -> std::array<T, 3> { return obj.center(); }, v);
            for (std::size_t i = 0; i < 3; ++i) {
                cnt[i] += c[i];
            }
        }

        if (count > 0) {
            for (std::size_t i = 0; i < cnt.size(); ++i)
                cnt[i] /= count;
        }
        return cnt;
    }
    std::array<T, 6> AABB() const
    {
        std::array<T, 6> aabb {
            std::numeric_limits<T>::max(),
            std::numeric_limits<T>::max(),
            std::numeric_limits<T>::max(),
            std::numeric_limits<T>::lowest(),
            std::numeric_limits<T>::lowest(),
            std::numeric_limits<T>::lowest()
        };

        for (auto& v : m_data) {
            auto c = std::visit([&aabb](const auto& obj) -> std::array<T, 6> {
                auto aabb_obj = obj.AABB();
                for (std::size_t i = 0; i < 3; ++i) {
                    aabb[i] = std::min(aabb[i], aabb_obj[i]);
                }
                for (std::size_t i = 3; i < 6; ++i) {
                    aabb[i] = std::max(aabb[i], aabb_obj[i]);
                }
                return aabb_obj;
            },
                v);
        }
        return aabb;
    }

    std::optional<T> intersect(const Particle<T>& particle) const
    {
        return m_node.intersect(m_data, particle, m_aabb);
    }

private:
    std::array<T, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };
    std::vector<std::variant<Us...>> m_data;
    KDTreeNode<T, Us...> m_node;
};

}