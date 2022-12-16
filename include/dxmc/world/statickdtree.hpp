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
    KDTreeNode(const std::vector<std::variant<Us...>>& data)
    {                
        buildTree(data);
    }
    KDTreeNode(const std::vector<std::variant<Us...>>& data, const std::vector<std::size_t>& idx)
    {
        buildTree(data, idx);
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

    void buildTree(const std::vector<std::variant<Us...>>& data)
    {
        std::vector<std::size_t> idx(data.size());
        std::iota(idx.begin(), idx.end(), std::size_t { 0 });
        buildTree(data, idx);
    }
    void buildTree(const std::vector<std::variant<Us...>>& data, std::vector<std::size_t>& idx)
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

        m_D = vectormath::argmax3<std::uint_fast32_t, T>(extent);

        // finding split
        std::sort(idx.begin(), idx.end(), [&](auto lh, auto rh) {
            auto func = [&](const auto& obj) -> T { return obj.center()[m_D]; };
            const auto lh_value = std::visit(func, data[lh]);
            const auto rh_value = std::visit(func, data[rh]);
            return lh_value < rh_value;
        });
        m_plane = 0;
    }

protected:
    template <AnyKDTreeType<Us...> U>
    static int figureOfMerit(const std::tuple<std::vector<Us*>...>& data, const T planesep, std::atomic_uint_fast32_t D)
    {
    }

    template <AnyKDTreeType<Us...> U>
    static int planeSide(const U obj, const T plane, std::uint_fast32_t D)
    {
        T max = std::numeric_limits<T>::lowest();
        T min = std::numeric_limits<T>::max();
        const auto aabb = obj->AABB();
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
    }

    template <AnyKDTreeType<Us...> U>
    void insert(U item)
    {
        m_data.emplace_back(std::variant<Us...>(item));
    }
    template <AnyKDTreeType<Us...> U>
    void insert(const std::vector<U>& items)
    {
        for (const auto& item : items) {
            m_data.emplace_back(std::variant<Us...>(item));
        }
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

private:
    std::vector<std::variant<Us...>> m_data;
    KDTreeNode<T, Us...> m_node;
};

}