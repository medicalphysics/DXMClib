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
#include <optional>
#include <tuple>
#include <vector>
#include <memory>

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
    KDTreeNode(const std::tuple<std::vector<Us*>...>& data, std::array<std::vector<std::uint_fast32_t>, std::tuple_size_v<std::tuple<std::vector<Us>...>>> idx)
    {
    }

    void clear() {
        if (m_left)
            m_left->clear();
        if (m_right)
            m_right->clear();
        m_data.clear();
        m_left = nullptr;
        m_right = nullptr;
    }

protected:
    static int figureOfMerit(const std::tuple<std::vector<Us*>...>& data, const T planesep, std::atomic_uint_fast32_t D)
    {
        int fom = 0;
        int shared = 0;

        auto func = [&](const auto& objects) -> void {
            for (auto obj : objects) {
                const auto side = KDTreeNode<T, Us...>::planeSide(obj, planesep, D);
                fom += side;
                if (side == 0) {
                    shared++;
                }
            }
        };
        std::apply([&func](auto&... objects) {
            (func(objects), ...);
        },
            data);

        return std::abs(fom) + shared;
    }

    template <AnyKDTreeType<Us...> U>
    static int planeSide(const U* obj, const T plane, std::uint_fast32_t D)
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
    std::tuple<std::vector<Us>...> m_data {};

    std::unique_ptr<KDTreeNode<T, Us...>> m_left = nullptr;
    std::unique_ptr<KDTreeNode<T, Us...>> m_right = nullptr;    
};

template <Floating T, KDTreeType<T>... Us>
class KDTree {
public:
    void build()
    {
        std::tuple<std::vector<Us>...> data_ptr;
        /*
        auto func = [&data_ptr](auto& objects) -> void {
            for (auto& obj : objects)

        };
        std::apply([this, &trans_func](auto&... objects) {
            (trans_func(objects), ...);
        },
            m_data);
            */

        KDTreeNode<T, Us...> m_node;
    }

    template <AnyKDTreeType<Us...> U>
    void insert(U item)
    {
        std::get<std::vector<U>>(m_data).push_back(item);
    }
    void translate(const std::array<T, 3>& vec)
    {
        auto trans_func = [&](auto& objects) -> void {
            for (auto& obj : objects)
                obj.translate(vec);
        };
        std::apply([this, &trans_func](auto&... objects) {
            (trans_func(objects), ...);
        },
            m_data);
    }
    std::array<T, 3> center() const
    {
        std::array<T, 3> cnt { 0, 0, 0 };
        std::size_t count = 0;
        auto center_func = [&](auto& objects) -> void {
            for (auto& obj : objects) {
                auto obj_cnt = obj.center();
                count++;
                for (std::size_t i = 0; i < cnt.size(); ++i)
                    cnt[i] += obj_cnt[i];
            }
        };
        std::apply([this, &center_func](auto&... objects) {
            (center_func(objects), ...);
        },
            m_data);

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

        auto aabb_func = [&aabb](auto& objects) -> void {
            for (auto& obj : objects) {
                auto aabb_obj = obj.AABB();
                for (std::size_t i = 0; i < 3; ++i)
                    aabb[i] = std::min(aabb[i], aabb_obj[i]);
                for (std::size_t i = 3; i < 6; ++i)
                    aabb[i] = std::max(aabb[i], aabb_obj[i]);
            }
        };
        std::apply([this, &aabb_func](auto&... objects) {
            (aabb_func(objects), ...);
        },
            m_data);
        return aabb;
    }

private:
    std::tuple<std::vector<Us>...> m_data {};
    KDTreeNode<T, Us...> m_node;
};

}