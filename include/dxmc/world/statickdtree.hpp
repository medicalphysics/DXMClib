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

namespace dxmc {

template <typename U, typename T>
concept KDTreeType = requires(U u, Particle<T> p, std::array<T, 3> vec) {
                         Floating<T>;
                         //u <=> u;
                         { u.intersect(p) } -> std::same_as<std::optional<T>>;
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
class KDTree {
public:
    template <AnyKDTreeType<Us...> U>
    void insert(U item)
    {
        std::get<std::vector<U>>(m_data).push_back(item);
    }

private:
    std::tuple<std::vector<Us>...> m_data;
};

}