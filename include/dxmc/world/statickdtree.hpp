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

template <typename T, typename U>
concept StaticKDTreeType = requires(U u, Particle<T> p, std::array<T, 3> vec) {
                              // Floating<T>;
                               //u <=> u;
                               u.translate(vec);
                               /* {
                                   u.intersect(p)
                                   } -> std::same_as<std::optional<T>>;

                               {
                                   u.center()
                                   } -> std::same_as<std::array< T, 3>>;

                               {
                                   u.AABB()
                                   } -> std::same_as<std::array< T, 6>>;
                           */ };

template <typename T, typename... Ts>
concept same_as_any = (... or std::same_as<T, Ts>);


template <StaticKDTreeType U>
struct Test {
    Test(float n, U o):num(n),obj(o) {
    }
    float num;
    U obj;
};


/*
template <Floating T, StaticKDTreeType<T>... Us>
//template <Floating T, typename... Us>
class StaticKDTree {
public:
    StaticKDTree(std::int_fast32_t max_depth = 8)
        : m_maxDepth(max_depth)
    {
    }
    //template <typename U, typename... Us>
    template <typename U>
        requires(std::same_as<U, Us>)
    void insert(U value)
    {
        std::get<std::vector<U>>(m_objects).push_back(value);
    }

protected:
private:
    std::int_fast32_t m_maxDepth = 0;
    std::tuple<std::vector<Us>...> m_objects;
};
*/

}