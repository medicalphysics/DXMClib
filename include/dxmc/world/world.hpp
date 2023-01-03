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

Copyright 2023 Erlend Andersen
*/

#pragma once

#include "dxmc/floating.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/worlditembase.hpp"

#include <concepts>

namespace dxmc {

template <typename U, typename T>
concept WorldItemType = std::derived_from<WorldItemBase<T>, U>;

template <typename U, typename... Us>
concept AnyWorldItemType = (... or std::same_as<U, Us>);

template <Floating T, WorldItemType<T>... Us>
class World2 {
public:
    World2()
    {
    }

private:
    std::tuple<std::vector<Us>...> m_items;
};
}
