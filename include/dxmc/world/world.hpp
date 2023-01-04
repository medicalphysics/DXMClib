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
#include "dxmc/world/kdtree.hpp"
#include "dxmc/world/worlditembase.hpp"

#include <concepts>
#include <tuple>
#include <vector>

namespace dxmc {

template <typename U, typename T>
concept WorldItemType = (std::derived_from<U, WorldItemBase<T>>);

template <typename U, typename... Us>
concept AnyWorldItemType = (... or std::same_as<U, Us>);

template <Floating T, WorldItemType<T>... Us>
class World2 {
public:
    World2()
    {
    }

    template <AnyWorldItemType<Us...> U>
    void addItem(U item)
    {
        std::get<std::vector<U>>(m_items).push_back(item);
    }

    template <AnyWorldItemType<Us...> U>
    const auto& getItems() const
    {
        return std::get<std::vector<U>>(m_items);
    }

    void build()
    {
        auto ptrs = itemPointers();
        m_kdtree = KDTree(ptrs);
        m_aabb = m_kdtree.AABB();
    }

    auto intersect(const Particle<T>& p) const
    {
        return m_kdtree.intersect(p, m_aabb);
    }

    // protected:
    std::vector<const WorldItemBase<T>*> itemPointers() const
    {
        std::vector<const WorldItemBase<T>*> ptrs;

        auto iter = [&ptrs](const auto& v) {
            for (const auto& item : v)
                ptrs.push_back(&item);
        };

        std::apply([&iter](const auto&... vec) {
            (iter(vec), ...);
        },
            m_items);

        return ptrs;
    }

private:
    std::array<T, 6> m_aabb;
    std::tuple<std::vector<Us>...> m_items;
    KDTree<T> m_kdtree;
};
}
