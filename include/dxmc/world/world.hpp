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

#include "dxmc/dxmcrandom.hpp"
#include "dxmc/floating.hpp"
#include "dxmc/material/material.hpp"
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
        : m_fillMaterial(Material2<T>::byNistName("Air, Dry (near sea level)").value())
    {
    }

    template <AnyWorldItemType<Us...> U>
    void addItem(U& item)
    {
        std::get<std::vector<U>>(m_items).push_back(item);
    }

    template <AnyWorldItemType<Us...> U>
    void addItem(U&& item)
    {
        std::get<std::vector<U>>(m_items).push_back(std::move(item));
    }

    template <AnyWorldItemType<Us...> U>
    const auto& getItems() const
    {
        return std::get<std::vector<U>>(m_items);
    }

    void build()
    {
        auto ptrs = getItemPointers();
        m_kdtree = KDTree(ptrs);
        m_aabb = m_kdtree.AABB();
    }

    auto intersect(const Particle<T>& p)
    {
        return m_kdtree.intersect(p, m_aabb);
    }
    void transport(Particle<T>& p, RandomState& state)
    {

        bool continueSampling = true;
        bool updateAttenuation = false;

        if (!pointInsideAABB(p.pos, m_aabb)) {
            // transport particle to aabb
            const auto t = WorldItemBase<T>::intersectAABB(p, m_aabb);
            if (t) {
                p.translate(std::nextafter(t.value()[0], std::numeric_limits<T>::max()));
            } else {
                continueSampling = false;
            }
        }

        T attenuationTotalInv = T { 1 } / (m_fillMaterial.attenuationValues(p.energy).sum() * m_fillMaterialDensity);

        while (continueSampling) {
            if (updateAttenuation) {
                attenuationTotalInv = T { 1 } / (m_fillMaterial.attenuationValues(p.energy).sum() * m_fillMaterialDensity);
            }

            const auto r1 = state.randomUniform<T>();
            const auto stepLenght = -std::log(r1) * attenuationTotalInv; // cm

            // where do we hit an object
            auto intersection = m_kdtree.intersect(p, m_aabb);
            const auto intersectionLenght = intersection.valid() ? intersection.intersection : std::numeric_limits<T>::max();

            if (intersectionLenght < stepLenght) {
                p.translate(std::nextafter(intersectionLenght, std::numeric_limits<T>::max()));
                intersection.item->transport(p, state);
            } else {
                p.translate(std::nextafter(stepLenght, std::numeric_limits<T>::max()));
                if (pointInsideAABB(p.pos, m_aabb)) {
                    const auto att = m_fillMaterial.attenuationValues(p.energy);
                    const auto r2 = state.randomUniform<T>(att.sum());
                    if (r2 < att.photoelectric) {

                    } else if (r2 < (att.photoelectric + att.incoherent)) {

                    } else {
                    }

                } else {
                    continueSampling = false;
                }
            }
        }
    }

protected:
    static bool pointInsideAABB(const std::array<T, 3>& p, const std::array<T, 6>& aabb)
    {
        bool inside = aabb[0] <= p[0] && p[0] <= aabb[3];
        inside = inside && aabb[1] <= p[1] && p[1] <= aabb[4];
        inside = inside && aabb[2] <= p[2] && p[2] <= aabb[5];
        return inside;
    }

    std::vector<WorldItemBase<T>*> getItemPointers()
    {
        std::vector<WorldItemBase<T>*> ptrs;

        auto iter = [&ptrs](auto& v) {
            for (auto& item : v)
                ptrs.push_back(&item);
        };

        std::apply([&iter](auto&... vec) {
            (iter(vec), ...);
        },
            m_items);

        return ptrs;
    }

private:
    std::array<T, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };
    std::tuple<std::vector<Us>...> m_items;
    KDTree<T> m_kdtree;
    Material2<T> m_fillMaterial;
    T m_fillMaterialDensity = T { 0.001225 };
    ResultObject<T> m_dose;
};
}
