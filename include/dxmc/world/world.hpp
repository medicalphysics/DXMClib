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
#include "dxmc/interactions.hpp"
#include "dxmc/material/material.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/basicshapes/aabb.hpp"
#include "dxmc/world/kdtree.hpp"
#include "dxmc/world/worlditems/worlditembase.hpp"

#include <concepts>
#include <tuple>
#include <vector>

namespace dxmc {

template <typename U, typename T>
concept WorldItemType = (std::derived_from<U, WorldItemBase<T>>);

template <typename U, typename... Us>
concept AnyWorldItemType = (... or std::same_as<U, Us>);

template <typename T, WorldItemType<T>... Us>
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

    void build(T AABB_padding = 10)
    {
        auto ptrs = getItemPointers();
        m_kdtree = KDTree(ptrs);
        m_aabb = m_kdtree.AABB();

        // adding padding
        const auto padding = std::max(AABB_padding, T { 0.1 }); // always at least 1 mm padding
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] = m_aabb[i] - padding;
            m_aabb[i + 3] = m_aabb[i + 3] + padding;
        }
    }

    const std::array<T, 6>& AABB() const
    {
        return m_aabb;
    }

    inline auto intersect(const Particle<T>& p)
    {
        return m_kdtree.intersectForward(p, m_aabb);
    }

    void transport(Particle<T>& p, RandomState& state)
    {
        bool continueSampling = true;
        bool updateAttenuation = true;

        {
            // transport particle to aabb
            const auto t = basicshape::AABB::intersect(p, m_aabb);
            if (t.valid()) {
                if (!t.rayOriginIsInsideItem) {
                    p.border_translate(t.intersection);
                }
            } else {
                continueSampling = false;
            }
        }

        T attenuationTotalInv;
        AttenuationValues<T> att;
        while (continueSampling) {
            if (updateAttenuation) {
                att = m_fillMaterial.attenuationValues(p.energy);
                attenuationTotalInv = T { 1 } / (att.sum() * m_fillMaterialDensity);
                updateAttenuation = false;
            }

            const auto r1 = state.randomUniform<T>();
            const auto stepLenght = -std::log(r1) * attenuationTotalInv; // cm

            // where do we hit an object
            const auto intersection = m_kdtree.intersect(p, m_aabb);
            const auto intersectionLenght = intersection.valid() ? intersection.intersection : std::numeric_limits<T>::max();

            if (intersectionLenght < stepLenght) {
                if (!intersection.rayOriginIsInsideItem) {
                    p.border_translate(intersectionLenght);
                }
                intersection.item->transport(p, state);
                continueSampling = p.energy > 0;
            } else {
                p.translate(stepLenght);
                if (basicshape::AABB::pointInside(p.pos, m_aabb)) {
                    const auto interactionResult = interactions::interact(att, p, m_fillMaterial, state);
                    updateAttenuation = interactionResult.particleEnergyChanged;
                    continueSampling = interactionResult.particleAlive;
                    m_dose.scoreEnergy(interactionResult.energyImparted);
                } else {
                    continueSampling = false;
                }
            }
        }
    }

private:
    std::array<T, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };
    std::tuple<std::vector<Us>...> m_items;
    KDTree<T> m_kdtree;
    Material2<T> m_fillMaterial;
    T m_fillMaterialDensity = T { 0.001225 };
    DoseScore<T> m_dose;
};
}
