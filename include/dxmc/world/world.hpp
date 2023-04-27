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

template <Floating T, WorldItemType<T>... Us>
class World2 {
public:
    World2()
        : m_fillMaterial(Material2<T>::byNistName("Air, Dry (near sea level)").value())
    {
    }

    void setMaterial(const Material2<T>& mat)
    {
        m_fillMaterial = mat;
    }
    void setMaterial(const Material2<T>& mat, T dens)
    {
        m_fillMaterial = mat;
        m_fillMaterialDensity = std::abs(dens);
    }
    void setMaterialDensity(T dens)
    {
        m_fillMaterialDensity = std::abs(dens);
    }

    template <AnyWorldItemType<Us...> U>
    auto& addItem(U& item)
    {
        std::get<std::vector<U>>(m_items).push_back(item);
        return std::get<std::vector<U>>(m_items).back();
    }

    template <AnyWorldItemType<Us...> U>
    auto& addItem(U&& item)
    {
        std::get<std::vector<U>>(m_items).push_back(std::move(item));
        return std::get<std::vector<U>>(m_items).back();
    }

    template <AnyWorldItemType<Us...> U>
    auto& addItem()
    {
        U item;
        std::get<std::vector<U>>(m_items).push_back(std::move(item));
        return std::get<std::vector<U>>(m_items).back();
    }

    template <AnyWorldItemType<Us...> U>
    const auto& getItems() const
    {
        return std::get<std::vector<U>>(m_items);
    }

    template <AnyWorldItemType<Us...> U>
    auto& getItems()
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

    void clearDose()
    {
        m_dose.clear();
        for (auto item : getItemPointers()) {
            item->clearDose();
        }
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

    std::array<T, 3> center() const
    {
        const auto [l, r] = vectormath::splice(m_aabb);
        return vectormath::scale(T { 0.5 }, vectormath::add(l, r));
    }

    inline auto intersect(const Particle<T>& p)
    {
        return m_kdtree.intersect(p, m_aabb);
    }

    inline auto intersectVisualization(const Particle<T>& p)
    {
        return m_kdtree.intersectVisualization(p, m_aabb);
    }

    inline bool transportParticleToWorld(Particle<T>& p)
    {
        if (!basicshape::AABB::pointInside(p.pos, m_aabb)) {
            const auto t = basicshape::AABB::intersect(p, m_aabb);
            if (t.valid()) {
                p.border_translate(t.intersection);
                return true;
            } else {
                return false;
            }
        }
        return true;
    }

    void transport(Particle<T>& p, RandomState& state)
    {
        bool continueSampling = transportParticleToWorld(p);
        bool updateAttenuation = true;

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
            const KDTreeIntersectionResult<T, WorldItemBase<T>> intersection = m_kdtree.intersect(p, m_aabb);

            if (intersection.valid()) { // Do we intersect anything?
                if (intersection.intersection < stepLenght) {
                    // We object is closer than free path.
                    if (!intersection.rayOriginIsInsideItem) { // if we are not already inside the object (we seldom are)
                        p.border_translate(intersection.intersection);
                    }
                    intersection.item->transport(p, state);
                    continueSampling = p.energy > 0;
                } else { // Free path is closer than object, we interact in the world empty space
                    p.translate(stepLenght);
                    const auto interactionResult = interactions::template interact<T, 5, 1>(att, p, m_fillMaterial, state);
                    updateAttenuation = interactionResult.particleEnergyChanged;
                    continueSampling = interactionResult.particleAlive;
                    m_dose.scoreEnergy(interactionResult.energyImparted);
                }
            } else { // We do not intersect any object
                p.translate(stepLenght);
                if (basicshape::AABB::pointInside(p.pos, m_aabb)) { // Are we still inside world?
                    const auto interactionResult = interactions::template interact<T, 5, 1>(att, p, m_fillMaterial, state);
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
