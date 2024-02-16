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
#include <variant>
#include <vector>

namespace dxmc {

template <typename U>
concept WorldItemType = (std::derived_from<U, WorldItemBase>);

template <typename U, typename... Us>
concept AnyWorldItemType = (... or std::same_as<U, Us>);

template <WorldItemType... Us>
class World {
public:
    World()
        : m_fillMaterial(Material<double>::byNistName("Air, Dry (near sea level)").value())
    {
    }

    World(std::size_t reserveNumberOfWorldItems)
        : m_fillMaterial(Material<double>::byNistName("Air, Dry (near sea level)").value())
    {
        m_items.reserve(reserveNumberOfWorldItems);
    }

    void setMaterial(const Material<double>& mat)
    {
        m_fillMaterial = mat;
    }
    void setMaterial(const Material<double>& mat, double dens)
    {
        m_fillMaterial = mat;
        m_fillMaterialDensity = std::abs(dens);
    }
    void setMaterialDensity(double dens)
    {
        m_fillMaterialDensity = std::abs(dens);
    }

    void reserveNumberOfItems(std::size_t size)
    {
        m_items.reserve(size);
    }

    template <AnyWorldItemType<Us...> U>
    auto& addItem(U& item)
    {
        m_items.push_back(item);
        return std::get<U>(m_items.back());
    }

    template <AnyWorldItemType<Us...> U>
    auto& addItem(U&& item)
    {
        m_items.push_back(std::move(item));
        return std::get<U>(m_items.back());
    }

    template <AnyWorldItemType<Us...> U>
    auto& addItem()
    {
        U item;
        m_items.push_back(std::move(item));
        return std::get<U>(m_items.back());
    }

    const auto& getItems() const
    {
        return m_items;
    }

    auto& getItems()
    {
        return m_items;
    }

    std::vector<WorldItemBase*> getItemPointers()
    {
        std::vector<WorldItemBase*> ptrs(m_items.size());
        std::transform(m_items.begin(), m_items.end(), ptrs.begin(), [](auto& v) {
            return std::visit([](auto&& arg) -> WorldItemBase* { return &arg; }, v);
        });
        return ptrs;
    }

    std::vector<const WorldItemBase*> getItemPointers() const
    {
        std::vector<const WorldItemBase*> ptrs(m_items.size());
        std::transform(m_items.begin(), m_items.end(), ptrs.begin(), [](auto& v) {
            return std::visit([](auto&& arg) -> const WorldItemBase* { return &arg; }, v);
        });
        return ptrs;
    }

    void clearEnergyScored()
    {
        m_energyScored.clear();
        for (auto& v : m_items) {
            std::visit([](auto&& arg) { arg.clearEnergyScored(); }, v);
        }
    }

    void clearDoseScored()
    {
        for (auto& v : m_items) {
            std::visit([](auto&& arg) { arg.clearDoseScored(); }, v);
        }
    }

    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        for (auto& v : m_items) {
            std::visit([calibration_factor](auto&& arg) { arg.addEnergyScoredToDoseScore(calibration_factor); }, v);
        }
    }

    void build(double AABB_padding = 10)
    {
        auto ptrs = getItemPointers();
        m_kdtree = KDTree(ptrs);
        m_aabb = m_kdtree.AABB();

        // adding padding
        const auto padding = std::max(AABB_padding, 0.1); // always at least 1 mm padding
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] = m_aabb[i] - padding;
            m_aabb[i + 3] = m_aabb[i + 3] + padding;
        }
    }

    const std::array<double, 6>& AABB() const
    {
        return m_aabb;
    }

    std::array<double, 3> center() const
    {
        const auto [l, r] = vectormath::splice(m_aabb);
        return vectormath::scale(0.5, vectormath::add(l, r));
    }

    void translate(const std::array<double, 3> dist)
    {
        for (auto& v : m_items)
            std::visit([&dist](auto&& arg) { arg.translate(dist); }, v);

        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] += dist[i];
            m_aabb[i + 3] += dist[i];
        }
    }

    inline auto intersect(const Particle& p)
    {
        return m_kdtree.intersect(p, m_aabb);
    }

    inline auto intersectVisualization(const Particle& p)
    {
        return m_kdtree.intersectVisualization(p, m_aabb);
    }

    inline bool transportParticleToWorld(Particle& p) const
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

    void transport(Particle& p, RandomState& state)
    {
        bool continueSampling = transportParticleToWorld(p);
        bool updateAttenuation = true;

        double attenuationTotalInv;
        AttenuationValues<double> att;
        while (continueSampling) {
            if (updateAttenuation) {
                att = m_fillMaterial.attenuationValues(p.energy);
                attenuationTotalInv = 1 / (att.sum() * m_fillMaterialDensity);
                updateAttenuation = false;
            }

            const auto r1 = state.randomUniform();
            const auto stepLength = -std::log(r1) * attenuationTotalInv; // cm

            // where do we hit an object
            const auto intersection = m_kdtree.intersect(p, m_aabb);

            if (intersection.valid()) { // Do we intersect anything?
                if (intersection.intersection < stepLength) {
                    // Object is closer than free path.
                    if (!intersection.rayOriginIsInsideItem) { // if we are not already inside the object (we seldom are)
                        p.border_translate(intersection.intersection);
                    }
                    intersection.item->transport(p, state);
                    continueSampling = p.energy > 0;
                } else { // Free path is closer than object, we interact in the world empty space
                    p.translate(stepLength);
                    const auto interactionResult = interactions::template interact<5, 1>(att, p, m_fillMaterial, state);
                    updateAttenuation = interactionResult.particleEnergyChanged;
                    continueSampling = interactionResult.particleAlive;
                    m_energyScored.scoreEnergy(interactionResult.energyImparted);
                }
            } else { // We do not intersect any object
                p.translate(stepLength);
                if (basicshape::AABB::pointInside(p.pos, m_aabb)) { // Are we still inside world?
                    const auto interactionResult = interactions::template interact<5, 1>(att, p, m_fillMaterial, state);
                    updateAttenuation = interactionResult.particleEnergyChanged;
                    continueSampling = interactionResult.particleAlive;
                    m_energyScored.scoreEnergy(interactionResult.energyImparted);
                } else {
                    continueSampling = false;
                }
            }
        }
    }

private:
    std::array<double, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };
    std::vector<std::variant<Us...>> m_items;
    KDTree m_kdtree;
    Material<double, 5> m_fillMaterial;
    double m_fillMaterialDensity = 0.001225;
    EnergyScore m_energyScored;
};
}
