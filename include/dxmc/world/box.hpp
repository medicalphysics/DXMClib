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

#include "dxmc/dxmcrandom.hpp"
#include "dxmc/floating.hpp"
#include "dxmc/interactions.hpp"
#include "dxmc/material/material.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/worlditembase.hpp"
#include "dxmc/world/intersection.hpp"

#include <limits>
#include <optional>

namespace dxmc {

template <Floating T>
class Box final : public WorldItemBase<T> {
public:
    Box(const std::array<T, 6>& aabb = { -1, -1, -1, 1, 1, 1 })
        : WorldItemBase<T>()
        , m_aabb(aabb)
        , m_material(Material2<T>::byNistName("Air, Dry (near sea level)").value())
    {
        m_materialDensity = NISTMaterials<T>::density("Air, Dry (near sea level)");
    }
    Box(T aabb_size, std::array<T, 3> pos = { 0, 0, 0 })
        : WorldItemBase<T>()
        , m_material(Material2<T>::byNistName("Air, Dry (near sea level)").value())
    {
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] = -std::abs(aabb_size) + pos[i];
            m_aabb[i + 3] = std::abs(aabb_size) + pos[i];
        }
        m_materialDensity = NISTMaterials<T>::density("Air, Dry (near sea level)");
    }

    void translate(const std::array<T, 3>& dist) noexcept override
    {
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] += dist[i];
            m_aabb[i + 3] += dist[i];
        }
    }
    std::array<T, 3> center() const noexcept override
    {
        std::array<T, 3> c {
            (m_aabb[0] + m_aabb[3]) * T { 0.5 },
            (m_aabb[1] + m_aabb[4]) * T { 0.5 },
            (m_aabb[2] + m_aabb[5]) * T { 0.5 },
        };
        return c;
    }

    std::array<T, 6> AABB() const noexcept override
    {
        return m_aabb;
    }
    std::optional<T> intersect(const Particle<T>& p) const noexcept override
    {
        const auto t = intersectAABB(p, m_aabb);
        if (t) {
            return (*t)[0] < T { 0 } ? std::make_optional((*t)[1]) : std::make_optional((*t)[0]);
        }
        return std::nullopt;
    }

    void transport(Particle<T>& p, RandomState& state) noexcept override
    {
        bool cont = pointInsideAABB(p.pos, m_aabb);
        bool updateAtt = false;
        auto att = m_material.attenuationValues(p.energy);
        auto attSumInv = 1 / (att.sum() * m_materialDensity);
        while (cont) {
            if (updateAtt) {
                att = m_material.attenuationValues(p.energy);
                attSumInv = 1 / (att.sum() * m_materialDensity);
            }
            const auto stepLen = -std::log(state.randomUniform<T>()) * attSumInv; // cm
            const auto intLen = intersect(p).value(); // this must be not nullopt

            if (stepLen < intLen) {
                // interaction happends
                p.translate(stepLen);
                const auto r1 = state.randomUniform(att.sum());
                if (r1 < att.photoelectric) {
                    const auto E = interactions::photoelectricEffect(att.photoelectric, p, m_material, state);
                    m_dose.scoreEnergy(E);
                    if (p.energy == T { 0 }) {
                        cont = false;
                    }
                } else if (r1 < att.photoelectric + att.incoherent) {
                    const auto E = interactions::comptonScatter(p, m_material, state);
                    m_dose.scoreEnergy(E);
                    if (p.energy == T { 0 }) {
                        cont = false;
                    }
                } else {
                    interactions::rayleightScatter(p, m_material, state);
                }
            } else {
                // transport to border
                p.translate(std::nextafter(intLen, std::numeric_limits<T>::max()));
                cont = false;
            }
        }
    }

    const DoseScore<T>& dose(std::size_t index = 0) const override
    {
        return m_dose;
    }

protected:
    static bool pointInsideAABB(const std::array<T, 3>& p, const std::array<T, 6>& aabb)
    {
        bool inside = aabb[0] <= p[0] && p[0] <= aabb[3];
        inside = inside && aabb[1] <= p[1] && p[1] <= aabb[4];
        inside = inside && aabb[2] <= p[2] && p[2] <= aabb[5];
        return inside;
    }

private:
    std::array<T, 6> m_aabb;
    Material2<T> m_material;
    T m_materialDensity = 1;
    DoseScore<T> m_dose;
};

}
