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
#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/worlditembase.hpp"

#include <limits>
#include <optional>

namespace dxmc {

template <Floating T>
class WorldSphere final : public WorldItemBase<T> {
public:
    WorldSphere(T radius = T { 16 }, const std::array<T, 3>& pos = { 0, 0, 0 })
        : WorldItemBase<T>()
        , m_radius(radius)
        , m_center(pos)
    {
    }

    void translate(const std::array<T, 3>& dist) override
    {
        m_center = vectormath::add(m_center, dist);
    }
    std::array<T, 3> center() const override
    {
        return m_center;
    }

    std::array<T, 6> AABB() const override
    {
        std::array<T, 6> aabb {
            m_center[0] - m_radius,
            m_center[1] - m_radius,
            m_center[2] - m_radius,
            m_center[0] + m_radius,
            m_center[1] + m_radius,
            m_center[2] + m_radius
        };
        return aabb;
    }
    std::optional<T> intersect(const Particle<T>& p) const override
    {
        const auto tbox = WorldItemBase<T>::intersectAABB(p, AABB());
        if (tbox)
            return intersectSphere(p, m_center, m_radius, *tbox);
        return std::nullopt;
    }

    T transport(Particle<T>& p, RandomState& state) override
    {
        return 0;
    }

protected:
private:
    T m_radius = 0;
    std::array<T, 3> m_center;
};

}
