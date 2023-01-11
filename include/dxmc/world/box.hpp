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
class Box final : public WorldItemBase<T> {
public:
    Box(const std::array<T, 6>& aabb = { -1, -1, -1, 1, 1, 1 })
        : WorldItemBase<T>()
        , m_aabb(aabb)
    {
    }

    void translate(const std::array<T, 3>& dist) override
    {
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] += dist[i];
            m_aabb[i + 3] += dist[i];
        }
    }
    std::array<T, 3> center() const override
    {
        std::array<T, 3> c {
            (m_aabb[0] + m_aabb[3]) * T { 0.5 },
            (m_aabb[1] + m_aabb[4]) * T { 0.5 },
            (m_aabb[2] + m_aabb[5]) * T { 0.5 },
        };
        return c;
    }

    std::array<T, 6> AABB() const override
    {
        return m_aabb;
    }
    std::optional<T> intersect(const Particle<T>& p) const override
    {
        const auto t = WorldItemBase<T>::intersectAABB<0>(p, m_aabb);
        if (t) {
            return (*t)[0] < T { 0 } ? std::make_optional((*t)[1]) : std::make_optional((*t)[0]);
        }
        return std::nullopt;
    }

    T transport(Particle<T>& p, RandomState& state) override
    {
        return 0;
    }

protected:
private:
    std::array<T, 6> m_aabb;
};

}
