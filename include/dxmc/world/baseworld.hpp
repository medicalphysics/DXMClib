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
#include "dxmc/world/kdtree.hpp"
//#include "dxmc/material.hpp"
#include "dxmc/attenuationlut.hpp"

#include <array>
#include <optional>
#include <vector>

namespace dxmc {

template <Floating T>
struct WorldResult {
};

template <Floating T>
class BaseWorld {
public:
    // virtual void setMaxPhotonEnergy(const T maxEnergy) = 0;
    template <int FORWARD = 1>
    std::optional<T> intersect(const Particle<T>& p) const
    {
        std::array<T, 2> t {
            std::numeric_limits<T>::lowest(),
            std::numeric_limits<T>::max()
        };
        for (std::size_t i = 0; i < 3; i++) {
            if (std::abs(p.dir[i]) > std::numeric_limits<T>::epsilon()) {
                const auto d_inv = T { 1 } / p.dir[i];
                const auto t1 = (m_aabb[i] - p.pos[i]) * d_inv;
                const auto t2 = (m_aabb[i + 3] - p.pos[i]) * d_inv;
                const auto t_min_cand = std::min(t1, t2);
                const auto t_max_cand = std::max(t1, t2);
                t[0] = std::max(t[0], t_min_cand);
                t[1] = std::min(t[1], t_max_cand);
            }
        }
        if (t[0] > t[1])
            return std::nullopt;
        else {
            if constexpr (FORWARD == 1) {
                if (t[0] >= T { 0 } && t[1] >= T { 0 }) {
                    return std::min(t[0], t[1])
                } else {
                    const auto t_cand = std::max(t[0], t[1]);
                    return t_cand < T { 0 } ? std::nullopt : std::make_optional(t_cand);
                }
            } else if constexpr (FORWARD == -1) {
                if (t[0] <= T { 0 } && t[1] <= T { 0 }) {
                    return std::max(t[0], t[1])
                } else {
                    const auto t_cand = std::min(t[0], t[1]);
                    return t_cand > T { 0 } ? std::nullopt : std::make_optional(t_cand);
                }
            } else {
                return std::abs(t[0]) < std::abs(t[1]) ? std::make_optional(t[0]) : std::make_optional(t[1]);
            }
        }
    }

protected:
private:
    std::array<T, 6> m_aabb = { -1, -1, -1, 1, 1, 1 };
    // AttenuationLut<T> m_attenuationLut;
    // std::vector<Material<T>> m_materials;
};
}