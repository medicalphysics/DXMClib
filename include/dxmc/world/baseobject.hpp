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
#include "dxmc/world/kdtree.hpp"

#include "dxmc/attenuationlut.hpp"

#include <array>
#include <optional>
#include <vector>

namespace dxmc {

template <Floating T>
struct ResultObject {
};

template <Floating T>
class BaseObject {
public:
    using Type = T;

    template <int FORWARD = 1>
    std::optional<T> intersectAABB(const Particle<T>& p) const
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
                    return std::min(t[0], t[1]);
                } else {
                    const auto t_cand = std::max(t[0], t[1]);
                    return t_cand < T { 0 } ? std::nullopt : t_cand;
                }
            } else if constexpr (FORWARD == -1) {
                if (t[0] <= T { 0 } && t[1] <= T { 0 }) {
                    return std::max(t[0], t[1]);
                } else {
                    const auto t_cand = std::min(t[0], t[1]);
                    return t_cand > T { 0 } ? std::nullopt : t_cand;
                }
            } else {
                return std::abs(t[0]) < std::abs(t[1]) ? t[0] : t[1];
            }
        }
    }
    virtual void translate(const std::array<T, 3>& dist) = 0;
    virtual std::array<T, 3> center() const = 0;
    const std::array<T, 6>& AABB() const { return m_aabb; }
    std::array<T, 6> AABB() { return m_aabb; }
    virtual std::optional<T> intersect(const Particle<T>& p) const = 0;
    virtual std::optional<T> intersect(const Particle<T>& p, const std::array<T, 2>& tbox) const = 0;
    virtual T transport(Particle<T>& p, RandomState& state) = 0;

protected:
    // virtual T stepLenght(const Particle<T>& p, RandomState& state) const = 0;
    std::array<T, 6> m_aabb { -1, -1, -1, 1, 1, 1 };
    ResultObject<T> m_result;

private:
};
}