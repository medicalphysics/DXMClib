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

#include <array>
#include <optional>

namespace dxmc {
namespace basicshape {
    namespace AABB {
        template <Floating F>
        static std::optional<std::array<F, 2>> intersect(const Particle<F>& p, const std::array<F, 6>& aabb)
        {
            std::array<F, 2> t;

            const auto dx = F { 1 } / p.dir[0];
            if (dx >= 0) {
                t[0] = (aabb[0] - p.pos[0]) * dx;
                t[1] = (aabb[3] - p.pos[0]) * dx;
            } else {
                t[0] = (aabb[3] - p.pos[0]) * dx;
                t[1] = (aabb[0] - p.pos[0]) * dx;
            }
            F tymin, tymax;
            const auto dy = F { 1 } / p.dir[1];
            if (dy >= 0) {
                tymin = (aabb[1] - p.pos[1]) * dy;
                tymax = (aabb[4] - p.pos[1]) * dy;
            } else {
                tymin = (aabb[4] - p.pos[1]) * dy;
                tymax = (aabb[1] - p.pos[1]) * dy;
            }
            if ((t[0] > tymax) || (tymin > t[1]))
                return std::nullopt;

            if (tymin > t[0])
                t[0] = tymin;
            if (tymax < t[1])
                t[1] = tymax;

            F tzmin, tzmax;
            const auto dz = F { 1 } / p.dir[2];
            if (dz >= 0) {
                tzmin = (aabb[2] - p.pos[2]) * dz;
                tzmax = (aabb[5] - p.pos[2]) * dz;
            } else {
                tzmin = (aabb[5] - p.pos[2]) * dz;
                tzmax = (aabb[2] - p.pos[2]) * dz;
            }

            if ((t[0] > tzmax) || (tzmin > t[1]))
                return std::nullopt;
            if (tzmin > t[0])
                t[0] = tzmin;
            if (tzmax < t[1])
                t[1] = tzmax;
            return std::make_optional(t);
        }

        template <Floating F>
        static bool pointInsideAABB(const std::array<F, 3>& p, const std::array<F, 6>& aabb)
        {
            return aabb[0] <= p[0] && p[0] <= aabb[3] && aabb[1] <= p[1] && p[1] <= aabb[4] && aabb[2] <= p[2] && p[2] <= aabb[5];
        }

        template <Floating F>
        static std::optional<F> intersectForward(const Particle<F>& p, const std::array<F, 6>& aabb)
        {
            const auto t_cand = intersect(p, aabb);
            if (t_cand) {
                const auto& t = t_cand.value();
                return pointInsideAABB(p.pos, aabb) ? std::make_optional(t[1]) : std::make_optional(t[0]);
            }
            return std::nullopt;
        }

    }
}
}