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
#include "dxmc/world/worldintersectionresult.hpp"

#include <algorithm>
#include <array>
#include <optional>

namespace dxmc {
namespace basicshape {
    namespace AABB {

        template <Floating T>
        bool pointInside(const std::array<T, 3>& p, const std::array<T, 6>& aabb)
        {
            return aabb[0] <= p[0] && p[0] <= aabb[3] && aabb[1] <= p[1] && p[1] <= aabb[4] && aabb[2] <= p[2] && p[2] <= aabb[5];
        }

        template <Floating T>
        std::optional<std::array<T, 2>> intersectForwardInterval(const Particle<T>& p, const std::array<T, 6>& aabb)
        {
            std::array t = { T { 0 }, std::numeric_limits<T>::max() };

            for (std::size_t i = 0; i < 3; ++i) {
                const auto d = T { 1 } / p.dir[i];
                const auto t0 = (aabb[i] - p.pos[i]) * d;
                const auto t1 = (aabb[i + 3] - p.pos[i]) * d;
                if (d > T { 0 }) {
                    t[0] = std::max(t0, t[0]);
                    t[1] = std::min(t1, t[1]);
                } else {
                    t[0] = std::max(t1, t[0]);
                    t[1] = std::min(t0, t[1]);
                }
            }
            return t[0] < t[1] && t[1] > T { 0 } ? std::make_optional(t) : std::nullopt;
        }
        template <Floating T>
        std::optional<std::array<T, 2>> intersectForwardInterval2(const Particle<T>& p, const std::array<T, 6>& aabb)
        {
            std::array<T, 2> t;

            const auto dx = 1 / p.dir[0];
            if (dx >= 0) {
                t[0] = (aabb[0] - p.pos[0]) * dx;
                t[1] = (aabb[3] - p.pos[0]) * dx;
            } else {
                t[0] = (aabb[3] - p.pos[0]) * dx;
                t[1] = (aabb[0] - p.pos[0]) * dx;
            }
            T tymin, tymax;
            const auto dy = 1 / p.dir[1];
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

            T tzmin, tzmax;
            const auto dz = 1 / p.dir[2];
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
            if (t[1] < 0)
                return std::nullopt;
            if (t[0] < T { 0 })
                t[0] = T { 0 };

            return std::make_optional(t);
        }

        template <Floating T>
        WorldIntersectionResult<T> intersect(const Particle<T>& p, const std::array<T, 6>& aabb)
        {
            T tmin = std::numeric_limits<T>::lowest();
            T tmax = std::numeric_limits<T>::max();
            for (std::size_t i = 0; i < 3; ++i) {
                if (p.dir[i] < T { 0 } || p.dir[i] > T { 0 }) {
                    const auto d = 1 / p.dir[i];
                    const auto t0 = (aabb[i] - p.pos[i]) * d;
                    const auto t1 = (aabb[i + 3] - p.pos[i]) * d;
                    if (d > T { 0 }) {
                        tmin = std::max(t0, tmin);
                        tmax = std::min(t1, tmax);
                    } else {
                        tmin = std::max(t1, tmin);
                        tmax = std::min(t0, tmax);
                    }
                }
            }
            const auto particleInside = tmin < T { 0 } && tmax > T { 0 };

            WorldIntersectionResult<T> res = {
                .intersection = particleInside ? tmax : tmin,
                .rayOriginIsInsideItem = particleInside,
                .intersectionValid = tmax > 0 && tmax > tmin
            };
            return res;
        }

        // Branched version with early exits
        template <Floating T>
        WorldIntersectionResult<T> intersect2(const Particle<T>& p, const std::array<T, 6>& aabb)
        {
            WorldIntersectionResult<T> res;
            std::array<T, 2> t;

            const auto dx = 1 / p.dir[0];
            if (dx >= 0) {
                t[0] = (aabb[0] - p.pos[0]) * dx;
                t[1] = (aabb[3] - p.pos[0]) * dx;
            } else {
                t[0] = (aabb[3] - p.pos[0]) * dx;
                t[1] = (aabb[0] - p.pos[0]) * dx;
            }
            T tymin, tymax;
            const auto dy = 1 / p.dir[1];
            if (dy >= 0) {
                tymin = (aabb[1] - p.pos[1]) * dy;
                tymax = (aabb[4] - p.pos[1]) * dy;
            } else {
                tymin = (aabb[4] - p.pos[1]) * dy;
                tymax = (aabb[1] - p.pos[1]) * dy;
            }
            if ((t[0] > tymax) || (tymin > t[1]))
                return res;

            if (tymin > t[0])
                t[0] = tymin;
            if (tymax < t[1])
                t[1] = tymax;

            T tzmin, tzmax;
            const auto dz = 1 / p.dir[2];
            if (dz >= 0) {
                tzmin = (aabb[2] - p.pos[2]) * dz;
                tzmax = (aabb[5] - p.pos[2]) * dz;
            } else {
                tzmin = (aabb[5] - p.pos[2]) * dz;
                tzmax = (aabb[2] - p.pos[2]) * dz;
            }

            if ((t[0] > tzmax) || (tzmin > t[1]))
                return res;
            if (tzmin > t[0])
                t[0] = tzmin;
            if (tzmax < t[1])
                t[1] = tzmax;
            if (t[1] < 0)
                return res;

            res.rayOriginIsInsideItem = t[0] < 0;
            res.intersection = res.rayOriginIsInsideItem ? t[1] : t[0];
            res.intersectionValid = true;

            return res;
        }

    }
}
}