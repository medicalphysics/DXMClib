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
#include "dxmc/world/visualizationintersectionresult.hpp"
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

        // Branched version with early exits
        template <Floating T>
        WorldIntersectionResult<T> intersect(const Particle<T>& p, const std::array<T, 6>& aabb)
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
        // Branched version with early exits
        template <Floating T>
        VisualizationIntersectionResult<T> intersectVisualization(const Particle<T>& p, const std::array<T, 6>& aabb)
        {
            VisualizationIntersectionResult<T> res;
            int axis_min = 0;
            int axis_max = 0;

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

            if (tymin > t[0]) {
                t[0] = tymin;
                axis_min = 1;
            }
            if (tymax < t[1]) {
                t[1] = tymax;
                axis_max = 1;
            }

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
            if (tzmax < t[1]) {
                t[1] = tzmax;
                axis_min = 2;
            }
            if (t[1] < 0) {
                return res;
                axis_max = 2;
            }

            res.rayOriginIsInsideItem = t[0] < 0;
            res.intersection = res.rayOriginIsInsideItem ? t[1] : t[0];
            res.intersectionValid = true;
            const int axis = res.rayOriginIsInsideItem ? axis_max : axis_min;
            res.normal[axis] = 1;
            return res;
        }
    }
}
}