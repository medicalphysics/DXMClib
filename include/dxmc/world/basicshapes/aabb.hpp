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
        template <Floating T=double>
        inline constexpr bool pointInside(const std::array<T, 3>& p, const std::array<T, 6>& aabb)
        {
            return aabb[0] <= p[0] && p[0] <= aabb[3] && aabb[1] <= p[1] && p[1] <= aabb[4] && aabb[2] <= p[2] && p[2] <= aabb[5];
        }

        template <Floating T=double>
        inline constexpr bool overlap(const std::array<T, 6>& a, const std::array<T, 6>& b)
        {
            const bool x = a[0] <= b[3] && a[3] >= b[0];
            const bool y = a[1] <= b[4] && a[4] >= b[1];
            const bool z = a[2] <= b[5] && a[5] >= b[2];
            return x && y && z;
        }

        template <Floating T = double, bool FORWARD = true>
        std::optional<std::array<T, 2>> intersectForwardInterval(const Particle<T>& p, const std::array<T, 6>& aabb)
        {
            const std::array<T, 3> pdir_inv = { 1 / p.dir[0], 1 / p.dir[1], 1 / p.dir[2] };
            T t1 = (aabb[0] - p.pos[0]) * pdir_inv[0];
            T t2 = (aabb[3] - p.pos[0]) * pdir_inv[0];

            std::array<T, 2> tm = { std::min(t1, t2), std::max(t1, t2) };

            for (int i = 1; i < 3; ++i) {
                t1 = (aabb[i] - p.pos[i]) * pdir_inv[i];
                t2 = (aabb[i + 3] - p.pos[i]) * pdir_inv[i];

                tm[0] = std::max(tm[0], std::min(std::min(t1, t2), tm[1]));
                tm[1] = std::min(tm[1], std::max(std::max(t1, t2), tm[0]));
            }
            if constexpr (FORWARD)
                tm[0] = std::max(tm[0], T { 0 });
            return tm[1] > tm[0] ? std::make_optional(tm) : std::nullopt;
        }

        template <Floating T = double>
        WorldIntersectionResult<T> intersect(const Particle<T>& p, const std::array<T, 6>& aabb)
        {
            WorldIntersectionResult<T> res;
            if (const auto t_cand = intersectForwardInterval<T, true>(p, aabb); t_cand) {
                const auto& t = *t_cand;
                res.rayOriginIsInsideItem = t[0] <= 0;
                res.intersection = res.rayOriginIsInsideItem ? t[1] : t[0];
                res.intersectionValid = true;
            }
            return res;
        }

        template <Floating T = double, typename U>
        VisualizationIntersectionResult<T, U> intersectVisualization(const Particle<T>& p, const std::array<T, 6>& aabb)
        {
            VisualizationIntersectionResult<T, U> res;
            if (const auto t_cand = intersectForwardInterval<T, true>(p, aabb); t_cand) {
                const auto& t = *t_cand;
                res.rayOriginIsInsideItem = t[0] <= 0;
                res.intersection = res.rayOriginIsInsideItem ? t[1] : t[0];
                res.intersectionValid = true;
                std::array<T, 6> hit_merit;
                for (int i = 0; i < 3; ++i) {
                    const auto hit = p.pos[i] + res.intersection * p.dir[i];
                    hit_merit[i] = std::abs(aabb[i] - hit);
                    hit_merit[i + 3] = std::abs(aabb[i + 3] - hit);
                }
                const auto pos = std::distance(hit_merit.cbegin(), std::min_element(hit_merit.cbegin(), hit_merit.cend()));
                if (pos < 3)
                    res.normal[pos] = -1;
                else
                    res.normal[pos - 3] = 1;
            }
            return res;
        }
    }
}
}