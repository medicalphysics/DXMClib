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
    namespace cylinder {

        template <Floating T>
        bool pointInside(const std::array<T, 3>& pos, const std::array<T, 3>& center, const T radii, const T half_height)
        {
            const std::array<T, 2> dp = { pos[0] - center[0], pos[1] - center[1] };
            return (center[2] - half_height < pos[2]) && (pos[2] < center[2] + half_height) && ((dp[0] * dp[0] + dp[1] * dp[1]) < radii * radii);
        }

        template <Floating T>
        std::optional<std::array<T, 2>> intersectCylinderWall(const Particle<T>& p, const std::array<T, 3>& center, const T radii)
        {
            // nummeric stable ray sphere intersection in 2D
            const auto a = p.dir[0] * p.dir[0] + p.dir[1] * p.dir[1];
            if (a < std::numeric_limits<T>::epsilon())
                return std::nullopt;

            const auto r2 = radii * radii;
            const std::array f = { p.pos[0] - center[0], p.pos[1] - center[1] };

            // positive b mean center of sphere is in front of ray
            const auto b = -f[0] * p.dir[0] - f[1] * p.dir[1];

            // positive c means ray starts outside of sphere
            const auto c = f[0] * f[0] + f[1] * f[1] - r2;

            if ((c > 0) && (b < 0)) {
                // if ray starts outside sphere and center is begind ray
                // we exit early
                return std::nullopt;
            }

            const std::array delta1 = { f[0] + b * p.dir[0] / a, f[1] + b * p.dir[1] / a };

            const auto delta = r2 - (delta1[0] * delta1[0] + delta1[1] * delta1[1]);

            if (delta < 0) {
                // no solution to the quadratic equation (we miss)
                return std::nullopt;
            }

            const int sign = b > 0 ? 1 : -1;
            const auto q = b + sign * std::sqrt(a * delta);
            if (c < 0 && b < 0) {
                std::array t = { q / a, c / q };
                return std::make_optional(t);
            } else {
                std::array t = { c / q, q / a };
                return std::make_optional(t);
            }
        }
        template <Floating T>
        std::optional<T> intersectCylinderWallForward(const Particle<T>& p, const std::array<T, 3>& center, const T radii)
        {
            // nummeric stable ray sphere intersection in 2D
            const auto a = p.dir[0] * p.dir[0] + p.dir[1] * p.dir[1];
            if (a < std::numeric_limits<T>::epsilon())
                return std::nullopt;

            const auto r2 = radii * radii;
            const std::array f = { p.pos[0] - center[0], p.pos[1] - center[1] };

            // positive b mean center of sphere is in front of ray
            const auto b = -f[0] * p.dir[0] - f[1] * p.dir[1];

            // positive c means ray starts outside of sphere
            const auto c = f[0] * f[0] + f[1] * f[1] - r2;

            if ((c > 0) && (b < 0)) {
                // if ray starts outside sphere and center is begind ray
                // we exit early
                return std::nullopt;
            }

            const std::array delta1 = { f[0] + b * p.dir[0] / a, f[1] + b * p.dir[1] / a };

            const auto delta = r2 - (delta1[0] * delta1[0] + delta1[1] * delta1[1]);

            if (delta < 0) {
                // no solution to the quadratic equation (we miss)
                return std::nullopt;
            }

            const int sign = b > 0 ? 1 : -1;
            const auto q = b + sign * std::sqrt(a * delta);

            if (c < 0 && b > 0) {
                // inside sphere
                return std::make_optional(q / a);
            } else {
                return std::make_optional(c / q);
            }
        }
        template <Floating T>
        std::optional<T> intersectDiscZ(const Particle<T>& p, const std::array<T, 3>& center, const T radii)
        {
            if (std::abs(p.dir[2]) <= std::numeric_limits<T>::epsilon())
                return std::nullopt;
            const auto tz = (center[2] - p.pos[2]) / p.dir[2];

            const auto xz = p.pos[0] + p.dir[0] * tz - center[0];
            const auto yz = p.pos[1] + p.dir[1] * tz - center[1];
            if (xz * xz + yz * yz < radii * radii)
                return std::make_optional(tz);

            return std::nullopt;
        }

        template <Floating T>
        std::optional<T> intersectDiscZForward(const Particle<T>& p, const std::array<T, 3>& center, const T radii)
        {
            if (std::abs(p.dir[2]) <= std::numeric_limits<T>::epsilon())
                return std::nullopt;
            const auto tz = (center[2] - p.pos[2]) / p.dir[2];
            if (tz > 0) {
                const auto xz = p.pos[0] + p.dir[0] * tz - center[0];
                const auto yz = p.pos[1] + p.dir[1] * tz - center[1];
                if (xz * xz + yz * yz < radii * radii)
                    return std::make_optional(tz);
            }
            return std::nullopt;
        }

        template <Floating T>
        std::optional<std::array<T, 2>> intersect(const Particle<T>& p, const std::array<T, 3>& center, const T radii, const T half_height)
        {
            auto t_cand = intersectCylinderWall(p, center, radii);
            std::optional<T> tmin_cand, tmax_cand;
            if (t_cand) {
                const auto t = *t_cand;
                const auto pzmin = p.pos[2] + p.dir[2] * t[0];
                if ((center[2] - half_height) <= pzmin && pzmin <= (center[2] + half_height)) {
                    tmin_cand = std::make_optional(t[0]);
                }
                const auto pzmax = p.pos[2] + p.dir[2] * t[1];
                if ((center[2] - half_height) <= pzmax && pzmax <= (center[2] + half_height)) {
                    tmax_cand = std::make_optional(t[1]);
                }
            }
            std::array centerDisc = { center[0], center[1], center[2] - half_height };
            auto tdl_cand = intersectDiscZ(p, centerDisc, radii);
            centerDisc[2] = center[2] + half_height;
            auto tdh_cand = intersectDiscZ(p, centerDisc, radii);

            if (tmin_cand || tmax_cand || tdl_cand || tdh_cand) {
                auto min = std::min({ tmin_cand, tdl_cand, tdh_cand }, [](const auto& lh, const auto& rh) -> bool { return lh.value_or(std::numeric_limits<T>::max()) < rh.value_or(std::numeric_limits<T>::max()); });
                auto max = std::max({ tmax_cand, tdl_cand, tdh_cand }, [](const auto& lh, const auto& rh) -> bool { return lh.value_or(std::numeric_limits<T>::min()) < rh.value_or(std::numeric_limits<T>::min()); });

                return std::make_optional<std::array<T, 2>>({ *min, *max });
            }
            return std::nullopt;
        }

        template <Floating T>
        std::optional<T> intersectForward(const Particle<T>& p, const std::array<T, 3>& center, const T radii, const T half_height)
        {
            auto t_cand = intersectCylinderWallForward(p, center, radii);
            if (t_cand) {
                const auto pz = p.pos[2] + p.dir[2] * *t_cand;
                if ((pz < center[2] - half_height) || pz > (center[2] + half_height)) {
                    t_cand = std::nullopt;
                }
            }
            std::array centerDisc = { center[0], center[1], center[2] - half_height };
            auto tdisc1_cand = intersectDiscZForward(p, centerDisc, radii);
            centerDisc[2] = center[2] + half_height;
            auto tdisc2_cand = intersectDiscZForward(p, centerDisc, radii);
            return std::min({ t_cand, tdisc1_cand, tdisc2_cand }, [](const auto& lh, const auto& rh) -> bool { return lh.value_or(std::numeric_limits<T>::max()) < rh.value_or(std::numeric_limits<T>::max()); });
        }
    }
}
}