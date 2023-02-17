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
        std::optional<std::array<T, 2>> intersect(const Particle<T>& p, const std::array<T, 3>& center, const T radii, const T half_height)
        {
            auto t_cand = intersectCylinderWall(p, center, radii);
            const std::array centerDiscl = { center[0], center[1], center[2] - half_height };
            const std::array centerDisch = { center[0], center[1], center[2] + half_height };
            if (t_cand) {
                auto& t = t_cand.value();
                const auto pz0 = p.pos[2] + p.dir[2] * t[0];
                const auto pz1 = p.pos[2] + p.dir[2] * t[1];
                const bool cylwall0 = centerDiscl[2] < pz0 && pz0 < centerDisch[2];
                const bool cylwall1 = centerDiscl[2] < pz1 && pz1 < centerDisch[2];

                if (cylwall0 && cylwall1) {
                    // we intersect cylindar wall on valid z interval
                    return t_cand;
                }

                const auto tdl_cand = intersectDiscZ(p, centerDiscl, radii);
                if (tdl_cand) {
                    if (!cylwall0) {
                        t[0] = std::max(t[0], tdl_cand.value());
                    }
                    if (!cylwall1) {
                        t[1] = std::min(t[1], tdl_cand.value());
                    }
                }
                const auto tdh_cand = intersectDiscZ(p, centerDisch, radii);
                if (tdh_cand) {
                    if (!cylwall0) {
                        t[0] = std::max(t[0], tdh_cand.value());
                    }
                    if (!cylwall1) {
                        t[1] = std::min(t[1], tdh_cand.value());
                    }
                }
                return t[1] > 0 ? t_cand : std::nullopt;
            } else {
                // no cylindar wall intersection, test for ends
                const auto tdl_cand = intersectDiscZ(p, centerDiscl, radii);
                const auto tdh_cand = intersectDiscZ(p, centerDisch, radii);
                if (tdl_cand && tdh_cand) {
                    if (tdl_cand < tdh_cand) {
                        if (*tdh_cand > 0) {
                            std::array t = { tdl_cand.value(), tdh_cand.value() };
                            return std::make_optional(t);
                        }
                    } else {
                        if (*tdl_cand > 0) {
                            std::array t = { tdh_cand.value(), tdl_cand.value() };
                            return std::make_optional(t);
                        }
                    }
                }
            }
            return std::nullopt;
        }

        template <Floating T>
        std::optional<T> intersectForward(const Particle<T>& p, const std::array<T, 3>& center, const T radii, const T half_height)
        {
            auto t_cand = intersectCylinderWallForward(p, center, radii);
            if (t_cand) {
                // valid range?
                auto& t = t_cand.value();
                const auto pz = p.pos[2] + p.dir[2] * t;
                if (center[2] - half_height < pz && pz < center[2] - half_height) {
                    return t_cand;
                } else {
                    const std::array centerDiscl = { center[0], center[1], center[2] - half_height };
                    const auto tdl_cand = intersectDiscZ(p, centerDiscl, radii);
                    if (tdl_cand) {
                        t = tdl_cand.value();
                    }
                    const std::array centerDisch = { center[0], center[1], center[2] + half_height };
                    const auto tdh_cand = intersectDiscZ(p, centerDisch, radii);
                    if (tdh_cand) {
                        t = std::min(t, tdh_cand.value());
                    }

                    return t > 0 ? t_cand : std::nullopt;
                }
            } else {
                const std::array centerDiscl = { center[0], center[1], center[2] - half_height };
                const auto tdl_cand = intersectDiscZ(p, centerDiscl, radii);
                const std::array centerDisch = { center[0], center[1], center[2] + half_height };
                const auto tdh_cand = intersectDiscZ(p, centerDisch, radii);
                if (tdl_cand && tdh_cand) {
                    if (*tdl_cand > 0 && *tdh_cand > 0) {
                        return std::min(tdl_cand, tdh_cand);
                    } else {
                        const auto& td_cand = std::max(tdl_cand, tdh_cand);
                        return *td_cand > 0 ? td_cand : std::nullopt;
                    }
                } else {
                    return std::nullopt;
                    // this is unreacable
                }
            }
        }

    }
}
}