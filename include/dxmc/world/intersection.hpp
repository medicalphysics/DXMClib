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

#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"

#include <limits>
#include <optional>

namespace dxmc {

template <Floating T>
std::optional<std::array<T, 2>> intersectAABB(const Particle<T>& p, const std::array<T, 6>& aabb)
{
    std::array<T, 2> t;

    const auto dx = T { 1 } / p.dir[0];
    if (dx >= 0) {
        t[0] = (aabb[0] - p.pos[0]) * dx;
        t[1] = (aabb[3] - p.pos[0]) * dx;
    } else {
        t[0] = (aabb[3] - p.pos[0]) * dx;
        t[1] = (aabb[0] - p.pos[0]) * dx;
    }
    T tymin, tymax;
    const auto dy = T { 1 } / p.dir[1];
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
    const auto dz = T { 1 } / p.dir[2];
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

template <Floating T>
bool pointInsideAABB(const std::array<T, 3>& p, const std::array<T, 6>& aabb)
{
    return aabb[0] <= p[0] && p[0] <= aabb[3] && aabb[1] <= p[1] && p[1] <= aabb[4] && aabb[2] <= p[2] && p[2] <= aabb[5];
}

template <Floating T>
std::optional<T> intersectAABBForward(const Particle<T>& p, const std::array<T, 6>& aabb)
{
    const auto t_cand = intersectAABB(p, aabb);
    if (t_cand) {
        const auto& t = t_cand.value();
        return pointInsideAABB(p.pos, aabb) ? std::make_optional(t[1]) : std::make_optional(t[0]);
    }
    return std::nullopt;
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
std::optional<std::array<T, 2>> intersectCylinderZ(const Particle<T>& p, const std::array<T, 3>& center, const T radii, const T half_height)
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
        return t_cand;
    } else {
        // no cylindar wall intersection, test for ends
        const auto tdl_cand = intersectDiscZ(p, centerDiscl, radii);
        const auto tdh_cand = intersectDiscZ(p, centerDisch, radii);
        if (tdl_cand && tdh_cand) {
            if (tdl_cand < tdh_cand) {
                std::array t = { tdl_cand.value(), tdh_cand.value() };
                return std::make_optional(t);
            } else {
                std::array t = { tdh_cand.value(), tdl_cand.value() };
                return std::make_optional(t);
            }
        }
    }
    return std::nullopt;
}

template <Floating T>
std::optional<T> intersectCylinderZForward(const Particle<T>& p, const std::array<T, 3>& center, const T radii, const T half_height)
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
                t = std::min(t, tdl_cand.value());
            }

            const std::array centerDisch = { center[0], center[1], center[2] + half_height };
            const auto tdh_cand = intersectDiscZ(p, centerDisch, radii);
            if (tdh_cand) {
                t = std::min(t, tdh_cand.value());
            }

            return t_cand;
        }
    } else {
        const std::array centerDiscl = { center[0], center[1], center[2] - half_height };
        const auto tdl_cand = intersectDiscZ(p, centerDiscl, radii);
        const std::array centerDisch = { center[0], center[1], center[2] + half_height };
        if (!tdl_cand) {
            return intersectDiscZ(p, centerDisch, radii);
        }

        const auto tdh_cand = intersectDiscZ(p, centerDisch, radii);
        if (!tdh_cand)
            return tdl_cand;
        return std::min(tdl_cand, tdh_cand);
    }
}

template <Floating T>
std::optional<T> intersectSphereForward(const Particle<T>& p, const std::array<T, 3>& center, const T radii)
{

    // nummeric stable ray sphere intersection
    const auto r2 = radii * radii;
    const auto f = vectormath::subtract(p.pos, center);

    // positive b mean center of sphere is in front of ray
    const auto b = -vectormath::dot(f, p.dir);

    // positive c means ray starts outside of sphere
    const auto c = vectormath::lenght_sqr(f) - r2;

    if ((c > 0) && (b < 0)) {
        // if ray starts outside sphere and center is begind ray
        // we exit early
        return std::nullopt;
    }

    const auto delta1 = vectormath::lenght_sqr(vectormath::add(f, vectormath::scale(p.dir, b)));

    const auto delta = r2 - delta1;

    if (delta < 0) {
        // no solution to the quadratic equation (we miss)
        return std::nullopt;
    }

    const int sign = b > 0 ? 1 : -1;
    const auto q = b + sign * std::sqrt(delta);

    if (c < 0 && b > 0) {
        // inside sphere
        return std::make_optional(q);
    } else {
        return std::make_optional(c / q);
    }
}

template <Floating T>
std::optional<std::array<T, 2>> intersectSphere(const Particle<T>& p, const std::array<T, 3>& center, const T radii)
{
    // nummeric stable ray sphere intersection
    const auto r2 = radii * radii;
    const auto f = vectormath::subtract(p.pos, center);

    // positive b mean center of sphere is in front of ray
    const auto b = -vectormath::dot(f, p.dir);

    // positive c means ray starts outside of sphere
    const auto c = vectormath::lenght_sqr(f) - r2;

    if ((c > 0) && (b < 0)) {
        // if ray starts outside sphere and center is begind ray
        // we exit early
        return std::nullopt;
    }

    const auto delta1 = vectormath::lenght_sqr(vectormath::add(f, vectormath::scale(p.dir, b)));

    const auto delta = r2 - delta1;

    if (delta < 0) {
        // no solution to the quadratic equation (we miss)
        return std::nullopt;
    }

    const int sign = b > 0 ? 1 : -1;
    const auto q = b + sign * std::sqrt(delta);
    if (c < 0 && b < 0) {
        std::array t = { q, c / q };
        return std::make_optional(t);
    } else {
        std::array t = { c / q, q };
        return std::make_optional(t);
    }
}
}
