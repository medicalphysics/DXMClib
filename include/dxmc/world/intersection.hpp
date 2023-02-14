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
        tzmin = (aabb[2] - p.pos[2]) * dz;
        tzmax = (aabb[5] - p.pos[2]) * dz;
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
std::optional<T> intersectCylinderWall(const Particle<T>& p, const std::array<T, 3>& center, const T radii, const T zStart, const T zStop, const std::array<T, 2>& tbox)
{
    // Nummerical stable cylindar ray intersection (not stable enought) should attempt to normalize d to unity
    const auto dx = p.dir[0];
    const auto dy = p.dir[1];
    const auto a = dx * dx + dy * dy;
    if (std::abs(a) > std::numeric_limits<T>::epsilon()) {
        // ray not parallell to z axis
        const auto fx = p.pos[0] - center[0];
        const auto fy = p.pos[1] - center[1];
        const auto b = -(fx * dx + fy * dy);
        const auto ba = b / a;
        const auto p1x = fx + ba * dx;
        const auto p1y = fy + ba * dy;
        const auto r2 = radii * radii;
        const auto det = r2 - (p1x * p1x + p1y * p1y);
        if (det > T { 0 }) {
            // line intersect 2d sphere
            const auto c = fx * fx + fy * fy - r2;
            const int sign = (T { 0 } < b) - (b < T { 0 });
            const auto q = b + sign * std::sqrt(a * det);
            if (c > T { 0 }) {
                // ray starts outside sphere
                if (b > T { 0 }) {
                    // ray starts before sphere
                    const auto t_cand = c / q;
                    if (tbox[0] < t_cand && t_cand < tbox[1]) {
                        const auto int_pointZ = p.pos[2] + p.dir[2] * t_cand;
                        const auto valid = zStart < int_pointZ && int_pointZ < zStop;
                        if (valid) {
                            return std::make_optional<T>(t_cand);
                        }
                    }
                }
            } else {
                // ray is inside
                const auto t_cand = q / a;
                if (tbox[0] < t_cand && t_cand < tbox[1]) {
                    // intersection is forward
                    const auto int_pointZ = p.pos[2] + p.dir[2] * t_cand;
                    const auto valid = zStart < int_pointZ && int_pointZ < zStop;
                    if (valid) {
                        return std::make_optional<T>(t_cand);
                    }
                }
            }
        }
    }
    return std::nullopt;
}

template <Floating T>
std::optional<T> intersectDiscZ(const Particle<T>& p, const std::array<T, 3>& center, const T radii, const T z)
{
    if (std::abs(p.dir[2]) <= std::numeric_limits<T>::epsilon())
        return std::nullopt;
    const auto tz = (z - p.pos[2]) / p.dir[2];

    const auto xz = p.pos[0] + p.dir[0] * tz - center[0];
    const auto yz = p.pos[1] + p.dir[1] * tz - center[1];
    if (xz * xz + yz * yz < radii * radii)
        return std::make_optional(tz);

    return std::nullopt;
}

template <Floating T>
std::optional<T> intersectCylinderZ(const Particle<T>& p, const std::array<T, 3>& center, const T radii, const T half_height)
{
    std::array<T, 6> aabb = {
        center[0] - radii,
        center[1] - radii,
        center[2] - half_height,
        center[0] + radii,
        center[1] + radii,
        center[2] + half_height
    };

    const auto tbox_cand = intersectAABB(p, aabb);
    if (!tbox_cand)
        return std::nullopt;

    const auto& tbox = tbox_cand.value();

    auto t_wall = intersectCylinderWall(p, center, radii, aabb[2], aabb[5], tbox);
    std::optional<T> t_disc;

    if (p.dir[2] > T { 0 }) {
        if (p.pos[2] < aabb[2])
            t_disc = intersectDiscZ(p, center, radii, aabb[2]);
        else
            t_disc = intersectDiscZ(p, center, radii, aabb[5]);
    } else {
        if (p.pos[2] > aabb[5])
            t_disc = intersectDiscZ(p, center, radii, aabb[5]);
        else
            t_disc = intersectDiscZ(p, center, radii, aabb[2]);
    }

    if (t_wall && t_disc)
        return *t_wall < *t_disc ? t_wall : t_disc;
    else
        return t_disc ? t_disc : t_wall;
}

template <Floating T>
std::optional<T> intersectSphereForward(const Particle<T>& p, const std::array<T, 3>& center, const T radii)
{
    /* std::array<T, 6> aabb = {
        center[0] - radii,
        center[1] - radii,
        center[2] - radii,
        center[0] + radii,
        center[1] + radii,
        center[2] + radii
    };

    const auto tbox_cand = intersectAABB(p, aabb);
    if (!tbox_cand)
        return std::nullopt;
    const auto& tbox = tbox_cand.value();

    std::optional<T> t_cand = std::nullopt;

    const auto f = vectormath::subtract(p.pos, center);
    const auto b = -vectormath::dot(f, p.dir);
    const auto p1 = vectormath::add(f, vectormath::scale(p.dir, b));
    const auto r2 = radii * radii;
    const auto det = r2 - vectormath::dot(p1, p1);
    if (det > T { 0 }) {
        const auto c = vectormath::dot(f, f) - r2;
        if (c > T { 0 }) {
            const int sign = (T { 0 } < b) - (b < T { 0 });
            const auto q = b + sign * std::sqrt(det);
            const auto t = c / q;
            if ((tbox[0] < t) && (t < tbox[1]))
                t_cand = t;
        } else {
            const int sign = (T { 0 } < b) - (b < T { 0 });
            const auto q = b + sign * std::sqrt(det);
            const auto t = q;
            if ((tbox[0] < t) && (t < tbox[1]))
                t_cand = t;
        }
    }

    return t_cand;
    */

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

    if (c > 0) {
        //outside sphere
        return std::make_optional(c / q);
    } else {
        return std::make_optional(q);
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
    feil i rekkefølgen
    std::array t = { c / q, q };
    return std::make_optional(t);
}
}
