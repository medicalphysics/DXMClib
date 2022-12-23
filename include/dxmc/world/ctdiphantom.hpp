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
#include "dxmc/vectormath.hpp"

#include <execution>
#include <limits>
#include <optional>

namespace dxmc {

template <Floating T>
class CTDIPhantom {
public:
    CTDIPhantom(T radius = T { 16 }, const std::array<T, 3>& pos = { 0, 0, 0 }, T height = T { 15 })
        : m_radius(radius)
        , m_halfheight(height * T { 0.5 })
        , m_center(pos)
    {
    }

    void translate(const std::array<T, 3>& dist)
    {
        m_center = vectormath::add(m_center, dist);
    }
    const std::array<T, 3>& center() const
    {
        return m_center;
    }

    std::array<T, 6> AABB() const
    {
        std::array<T, 3> llc = vectormath::add(m_center, -m_radius);
        llc[2] = m_center[2] - m_halfheight;
        std::array<T, 3> urc = vectormath::add(m_center, m_radius);
        urc[2] = m_center[2] + m_halfheight;
        return vectormath::join(llc, urc);
    }

    std::optional<T> intersect(const Particle<T>& p) const
    {
        constexpr std::array<T, 2> tbox { 0, std::numeric_limits<T>::max() };
        return intersectCylinder(p, m_center, m_radius, m_center[2] - m_halfheight, m_center[2] + m_halfheight, tbox);
    }
    std::optional<T> intersect(const Particle<T>& p, const std::array<T, 2>& tbox) const
    {
        return intersectCylinder(p, m_center, m_radius, m_center[2] - m_halfheight, m_center[2] + m_halfheight, tbox);
    }
    T transport(Particle<T>& p, RandomState& state)
    {
        return 0;
    }

protected:
    static std::optional<T> intersectCylinder(const Particle<T>& p, const std::array<T, 3>& center, const T radii, const T zStart, const T zStop, const std::array<T, 2>& tbox)
    {
        bool has_value = false;
        T t_min;

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
                                t_min = c / q;
                                has_value = true;
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
                            t_min = t_cand;
                            has_value = true;
                        }
                    }
                }
            }
        }

        const auto t_disc1 = intersectDisc(p, center, radii, zStart, tbox);
        if (t_disc1) {
            const auto t_cand = t_disc1.value();
            t_min = has_value ? std::min(t_min, t_cand) : t_cand;
            has_value = true;
        }
        const auto t_disc2 = intersectDisc(p, center, radii, zStop, tbox);
        if (t_disc2) {
            const auto t_cand = t_disc2.value();
            t_min = has_value ? std::min(t_min, t_cand) : t_cand;
            has_value = true;
        }

        return has_value ? std::make_optional(t_min) : std::nullopt;
    }

    static std::optional<T> intersectDisc(const Particle<T>& p, const std::array<T, 3>& center, const T radii, const T z, const std::array<T, 2>& tbox)
    {
        if (std::abs(p.dir[2]) <= std::numeric_limits<T>::epsilon())
            return std::nullopt;
        const auto tz = (z - p.pos[2]) / p.dir[2];
        if (tbox[0] < tz && tz < tbox[1]) {
            const auto xz = p.pos[0] + p.dir[0] * tz - center[0];
            const auto yz = p.pos[1] + p.dir[1] * tz - center[1];
            if (xz * xz + yz * yz < radii * radii)
                return tz;
        }
        return std::nullopt;
    }

private:
    T m_radius = 0;
    T m_halfheight = 0;
    std::array<T, 3> m_center;
};

}
