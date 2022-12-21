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
        std::array<std::optional<T>, 4> t_cand;

        // double for precision, may get errors using float
        // Fix later numerical stable squared solver???
        const auto a = static_cast<double>(p.dir[0]) * static_cast<double>(p.dir[0]) + static_cast<double>(p.dir[1]) * static_cast<double>(p.dir[1]);
        const auto dx = static_cast<double>(p.pos[0]) - static_cast<double>(center[0]);
        const auto dy = static_cast<double>(p.pos[1]) - static_cast<double>(center[1]);
        const auto b = 2 * (dx * static_cast<double>(p.dir[0]) + dy * static_cast<double>(p.dir[1]));
        const auto c = dx * dx + dy * dy - static_cast<double>(radii) * static_cast<double>(radii);
        const auto det = b * b - 4 * a * c;

        // no intersection to 2d circle, treats tangents as no intersection
        if (det < T { 0 } || std::abs(a) <= std::numeric_limits<T>::epsilon()) {
            // test if we intersect end discs
            t_cand[0] = std::nullopt;
            t_cand[1] = std::nullopt;
        } else {
            const auto sqrt_det = std::sqrt(det);
            const auto den = double { 0.5 } / a;
            const auto x1 = b > 0 ? (-b - sqrt_det) * den : (-b + sqrt_det) * den;
            const T t1 = static_cast<T>(x1);
            if (tbox[0] < t1 && t1 < tbox[1]) {
                const auto int_pointZ = p.pos[2] + p.dir[2] * t1;
                const auto valid = zStart < int_pointZ && int_pointZ < zStop;
                if (valid)
                    t_cand[0] = t1;
            }
            const auto x2 = c / (a * x1);
            const T t2 = static_cast<T>(x2);
            if (tbox[0] < t2 && t2 < tbox[1]) {
                const auto int_pointZ = p.pos[2] + p.dir[2] * t2;
                const auto valid = zStart < int_pointZ && int_pointZ < zStop;
                if (valid)
                    t_cand[1] = t2;
            }
        }
        t_cand[2] = intersectDisc(p, center, radii, zStart, tbox);
        t_cand[3] = intersectDisc(p, center, radii, zStop, tbox);
        auto t = std::reduce(std::execution::unseq, t_cand.begin() + 1, t_cand.end(), t_cand[0], [](const std::optional<T>& lh, const std::optional<T>& rh) -> std::optional<T> {
            if (lh && rh)
                return lh < rh ? lh : rh;
            else if (lh)
                return lh;
            return rh;
        });
        if (t)
            return t;
        else
            return std::nullopt;
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
