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

#include "dxmc/floating.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/dxmcrandom.hpp"

#include <optional>
#include <execution>
#include <limits>

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
        return intersectCylinder(p, center(), m_radius, this->m_aabb[2], this->m_aabb[5], tbox);
    }
    std::optional<T> intersect(const Particle<T>& p, const std::array<T, 2>& tbox) const
    {
        return intersectCylinder(p, center(), m_radius, this->m_aabb[2], this->m_aabb[5], tbox);
    }
    T transport(Particle<T>& p, RandomState& state)
    {
        return 0;
    }

protected:
    static std::optional<T> intersectCylinder(const Particle<T>& p, const std::array<T, 3>& center, const T radii, const T zStart, const T zStop, const std::array<T, 2>& tbox)
    {
        const auto a = p.dir[0] * p.dir[0] + p.dir[1] * p.dir[1];
        const auto d0 = p.pos[0] - center[0];
        const auto d1 = p.pos[1] - center[1];
        const auto b = T { 2 } * (p.dir[0] * d0 + p.dir[1] * d1);
        const auto c = d0 * d0 + d1 * d1 - radii * radii;
        const auto det = b * b - T { 4 } * a * c;

        std::array<std::optional<T>, 4> t_cand;

        // no intersection to 2d circle, treats tangents as no intersection
        if (det <= std::numeric_limits<T>::epsilon()) {
            // test if we intersect end discs
            t_cand[0] = std::nullopt;
            t_cand[1] = std::nullopt;
        } else {
            const auto sqrt_det = std::sqrt(det);
            const auto den = T { 0.5 } / a;
            const auto t1 = (-b - sqrt_det) * den;
            if (tbox[0] < t1 && t1 < tbox[1]) {
                const auto int_pointZ = p.pos[2] + p.dir[2] * t1;
                const auto valid = zStart < int_pointZ && int_pointZ < zStop;
                if (valid)
                    t_cand[0] = t1;
            }
            const auto t2 = (-b + sqrt_det) * den;
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
        return t;
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
