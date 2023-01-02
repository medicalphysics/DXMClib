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
#include "dxmc/world/worlditembase.hpp"

#include <limits>
#include <optional>

namespace dxmc {

template <Floating T>
class Sphere final : public WorldItemBase<T> {
public:
    Sphere(T radius = T { 16 }, const std::array<T, 3>& pos = { 0, 0, 0 })
        : WorldItemBase<T>()
        , m_radius(radius)
        , m_center(pos)
    {
    }

    void translate(const std::array<T, 3>& dist) override
    {
        m_center = vectormath::add(m_center, dist);
    }
    std::array<T, 3> center() const override
    {
        return m_center;
    }

    std::array<T, 6> AABB() const override
    {
        std::array<T, 3> llc = vectormath::add(m_center, -m_radius);
        std::array<T, 3> urc = vectormath::add(m_center, m_radius);
        return vectormath::join(llc, urc);
    }
    IntersectionResult<T>
    intersect(const Particle<T>& p) const override
    {
        constexpr std::array<T, 2> tbox { 0, std::numeric_limits<T>::max() };
        return intersect(p, tbox);
    }
    IntersectionResult<T> intersect(const Particle<T>& p, const std::array<T, 2>& tbox) const override
    {

        const auto t = intersectSphere(p, m_center, m_radius, tbox);
        IntersectionResult<T> res;
        if (t) {
            res.item = this;
            res.intersection = t.value();
        }
        return res;
    }
    T transport(Particle<T>& p, RandomState& state) override
    {
        return 0;
    }

protected:
    static std::optional<T> intersectSphere(const Particle<T>& p, const std::array<T, 3>& center, const T radii, const std::array<T, 2>& tbox)
    {
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
    }

private:
    T m_radius = 0;
    std::array<T, 3> m_center;
};

}
