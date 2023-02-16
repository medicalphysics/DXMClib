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
    namespace sphere {

        template <Floating T>
        std::optional<T> intersectForward(const Particle<T>& p, const std::array<T, 3>& center, const T radii)
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
        std::optional<std::array<T, 2>> intersect(const Particle<T>& p, const std::array<T, 3>& center, const T radii)
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
}
}