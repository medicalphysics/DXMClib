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

#include <array>
#include <optional>

namespace dxmc {
namespace basicshape {
    namespace tetrahedron {

        template <Floating T>
        bool planeSidePositive(const std::array<T, 3>& po, const std::array<T, 3>& v0, const std::array<T, 3>& v1)
        {
            // v0 and v1 spans the plane
            // po is vector from ray pos and a point on the plane

            // (v0 x op) * (v0 x v1) = (v0*v0)*(op*v1) - (v0*v1)*(op*v0)
            //                       =       1*(op*v1) - (v0*v1)*(op*v0)
            return vectormath::dot(po, v1) - vectormath::dot(v0, v1) * vectormath::dot(po, v0) > 0;
        }

        template <Floating T>
        bool pointInside(const std::array<T, 3>& pos,
            const std::array<T, 3>& v0,
            const std::array<T, 3>& v1,
            const std::array<T, 3>& v2,
            const std::array<T, 3>& v3)
        {
            const auto v0v1 = vectormath::subtract(v1, v0);
            const auto v1v2 = vectormath::subtract(v2, v1);
            const auto v1v3 = vectormath::subtract(v3, v1);
            const auto v2v3 = vectormath::subtract(v3, v2);
            const auto v0v3 = vectormath::subtract(v3, v0);
            const auto v3v2 = vectormath::subtract(v2, v3); // negative duplicate

            const auto v0p = vectormath::subtract(pos, v0);
            const auto v1p = vectormath::subtract(pos, v0);

            const auto pd0 = !planeSidePositive(v0p, v0v1, v1v2);
            const auto pd1 = !planeSidePositive(v0p, v0v1, v1v3);
            const auto pd2 = !planeSidePositive(v1p, v1v2, v2v3);
            const auto pd3 = !planeSidePositive(v0p, v0v3, v3v2);
            return pd0 && pd1 && pd2 && pd3;
        }

        template<Floating T>
        inline T permuted_inner_product(const std::array<std::array<T, 3>, 2>& lh, const std::array<std::array<T, 3>, 2>& rh)
        {
            // triple product a*(b x c)
            return vectormath::dot(lh[0], rh[1]) + vectormath::dot(lh[1], rh[0]);
        }

        template <Floating T>
        WorldIntersectionResult<T> intersect(const Particle<T>& p,
            const std::array<T, 3>& v0,
            const std::array<T, 3>& v1,
            const std::array<T, 3>& v2,
            const std::array<T, 3>& v3)
        {
            // intersection uing plucker coordinates

            const std::array<std::array<T, 3>, 2> p = { particle.dir, vectormath::cross(particle.dir, particle.pos) };

            // Faces: F3=(V0 V1 V2), F2=(V1 V0 V3), F1=(V2 V3 V0), F0=(V3 V2 V1)
            // Edges: E30=(V1 V2), E31=(V2 V0), E32=(V0 V1)
            // Edges: E20=(V0 V3), E21=(V3 V1), E22=(V1 V0)
            // Edges: E10=(V3 V0), E11=(V0 V2), E12=(V2 V3)
            // Edges: E00=(V2 V1), E01=(V1 V3), E02=(V3 V2)
            
            // Plucker coordinates for edges
            const auto e30 = vectormath::subtract(v2, v1);
            const std::array<std::array<T, 3>, 2> pe30 = { e30, vectormath::cross(e30, v1) };



            WorldIntersectionResult<T> res;
            return res;
        }

        template <Floating T>
        std::optional<std::array<T, 2>> intersectForwardInterval(const Particle<T>& p, const std::array<T, 3>& center, const T radii)
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

        template <Floating T, typename U>
        VisualizationIntersectionResult<T, U> intersectVisualization(const Particle<T>& p, const std::array<T, 3>& center, const T radii)
        {
            VisualizationIntersectionResult<T, U> res;
            const auto t_opt = intersectForwardInterval(p, center, radii);
            if (t_opt) {
                const auto& t = t_opt.value();
                res.intersectionValid = true;
                res.rayOriginIsInsideItem = t[0] < 0;
                res.intersection = res.rayOriginIsInsideItem ? t[1] : t[0];

                const auto pos = vectormath::add(p.pos, vectormath::scale(p.dir, res.intersection));
                res.normal = vectormath::subtract(center, pos);
                vectormath::normalize(res.normal);
            }
            return res;
        }
    }
}
}