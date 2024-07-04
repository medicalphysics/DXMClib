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

#include "dxmc/constants.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/visualizationintersectionresult.hpp"
#include "dxmc/world/worldintersectionresult.hpp"

#include <array>
#include <optional>

namespace dxmc {
namespace basicshape {
    namespace cylinder {

        struct Cylinder {
            std::array<double, 3> center = { 0, 0, 0 };
            std::array<double, 3> direction = { 0, 0, 1 };
            double radius = 1;
            double half_height = 1;
            Cylinder() = default;
            Cylinder(const std::array<double, 3>& center_arr, const std::array<double, 3>& direction_arr, double radii, double half_height_wall)
                : center(center_arr)
                , radius(radii)
                , half_height(half_height_wall)
            {
                direction = vectormath::normalized(direction_arr);
            }
            double volume() const
            {
                return std::numbers::pi_v<double> * radius * radius * half_height * 2;
            }
        };

        static std::array<double, 6> cylinderAABB(const Cylinder& cyl)
        {
            // calculating disc extents
            const std::array<double, 3> e = {
                cyl.radius * std::sqrt(1 - cyl.direction[0] * cyl.direction[0]),
                cyl.radius * std::sqrt(1 - cyl.direction[1] * cyl.direction[1]),
                cyl.radius * std::sqrt(1 - cyl.direction[2] * cyl.direction[2])
            };
            const auto l = vectormath::scale(cyl.direction, cyl.half_height);
            const auto p0 = vectormath::subtract(cyl.center, l);
            const auto p1 = vectormath::add(cyl.center, l);
            std::array<double, 6> aabb = {
                std::min(p0[0], p1[0]) - e[0],
                std::min(p0[1], p1[1]) - e[1],
                std::min(p0[2], p1[2]) - e[2],
                std::max(p0[0], p1[0]) + e[0],
                std::max(p0[1], p1[1]) + e[1],
                std::max(p0[2], p1[2]) + e[2]
            };
            return aabb;
        }

        static constexpr bool isOverPlane(const std::array<double, 3>& planepoint, const std::array<double, 3>& planenormal, const std::array<double, 3>& point) noexcept
        {
            return vectormath::dot(vectormath::subtract(point, planepoint), planenormal) >= 0;
        }

        static constexpr bool pointInside(const std::array<double, 3>& pos, const Cylinder& cylinder)
        {
            const auto d = vectormath::cross(cylinder.direction, vectormath::subtract(pos, cylinder.center));
            // test if point inside infinite cylinder
            if (vectormath::length_sqr(d) <= cylinder.radius * cylinder.radius) {
                // test for side of two end planes
                const auto e = vectormath::scale(cylinder.direction, cylinder.half_height);
                const auto p0 = vectormath::subtract(cylinder.center, e);

                if (vectormath::dot(vectormath::subtract(pos, p0), cylinder.direction) >= 0) {
                    // (pos - p)*normal >= 0
                    const auto p1 = vectormath::add(cylinder.center, e);
                    if (vectormath::dot(vectormath::subtract(pos, p1), cylinder.direction) <= 0) {
                        return true;
                    }
                }
            }
            return false;
        }

        static std::optional<std::array<double, 2>> intersectInterval(const ParticleType auto& p, const Cylinder& cylinder)
        {
            // infinite cylinder wall
            std::array<double, 2> t = { 0, 0 };
            const auto oc = vectormath::subtract(p.pos, cylinder.center);
            const auto ococ = vectormath::length_sqr(oc);
            const auto card = vectormath::dot(cylinder.direction, p.dir);
            const auto a = 1 - card * card;

            const auto caoc = vectormath::dot(cylinder.direction, oc);

            const auto b = vectormath::dot(oc, p.dir) - caoc * card;
            const auto c = ococ - caoc * caoc - cylinder.radius * cylinder.radius;
            const auto h2 = b * b - a * c;

            if (h2 < 0.0) // no intersection on wall or cap
                return std::nullopt;

            const auto card_inv = 1 / card;
            if (a < GEOMETRIC_ERROR<double>()) { // parallell ray, no wall intersect
                // we do an easy cap test
                const auto tc_1 = (-caoc + cylinder.half_height) * card_inv;
                const auto tc_2 = (-caoc - cylinder.half_height) * card_inv;
                t[0] = std::min(tc_1, tc_2);
                t[1] = std::max(tc_1, tc_2);
                return t[1] > 0.0 ? std::make_optional(t) : std::nullopt;
            }

            const auto h = std::sqrt(h2);
            const auto a_inv = 1 / a;
            t[0] = (-b - h) * a_inv;
            t[1] = (-b + h) * a_inv;

            const auto y0 = caoc + card * t[0];
            const auto y1 = caoc + card * t[1];
            if (-cylinder.half_height <= y0 && y0 <= cylinder.half_height) {
                if (-cylinder.half_height <= y1 && y1 <= cylinder.half_height)
                    // Both hits are walls
                    return std::make_optional(t);
                else {
                    // y1 failed
                    if (y1 > 0.0) // upper plane
                        t[1] = (-caoc + cylinder.half_height) * card_inv;
                    else // lower plane
                        t[1] = (-caoc - cylinder.half_height) * card_inv;
                    return std::make_optional(t);
                }
            } else {
                if (-cylinder.half_height <= y1 && y1 <= cylinder.half_height) {
                    // y0 failed
                    if (y0 > 0.0) // upper plane
                        t[0] = (-caoc + cylinder.half_height) * card_inv;
                    else // lower plane
                        t[0] = (-caoc - cylinder.half_height) * card_inv;
                    return std::make_optional(t);
                } else {
                    // both failed, we most likely miss
                    const auto tc1 = (-caoc - cylinder.half_height) * card_inv;
                    if (tc1 < t[0] || tc1 > t[1])
                        return std::nullopt;
                    const auto tc2 = (-caoc + cylinder.half_height) * card_inv;
                    if (tc2 < t[0] || tc2 > t[1])
                        return std::nullopt;
                    t[0] = std::min(tc1, tc2);
                    t[1] = std::max(tc1, tc2);
                    return std::make_optional(t);
                }
            }
            return std::nullopt;
        }

        static std::optional<std::array<double, 2>> intersectForwardInterval(const ParticleType auto& p, const Cylinder& cylinder)
        {
            auto t_cand = intersectInterval(p, cylinder);
            if (t_cand) {
                auto& v = t_cand.value();
                if (v[1] <= 0)
                    return std::nullopt;
                v[0] = std::max(v[0], 0.0);
            }
            return t_cand;
        }

        static WorldIntersectionResult intersect(const ParticleType auto& p, const Cylinder& cylinder)
        {
            const auto t_cand = intersectInterval(p, cylinder);
            WorldIntersectionResult res;
            if (t_cand) {
                const auto& v = t_cand.value();
                if (v[1] > 0) {
                    res.rayOriginIsInsideItem = v[0] <= 0;
                    res.intersection = res.rayOriginIsInsideItem ? v[1] : v[0];
                    res.intersectionValid = true;
                }
            }
            return res;
        }

        static std::array<double, 3> normalOnPoint(const std::array<double, 3>& pos, const Cylinder& cylinder)
        {
            // finding normal
            const auto pa = vectormath::subtract(pos, cylinder.center);
            const auto ns = vectormath::dot(pa, cylinder.direction);
            const auto d = vectormath::subtract(pa, vectormath::scale(cylinder.direction, ns));
            const auto r = cylinder.radius * (1 - GEOMETRIC_ERROR());
            std::array<double, 3> normal;
            if (vectormath::length_sqr(d) < r * r) {
                // we hit plane
                if (vectormath::dot(pa, cylinder.direction) < 0) {
                    normal = vectormath::scale(cylinder.direction, -1.0);
                } else {
                    normal = cylinder.direction;
                }
            } else {
                normal = vectormath::normalized(d);
            }
            return normal;
        }

        template <typename U>
        VisualizationIntersectionResult<U> intersectVisualization(const ParticleType auto& p, const Cylinder& cylinder)
        {
            const auto t_cand = intersectInterval(p, cylinder);
            VisualizationIntersectionResult<U> res;
            if (t_cand) {
                const auto& v = t_cand.value();
                if (v[1] > 0) {
                    res.rayOriginIsInsideItem = v[0] <= 0;
                    res.intersection = res.rayOriginIsInsideItem ? v[1] : v[0];
                    res.intersectionValid = true;

                    // finding normal
                    const auto p0 = vectormath::add(p.pos, vectormath::scale(p.dir, res.intersection));
                    res.normal = normalOnPoint(p0, cylinder);
                    if (res.rayOriginIsInsideItem) {
                        res.normal = vectormath::scale(res.normal, -1.0);
                    }
                }
            }
            return res;
        }
    }
}
}