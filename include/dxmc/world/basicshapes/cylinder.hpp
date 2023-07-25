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
#include "dxmc/floating.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/visualizationintersectionresult.hpp"
#include "dxmc/world/worldintersectionresult.hpp"

#include <array>
#include <optional>

namespace dxmc {
namespace basicshape {
    namespace cylinder {

        template <Floating T>
        struct Cylinder {
            std::array<T, 3> center = { 0, 0, 0 };
            std::array<T, 3> direction = { 0, 0, 1 };
            T radius = 1;
            T half_height = 1;
            Cylinder() = default;
            Cylinder(std::array<T, 3> center_arr, std::array<T, 3> direction_arr, T radii, T half_height_wall)
                : radius(radii)
                , half_height(half_height_wall)
                , center(center_arr)
                , direction(direction_arr)
            {
            }
        };

        template <Floating T>
        constexpr inline bool isOverPlane(const std::array<T, 3>& planepoint, const std::array<T, 3>& planenormal, const std::array<T, 3>& point) noexcept
        {
            return vectormath::dot(vectormath::subtract(point, planepoint), planenormal) >= 0;
        }

        template <Floating T>
        bool pointInside(const std::array<T, 3>& pos, const Cylinder<T>& cylinder)
        {
            // test if point inside infinite cylinder

            const auto d = vectormath::cross(cylinder.direction, vectormath::subtract(pos, cylinder.center));
            if (vectormath::lenght_sqr(d) <= cylinder.radius * cylinder.radius) {
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

        template <Floating T>
        std::optional<std::array<T, 2>> intersectInterval(const Particle<T>& p, const Cylinder<T>& cylinder)
        {
            // return line segment cylinder wall intersection
            // intersection may be behind line start
            const auto e = vectormath::scale(cylinder.direction, cylinder.half_height);
            const auto p0 = vectormath::subtract(cylinder.center, e);
            const auto p1 = vectormath::add(cylinder.center, e);

            std::optional<T> tplane0;
            { // testing cylindar planes
                std::optional<T> tplane1;
                const auto planar = vectormath::dot(cylinder.direction, p.dir);
                if (std::abs(planar) > GEOMETRIC_ERROR<T>()) {
                    const auto den_inv = T { 1 } / planar;
                    const auto t0 = vectormath::dot(vectormath::subtract(p0, p.pos), cylinder.direction) * den_inv;
                    const auto tp0 = vectormath::add(p.pos, vectormath::scale(p.dir, t0));

                    const auto r2 = cylinder.radius * cylinder.radius;
                    if (vectormath::lenght_sqr(vectormath::subtract(tp0, p0)) <= r2)
                        tplane0 = t0;
                    const auto t1 = vectormath::dot(vectormath::subtract(p1, p.pos), cylinder.direction) * den_inv;
                    const auto tp1 = vectormath::add(p.pos, vectormath::scale(p.dir, t1));
                    if (vectormath::lenght_sqr(vectormath::subtract(tp1, p1)) <= r2)
                        tplane1 = t1;
                    // sorting hits
                    if (tplane0 && tplane1) {
                        // we hit both planes and early exits
                        if (tplane0 > tplane1) {
                            tplane0.swap(tplane1);
                        }
                        std::array<T, 2> t_planes = { tplane0.value(), tplane1.value() };
                        return std::make_optional(t_planes);
                    } else if (tplane1) {
                        // if not exit and one intersection, intersection is on tplane 0
                        tplane0.swap(tplane1);
                    }
                }
            }

            { // testing cylinder wall
                const auto v0 = vectormath::cross(vectormath::subtract(p.pos, cylinder.center), cylinder.direction);
                const auto v1 = vectormath::cross(p.dir, cylinder.direction);

                const auto a = vectormath::dot(v1, v1);
                if (a > GEOMETRIC_ERROR<T>()) {
                    const auto b = 2 * vectormath::dot(v0, v1);
                    const auto c = vectormath::dot(v0, v0) - cylinder.radius * cylinder.radius;
                    const auto den = b * b - 4 * a * c;
                    if (den > 0) {
                        const auto den_s = std::sqrt(den);
                        const auto a2inv = 1 / (2 * a);
                        std::array t_wall = { (-b - den_s) * a2inv, (-b + den_s) * a2inv };
                        // we have intersections on infinite cylinder wall, testing intersection points inside cylinder
                        // intersection points
                        const auto pt0 = vectormath::add(p.pos, vectormath::scale(p.dir, t_wall[0]));
                        const auto pt1 = vectormath::add(p.pos, vectormath::scale(p.dir, t_wall[1]));
                        if (!isOverPlane(p0, cylinder.direction, pt0) || isOverPlane(p1, cylinder.direction, pt0)) {
                            // we miss valid cylinder wall
                            if (tplane0) {
                                t_wall[0] = tplane0.value();
                                tplane0.reset();
                            } else { // no valid intersection
                                return std::nullopt;
                            }
                        }
                        if (!isOverPlane(p0, cylinder.direction, pt1) || isOverPlane(p1, cylinder.direction, pt1)) {
                            // we miss valid cylinder wall
                            if (tplane0) {
                                t_wall[1] = tplane0.value();
                            } else { // no valid intersection
                                return std::nullopt;
                            }
                        }
                        return std::make_optional(t_wall);
                    }
                }
            }
            return std::nullopt;
        }

        template <Floating T>
        std::optional<std::array<T, 2>> intersectForwardInterval(const Particle<T>& p, const Cylinder<T>& cylinder)
        {
            auto t_cand = intersectInterval(p, cylinder);
            if (t_cand) {
                auto& v = t_cand.value();
                if (v[1] < T { 0 })
                    return std::nullopt;
                v[0] = std::max(v[0], T { 0 });
            }
            return t_cand;
        }

        template <Floating T>
        WorldIntersectionResult<T> intersect(const Particle<T>& p, const Cylinder<T>& cylinder)
        {
            const auto t_cand = intersectInterval(p, cylinder);
            WorldIntersectionResult<T> res;
            if (t_cand) {
                const auto& v = t_cand.value();
                if (v[1] > T { 0 }) {
                    res.rayOriginIsInsideItem = v[0] < T { 0 };
                    res.intersection = res.rayOriginIsInsideItem ? v[1] : v[0];
                    res.intersectionValid = true;
                }
            }
            return res;
        }

        template <Floating T, typename U>
        VisualizationIntersectionResult<T, U> intersectVisualization(const Particle<T>& p, const Cylinder<T>& cylinder)
        {

            const auto t_cand = intersectInterval(p, cylinder);
            VisualizationIntersectionResult<T, U> res;
            if (t_cand) {
                const auto& v = t_cand.value();
                if (v[1] > T { 0 }) {
                    res.rayOriginIsInsideItem = v[0] < T { 0 };
                    res.intersection = res.rayOriginIsInsideItem ? v[1] : v[0];
                    res.intersectionValid = true;

                    // finding normal
                    const auto p0 = vectormath::add(p.pos, vectormath::scale(p.dir, res.intersection));
                    const auto pa = vectormath::subtract(cylinder.center, p0);
                    const auto d = vectormath::cross(pa, cylinder.direction);
                    const auto r = cylinder.radius * (1 - GEOMETRIC_ERROR<T>());
                    if (vectormath::lenght_sqr(d) < r * r) {
                        // we hit plane
                        if (vectormath::dot(pa, cylinder.direction) < 0) {
                            res.normal = vectormath::scale(cylinder.direction, T { -1 });
                        } else {
                            res.normal = cylinder.direction;
                        }
                    } else {
                        res.normal = vectormath::scale(d, T { -1 });
                    }
                }
            }
            return res;
        }
    }
}
}