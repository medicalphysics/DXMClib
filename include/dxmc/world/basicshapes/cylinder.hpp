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

        inline std::array<double, 6> cylinderAABB(const Cylinder& cyl)
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

        constexpr inline bool isOverPlane(const std::array<double, 3>& planepoint, const std::array<double, 3>& planenormal, const std::array<double, 3>& point) noexcept
        {
            return vectormath::dot(vectormath::subtract(point, planepoint), planenormal) >= 0;
        }

        bool pointInside(const std::array<double, 3>& pos, const Cylinder& cylinder)
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

        std::optional<double> intersectDisc(const ParticleType auto& p, const std::array<double, 3>& center, const std::array<double, 3>& normal, double radius)
        {
            const auto D = vectormath::dot(p.dir, normal);
            constexpr double minOrt = GEOMETRIC_ERROR();
            if (D < minOrt && D > -minOrt)
                return std::nullopt; // dir and normal is orthogonal, we exits
            auto t = vectormath::dot(vectormath::subtract(center, p.pos), normal) / D;

            // intersection point
            const auto p_int = vectormath::add(p.pos, vectormath::scale(p.dir, t));

            // distance from center
            const auto c_dist = vectormath::subtract(center, p_int);
            // check if distance from center is less than radius
            if (vectormath::dot(c_dist, c_dist) <= radius * radius) {
                return t;
            }
            return std::nullopt;
        }

        std::optional<std::array<double, 2>> intersectInterval(const ParticleType auto& p, const Cylinder& cylinder)
        {
            // return line segment cylinder wall intersection
            // intersection may be behind line start
            const auto e = vectormath::scale(cylinder.direction, cylinder.half_height);
            const auto p0 = vectormath::subtract(cylinder.center, e);
            const auto p1 = vectormath::add(cylinder.center, e);

            std::optional<double> tplane0;
            {
                const auto center0 = vectormath::add(cylinder.center, vectormath::scale(cylinder.direction, cylinder.half_height));
                const auto center1 = vectormath::add(cylinder.center, vectormath::scale(cylinder.direction, -cylinder.half_height));
                tplane0 = intersectDisc(p, center0, cylinder.direction, cylinder.radius);
                if (!tplane0) {
                    tplane0 = intersectDisc(p, center1, cylinder.direction, cylinder.radius);
                } else {
                    std::optional<double> tplane1 = intersectDisc(p, center1, cylinder.direction, cylinder.radius);
                    if (tplane0 && tplane1) {
                        // early exit
                        const auto [mi, ma] = std::minmax(*tplane0, *tplane1);
                        return std::make_optional<std::array<double, 2>>({ mi, ma });
                    }
                }
            }

            { // testing cylinder wall
                const auto v0 = vectormath::cross(vectormath::subtract(p.pos, cylinder.center), cylinder.direction);
                const auto v1 = vectormath::cross(p.dir, cylinder.direction);

                const auto a = vectormath::dot(v1, v1);
                if (a > GEOMETRIC_ERROR()) {
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

        std::optional<std::array<double, 2>> intersectForwardInterval(const ParticleType auto& p, const Cylinder& cylinder)
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

        WorldIntersectionResult intersect(const ParticleType auto& p, const Cylinder& cylinder)
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
                    const auto pa = vectormath::subtract(p0, cylinder.center);
                    const auto ns = vectormath::dot(pa, cylinder.direction);
                    const auto d = vectormath::subtract(pa, vectormath::scale(cylinder.direction, ns));
                    const auto r = cylinder.radius * (1 - GEOMETRIC_ERROR());
                    if (vectormath::length_sqr(d) < r * r) {
                        // we hit plane
                        if (vectormath::dot(pa, cylinder.direction) < 0) {
                            res.normal = vectormath::scale(cylinder.direction, -1.0);
                        } else {
                            res.normal = cylinder.direction;
                        }
                    } else {
                        res.normal = d;
                        vectormath::normalize(res.normal);
                    }
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