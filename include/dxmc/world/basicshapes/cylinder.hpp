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

        /*static std::array<std::optional<double>, 2> intersectDisc(const ParticleType auto& p, const Cylinder& cylinder)
        {
            std::array<std::optional<double>, 2> res;

            const auto denom = vectormath::dot(p.dir, cylinder.direction);
            const auto r2 = cylinder.radius * cylinder.radius;

            const auto c1 = vectormath::add(cylinder.center, vectormath::scale(cylinder.direction, cylinder.half_height));
            const auto dist1 = vectormath::subtract(c1, p.pos);
            const auto t1 = vectormath::dot(dist1, cylinder.direction) / denom;
            const auto hit_dist1 = vectormath::length_sqr(vectormath::subtract(vectormath::add(p.pos, vectormath::scale(p.dir, t1)), c1));
            res[0] = hit_dist1 <= r2 ? std::make_optional(t1) : std::nullopt;

            const auto c2 = vectormath::add(cylinder.center, vectormath::scale(cylinder.direction, -cylinder.half_height));
            const auto dist2 = vectormath::subtract(c2, p.pos);
            const auto t2 = vectormath::dot(dist2, cylinder.direction) / denom;
            const auto hit_dist2 = vectormath::length_sqr(vectormath::subtract(vectormath::add(p.pos, vectormath::scale(p.dir, t2)), c2));
            res[1] = hit_dist2 <= r2 ? std::make_optional(t2) : std::nullopt;

            return res;
        }

        static std::optional<std::array<double, 2>> intersectInterval(const ParticleType auto& p, const Cylinder& cylinder)
        {
            std::array<double, 2> t = { 0, 0 };

            const auto na = vectormath::cross(p.dir, cylinder.direction);

            if (vectormath::length_sqr(na) < GEOMETRIC_ERROR<>()) {
                // only test caps
                const auto res = intersectDisc(p, cylinder);
                if (res[0] && res[1]) {
                    t[0] = std::min(*res[0], *res[1]);
                    t[1] = std::max(*res[0], *res[1]);
                    return t;
                }
                return std::nullopt;
            }

            const auto b = vectormath::subtract(cylinder.center, p.pos);
            const auto ba = vectormath::cross(b, cylinder.direction);
            const auto naba = vectormath::dot(na, ba);
            const auto nana = vectormath::dot(na, na);
            const auto bna = vectormath::dot(b, na);

            const auto r2 = cylinder.radius * cylinder.radius;
            const auto rot_term = nana * r2 - bna * bna;
            if (rot_term <= 0.0)
                return std::nullopt;

            const auto term = std::sqrt(rot_term);
            const auto nana_inv = 1 / nana;

            // cylinder intersects
            const auto ct1 = (naba + term) * nana_inv;
            const auto ct2 = (naba - term) * nana_inv;
            t[0] = ct1 < ct2 ? ct1 : ct2;
            t[1] = ct1 > ct2 ? ct1 : ct2;

            const auto y1 = vectormath::dot(cylinder.direction, vectormath::subtract(vectormath::scale(p.dir, t[0]), b));
            const auto y2 = vectormath::dot(cylinder.direction, vectormath::subtract(vectormath::scale(p.dir, t[1]), b));
            if (-cylinder.half_height <= y1 && y1 <= cylinder.half_height && -cylinder.half_height <= y2 && y2 <= cylinder.half_height) {
                // we hit only cylinder return
                return t;
            } else if (-cylinder.half_height <= y1 && y1 <= cylinder.half_height) {
                const auto res = intersectDisc(p, cylinder);
                t[1] = res[0] ? *res[0] : *res[1];
                return t;
            } else if (-cylinder.half_height <= y2 && y2 <= cylinder.half_height) {
                const auto res = intersectDisc(p, cylinder);
                t[0] = res[0] ? *res[0] : *res[1];
                return t;
            } else {
                const auto res = intersectDisc(p, cylinder);
                if (res[0] && res[1]) {
                    t[0] = std::min(*res[0], *res[1]);
                    t[1] = std::max(*res[0], *res[1]);
                    return t;
                }
            }
            return std::nullopt;
        }*/

        static std::optional<std::array<double, 2>> intersectInterval(const ParticleType auto& p, const Cylinder& cylinder)
        {
            const auto b = vectormath::subtract(cylinder.center, p.pos);
            const auto na = vectormath::cross(p.dir, cylinder.direction);
            const auto nana = vectormath::length_sqr(na);
            const auto r2 = cylinder.radius * cylinder.radius;
            if (nana < GEOMETRIC_ERROR<double>()) {
                std::array<double, 2> t;

                const auto c1 = vectormath::add(b, vectormath::scale(cylinder.direction, -cylinder.half_height));
                t[0] = vectormath::dot(cylinder.direction, c1) / vectormath::dot(cylinder.direction, p.dir);
                const auto ndc1 = vectormath::subtract(vectormath::scale(p.dir, t[0]), c1);
                if (vectormath::dot(ndc1, ndc1) < r2) {
                    const auto c2 = vectormath::add(b, vectormath::scale(cylinder.direction, cylinder.half_height));
                    t[1] = vectormath::dot(cylinder.direction, c2) / vectormath::dot(cylinder.direction, p.dir);
                    const auto ndc2 = vectormath::subtract(vectormath::scale(p.dir, t[1]), c2);
                    if (vectormath::dot(ndc2, ndc2) < r2)
                        return t;
                }
                return std::nullopt;
                // endpoints only
            } else {
                const auto bna = vectormath::dot(b, na);
                const auto broot = nana * r2 - bna * bna;
                if (broot <= 0.0) {
                    // no intersection at all
                    return std::nullopt;
                }
                const auto naba = vectormath::dot(na, vectormath::cross(b, cylinder.direction));
                std::array<double, 2> t = {
                    (naba - broot) / nana,
                    (naba + broot) / nana
                };

                const std::array<double, 2> y = {
                    vectormath::dot(cylinder.direction, vectormath::subtract(vectormath::scale(p.dir, t[0]), b)),
                    vectormath::dot(cylinder.direction, vectormath::subtract(vectormath::scale(p.dir, t[1]), b))
                };
                if (-cylinder.half_height <= y[0] && y[0] <= cylinder.half_height && -cylinder.half_height <= y[1] && y[1] <= cylinder.half_height) {
                    // we are inside cylinder
                    return t;
                } else if (-cylinder.half_height <= y[0] && y[0] <= cylinder.half_height) {
                    // cap upper
                    if (std::abs(vectormath::dot(cylinder.direction, p.dir)) > GEOMETRIC_ERROR<double>()) {
                        const auto c = vectormath::add(b, vectormath::scale(cylinder.direction, cylinder.half_height));
                        t[1] = vectormath::dot(cylinder.direction, c) / vectormath::dot(cylinder.direction, p.dir);
                        return t;
                    }
                } else if (-cylinder.half_height <= y[1] && y[1] <= cylinder.half_height) {
                    if (std::abs(vectormath::dot(cylinder.direction, p.dir)) > GEOMETRIC_ERROR<double>()) {
                        const auto c = vectormath::add(b, vectormath::scale(cylinder.direction, -cylinder.half_height));
                        t[0] = vectormath::dot(cylinder.direction, c) / vectormath::dot(cylinder.direction, p.dir);
                        return t;
                    }
                }
                return std::nullopt;
            }
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