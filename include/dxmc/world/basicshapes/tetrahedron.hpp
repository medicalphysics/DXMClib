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

        template <Floating T>
        inline T permuted_inner_product(const std::array<std::array<T, 3>, 2>& lh, const std::array<std::array<T, 3>, 2>& rh)
        {
            // triple product a*(b x c)
            return vectormath::dot(lh[0], rh[1]) + vectormath::dot(lh[1], rh[0]);
        }

        template <Floating T>
        std::optional<std::array<T, 2>> intersectForwardInterval(
            const Particle<T>& particle,
            const std::array<T, 3>& v0,
            const std::array<T, 3>& v1,
            const std::array<T, 3>& v2,
            const std::array<T, 3>& v3)
        {
            // intersection using plucker coordinates

            const std::array<std::array<T, 3>, 2> p = { particle.dir, vectormath::cross(particle.dir, particle.pos) };

            // Faces: F3=(V0 V1 V2), F2=(V1 V0 V3), F1=(V2 V3 V0), F0=(V3 V2 V1)
            // Edges: E30=(V1 V2), E31=(V2 V0), E32=(V0 V1)
            // Edges: E20=(V0 V3), E21=(V3 V1), E22=(V1 V0)
            // Edges: E10=(V3 V0), E11=(V0 V2), E12=(V2 V3)
            // Edges: E00=(V2 V1), E01=(V1 V3), E02=(V3 V2)

            // Edges
            const auto e30 = vectormath::subtract(v2, v1);
            const auto e31 = vectormath::subtract(v0, v2);
            const auto e32 = vectormath::subtract(v1, v0);
            const auto e20 = vectormath::subtract(v3, v0);
            const auto e21 = vectormath::subtract(v1, v3);
            const auto e22 = vectormath::subtract(v0, v1);
            const auto e10 = vectormath::subtract(v0, v3);
            const auto e11 = vectormath::subtract(v2, v0);
            const auto e12 = vectormath::subtract(v3, v2);
            const auto e00 = vectormath::subtract(v1, v2);
            const auto e01 = vectormath::subtract(v3, v1);
            const auto e02 = vectormath::subtract(v2, v3);

            // Plucker koordinates for edges
            const std::array<std::array<T, 3>, 2> pe30 = { e30, vectormath::cross(e30, v1) };
            const std::array<std::array<T, 3>, 2> pe31 = { e31, vectormath::cross(e31, v2) };
            const std::array<std::array<T, 3>, 2> pe32 = { e32, vectormath::cross(e32, v0) };
            const std::array<std::array<T, 3>, 2> pe20 = { e20, vectormath::cross(e20, v0) };
            const std::array<std::array<T, 3>, 2> pe21 = { e21, vectormath::cross(e21, v3) };
            const std::array<std::array<T, 3>, 2> pe22 = { e22, vectormath::cross(e22, v1) };
            const std::array<std::array<T, 3>, 2> pe10 = { e10, vectormath::cross(e10, v3) };
            const std::array<std::array<T, 3>, 2> pe11 = { e11, vectormath::cross(e11, v0) };
            const std::array<std::array<T, 3>, 2> pe12 = { e12, vectormath::cross(e12, v2) };
            const std::array<std::array<T, 3>, 2> pe00 = { e00, vectormath::cross(e00, v2) };
            const std::array<std::array<T, 3>, 2> pe01 = { e01, vectormath::cross(e01, v1) };
            const std::array<std::array<T, 3>, 2> pe02 = { e02, vectormath::cross(e02, v3) };

            // Inner product
            const auto r30 = permuted_inner_product(p, pe30);
            const auto r31 = permuted_inner_product(p, pe31);
            const auto r32 = permuted_inner_product(p, pe32);
            const auto r20 = permuted_inner_product(p, pe20);
            const auto r21 = permuted_inner_product(p, pe21);
            const auto r22 = permuted_inner_product(p, pe22);
            const auto r10 = permuted_inner_product(p, pe10);
            const auto r11 = permuted_inner_product(p, pe11);
            const auto r12 = permuted_inner_product(p, pe12);
            const auto r00 = permuted_inner_product(p, pe00);
            const auto r01 = permuted_inner_product(p, pe01);
            const auto r02 = permuted_inner_product(p, pe02);

            bool search_enter = true;
            bool search_exit = true;

            std::array<T, 2> t;

            constexpr auto e = GEOMETRIC_ERROR<T>();

            // testing Face 3
            if (std::abs(r30) < e || std::abs(r31) < e || std::abs(r32) < e) {
                if (r30 >= 0 && r31 >= 0 && r32 >= 0) {
                    const auto wsum_inv = 1 / (r30 + r31 + r32);
                    const auto vIdx = vectormath::argmax3<std::uint_fast32_t>(particle.dir);
                    t[0] = (r30 * wsum_inv * v1[vIdx] + r31 * wsum_inv * v2[vIdx] * r32 * wsum_inv * v0[vIdx] - particle.pos[vIdx]) / particle.dir[vIdx];
                    search_enter = false;
                }
                if (r30 <= 0 && r31 <= 0 && r32 <= 0) {
                    const auto wsum_inv = 1 / (r30 + r31 + r32);
                    const auto vIdx = vectormath::argmax3<std::uint_fast32_t>(particle.dir);
                    t[1] = (r30 * wsum_inv * v1[vIdx] + r31 * wsum_inv * v2[vIdx] * r32 * wsum_inv * v0[vIdx] - particle.pos[vIdx]) / particle.dir[vIdx];
                    search_exit = false;
                }
            }
            // testing Face 2
            if (std::abs(r20) < e || std::abs(r21) < e || std::abs(r22) < e) {
                if (search_enter && r20 >= 0 && r21 >= 0 && r22 >= 0) {
                    const auto wsum_inv = 1 / (r20 + r21 + r22);
                    const auto vIdx = vectormath::argmax3<std::uint_fast32_t>(particle.dir);
                    t[0] = (r20 * wsum_inv * v0[vIdx] + r21 * wsum_inv * v3[vIdx] * r22 * wsum_inv * v1[vIdx] - particle.pos[vIdx]) / particle.dir[vIdx];
                    search_enter = false;
                }
                if (search_exit && r20 <= 0 && r21 <= 0 && r22 <= 0) {
                    const auto wsum_inv = 1 / (r20 + r21 + r22);
                    const auto vIdx = vectormath::argmax3<std::uint_fast32_t>(particle.dir);
                    t[1] = (r20 * wsum_inv * v0[vIdx] + r21 * wsum_inv * v3[vIdx] * r22 * wsum_inv * v1[vIdx] - particle.pos[vIdx]) / particle.dir[vIdx];
                    search_exit = false;
                }
            }
            // testing Face 1
            if (std::abs(r10) < e || std::abs(r11) < e || std::abs(r12) < e) {
                if (search_enter && r10 >= 0 && r11 >= 0 && r12 >= 0) {
                    const auto wsum_inv = 1 / (r10 + r11 + r12);
                    const auto vIdx = vectormath::argmax3<std::uint_fast32_t>(particle.dir);
                    t[0] = (r10 * wsum_inv * v3[vIdx] + r11 * wsum_inv * v0[vIdx] * r12 * wsum_inv * v2[vIdx] - particle.pos[vIdx]) / particle.dir[vIdx];
                    search_enter = false;
                }
                if (search_exit && r10 <= 0 && r11 <= 0 && r12 <= 0) {
                    const auto wsum_inv = 1 / (r10 + r11 + r12);
                    const auto vIdx = vectormath::argmax3<std::uint_fast32_t>(particle.dir);
                    t[1] = (r10 * wsum_inv * v3[vIdx] + r11 * wsum_inv * v0[vIdx] * r12 * wsum_inv * v2[vIdx] - particle.pos[vIdx]) / particle.dir[vIdx];
                    search_exit = false;
                }
            }
            if (search_enter && search_exit) {
                return std::nullopt;
            } else {
                // testing Face 0
                if (std::abs(r00) < e || std::abs(r01) < e || std::abs(r02) < e) {
                    if (search_enter && r00 >= 0 && r01 >= 0 && r02 >= 0) {
                        const auto wsum_inv = 1 / (r00 + r01 + r02);
                        const auto vIdx = vectormath::argmax3<std::uint_fast32_t>(particle.dir);
                        t[0] = (r00 * wsum_inv * v2[vIdx] + r01 * wsum_inv * v1[vIdx] * r02 * wsum_inv * v2[vIdx] - particle.pos[vIdx]) / particle.dir[vIdx];
                        search_enter = false;
                    }
                    if (search_exit && r00 <= 0 && r01 <= 0 && r02 <= 0) {
                        const auto wsum_inv = 1 / (r00 + r01 + r02);
                        const auto vIdx = vectormath::argmax3<std::uint_fast32_t>(particle.dir);
                        t[1] = (r00 * wsum_inv * v2[vIdx] + r01 * wsum_inv * v1[vIdx] * r02 * wsum_inv * v2[vIdx] - particle.pos[vIdx]) / particle.dir[vIdx];
                        search_exit = false;
                    }
                }
            }
            if (!search_enter && !search_exit) {
                // we have found both intersections
                return std::make_optional(t);
            }
            return std::nullopt;
        }

        template <Floating T>
        WorldIntersectionResult<T> intersect(
            const Particle<T>& p,
            const std::array<T, 3>& v0,
            const std::array<T, 3>& v1,
            const std::array<T, 3>& v2,
            const std::array<T, 3>& v3)
        {
            const auto t_opt = intersectForwardInterval(p, v0, v1, v2, v3);
            WorldIntersectionResult<T> res;
            if (t_opt) {
                res.intersectionValid = true;
                const auto& t = *t_opt;
                res.rayOriginIsInsideItem = t[0] < 0 && 0 < t[1];
                res.intersection = res.rayOriginIsInsideItem ? t[1] : t[0];
            }
            return res;
        }

        template <Floating T, typename U>
        VisualizationIntersectionResult<T, U> intersectVisualization(
            const Particle<T>& particle,
            const std::array<T, 3>& v0,
            const std::array<T, 3>& v1,
            const std::array<T, 3>& v2,
            const std::array<T, 3>& v3)
        {

            VisualizationIntersectionResult<T, U> res;
            // intersection using plucker coordinates

            const std::array<std::array<T, 3>, 2> p = { particle.dir, vectormath::cross(particle.dir, particle.pos) };

            // Faces: F3=(V0 V1 V2), F2=(V1 V0 V3), F1=(V2 V3 V0), F0=(V3 V2 V1)
            // Edges: E30=(V1 V2), E31=(V2 V0), E32=(V0 V1)
            // Edges: E20=(V0 V3), E21=(V3 V1), E22=(V1 V0)
            // Edges: E10=(V3 V0), E11=(V0 V2), E12=(V2 V3)
            // Edges: E00=(V2 V1), E01=(V1 V3), E02=(V3 V2)

            // Edges
            const auto e30 = vectormath::subtract(v2, v1);
            const auto e31 = vectormath::subtract(v0, v2);
            const auto e32 = vectormath::subtract(v1, v0);
            const auto e20 = vectormath::subtract(v3, v0);
            const auto e21 = vectormath::subtract(v1, v3);
            const auto e22 = vectormath::subtract(v0, v1);
            const auto e10 = vectormath::subtract(v0, v3);
            const auto e11 = vectormath::subtract(v2, v0);
            const auto e12 = vectormath::subtract(v3, v2);
            const auto e00 = vectormath::subtract(v1, v2);
            const auto e01 = vectormath::subtract(v3, v1);
            const auto e02 = vectormath::subtract(v2, v3);

            // Plucker koordinates for edges
            const std::array<std::array<T, 3>, 2> pe30 = { e30, vectormath::cross(e30, v1) };
            const std::array<std::array<T, 3>, 2> pe31 = { e31, vectormath::cross(e31, v2) };
            const std::array<std::array<T, 3>, 2> pe32 = { e32, vectormath::cross(e32, v0) };
            const std::array<std::array<T, 3>, 2> pe20 = { e20, vectormath::cross(e20, v0) };
            const std::array<std::array<T, 3>, 2> pe21 = { e21, vectormath::cross(e21, v3) };
            const std::array<std::array<T, 3>, 2> pe22 = { e22, vectormath::cross(e22, v1) };
            const std::array<std::array<T, 3>, 2> pe10 = { e10, vectormath::cross(e10, v3) };
            const std::array<std::array<T, 3>, 2> pe11 = { e11, vectormath::cross(e11, v0) };
            const std::array<std::array<T, 3>, 2> pe12 = { e12, vectormath::cross(e12, v2) };
            const std::array<std::array<T, 3>, 2> pe00 = { e00, vectormath::cross(e00, v2) };
            const std::array<std::array<T, 3>, 2> pe01 = { e01, vectormath::cross(e01, v1) };
            const std::array<std::array<T, 3>, 2> pe02 = { e02, vectormath::cross(e02, v3) };

            // Inner product
            const auto r30 = permuted_inner_product(p, pe30);
            const auto r31 = permuted_inner_product(p, pe31);
            const auto r32 = permuted_inner_product(p, pe32);
            const auto r20 = permuted_inner_product(p, pe20);
            const auto r21 = permuted_inner_product(p, pe21);
            const auto r22 = permuted_inner_product(p, pe22);
            const auto r10 = permuted_inner_product(p, pe10);
            const auto r11 = permuted_inner_product(p, pe11);
            const auto r12 = permuted_inner_product(p, pe12);
            const auto r00 = permuted_inner_product(p, pe00);
            const auto r01 = permuted_inner_product(p, pe01);
            const auto r02 = permuted_inner_product(p, pe02);

            bool search_enter = true;
            bool search_exit = true;
            std::array<T, 3> normal_enter, normal_exit;
            std::array<T, 2> t;

            constexpr auto e = GEOMETRIC_ERROR<T>();

            // testing Face 3
            if (std::abs(r30) < e || std::abs(r31) < e || std::abs(r32) < e) {
                if (r30 >= 0 && r31 >= 0 && r32 >= 0) {
                    const auto wsum_inv = 1 / (r30 + r31 + r32);
                    const auto vIdx = vectormath::argmax3<std::uint_fast32_t>(particle.dir);
                    t[0] = (r30 * wsum_inv * v1[vIdx] + r31 * wsum_inv * v2[vIdx] * r32 * wsum_inv * v0[vIdx] - particle.pos[vIdx]) / particle.dir[vIdx];
                    search_enter = false;
                    normal_enter = vectormath::cross(v1, v2);
                    vectormath::normalize(normal_enter);
                }
                if (r30 <= 0 && r31 <= 0 && r32 <= 0) {
                    const auto wsum_inv = 1 / (r30 + r31 + r32);
                    const auto vIdx = vectormath::argmax3<std::uint_fast32_t>(particle.dir);
                    t[1] = (r30 * wsum_inv * v1[vIdx] + r31 * wsum_inv * v2[vIdx] * r32 * wsum_inv * v0[vIdx] - particle.pos[vIdx]) / particle.dir[vIdx];
                    search_exit = false;
                    normal_exit = vectormath::cross(v1, v2);
                    vectormath::normalize(normal_exit);
                }
            }
            // testing Face 2
            if (std::abs(r20) < e || std::abs(r21) < e || std::abs(r22) < e) {
                if (search_enter && r20 >= 0 && r21 >= 0 && r22 >= 0) {
                    const auto wsum_inv = 1 / (r20 + r21 + r22);
                    const auto vIdx = vectormath::argmax3<std::uint_fast32_t>(particle.dir);
                    t[0] = (r20 * wsum_inv * v0[vIdx] + r21 * wsum_inv * v3[vIdx] * r22 * wsum_inv * v1[vIdx] - particle.pos[vIdx]) / particle.dir[vIdx];
                    search_enter = false;
                    normal_enter = vectormath::cross(v0, v3);
                    vectormath::normalize(normal_enter);
                }
                if (search_exit && r20 <= 0 && r21 <= 0 && r22 <= 0) {
                    const auto wsum_inv = 1 / (r20 + r21 + r22);
                    const auto vIdx = vectormath::argmax3<std::uint_fast32_t>(particle.dir);
                    t[1] = (r20 * wsum_inv * v0[vIdx] + r21 * wsum_inv * v3[vIdx] * r22 * wsum_inv * v1[vIdx] - particle.pos[vIdx]) / particle.dir[vIdx];
                    search_exit = false;
                    normal_exit = vectormath::cross(v0, v3);
                    vectormath::normalize(normal_exit);
                }
            }
            // testing Face 1
            if (std::abs(r10) < e || std::abs(r11) < e || std::abs(r12) < e) {
                if (search_enter && r10 >= 0 && r11 >= 0 && r12 >= 0) {
                    const auto wsum_inv = 1 / (r10 + r11 + r12);
                    const auto vIdx = vectormath::argmax3<std::uint_fast32_t>(particle.dir);
                    t[0] = (r10 * wsum_inv * v3[vIdx] + r11 * wsum_inv * v0[vIdx] * r12 * wsum_inv * v2[vIdx] - particle.pos[vIdx]) / particle.dir[vIdx];
                    search_enter = false;
                    normal_enter = vectormath::cross(v3, v0);
                    vectormath::normalize(normal_enter);
                }
                if (search_exit && r10 <= 0 && r11 <= 0 && r12 <= 0) {
                    const auto wsum_inv = 1 / (r10 + r11 + r12);
                    const auto vIdx = vectormath::argmax3<std::uint_fast32_t>(particle.dir);
                    t[1] = (r10 * wsum_inv * v3[vIdx] + r11 * wsum_inv * v0[vIdx] * r12 * wsum_inv * v2[vIdx] - particle.pos[vIdx]) / particle.dir[vIdx];
                    search_exit = false;
                    normal_exit = vectormath::cross(v3, v0);
                    vectormath::normalize(normal_exit);
                }
            }
            if (search_enter && search_exit) {
                return res;
            } else {
                // testing Face 0
                if (std::abs(r00) < e || std::abs(r01) < e || std::abs(r02) < e) {
                    if (search_enter && r00 >= 0 && r01 >= 0 && r02 >= 0) {
                        const auto wsum_inv = 1 / (r00 + r01 + r02);
                        const auto vIdx = vectormath::argmax3<std::uint_fast32_t>(particle.dir);
                        t[0] = (r00 * wsum_inv * v2[vIdx] + r01 * wsum_inv * v1[vIdx] * r02 * wsum_inv * v2[vIdx] - particle.pos[vIdx]) / particle.dir[vIdx];
                        search_enter = false;
                        normal_enter = vectormath::cross(v2, v1);
                        vectormath::normalize(normal_enter);
                    }
                    if (search_exit && r00 <= 0 && r01 <= 0 && r02 <= 0) {
                        const auto wsum_inv = 1 / (r00 + r01 + r02);
                        const auto vIdx = vectormath::argmax3<std::uint_fast32_t>(particle.dir);
                        t[1] = (r00 * wsum_inv * v2[vIdx] + r01 * wsum_inv * v1[vIdx] * r02 * wsum_inv * v2[vIdx] - particle.pos[vIdx]) / particle.dir[vIdx];
                        search_exit = false;
                        normal_exit = vectormath::cross(v2, v1);
                        vectormath::normalize(normal_exit);
                    }
                }
            }
            if (!search_enter && !search_exit) {
                res.intersectionValid = true;
                res.rayOriginIsInsideItem = t[0] < 0 && 0 < t[1];
                res.intersection = res.rayOriginIsInsideItem ? t[1] : t[0];
                res.normal = res.rayOriginIsInsideItem ? normal_exit : normal_enter;
            }
            return res;
        }
    }
}
}
