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
#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <vector>

#include <memory>

namespace dxmc {

template <Floating T, int FORWARD>
std::optional<T> intersectTriangle(const Particle<T>& particle, const std::array<T, 3>& v1, const std::array<T, 3>& v2, const std::array<T, 3>& v3)
{
    constexpr T epsilon = std::numeric_limits<T>::epsilon();
    const std::array<T, 3> v1v2 { v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2] };
    const std::array<T, 3> v1v3 { v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2] };
    const auto& pvec = vectormath::cross(particle.dir, v1v3);
    const auto det = vectormath::dot(v1v2, pvec);

    // Negative det gives intersection on backside of face
    if constexpr (FORWARD == 1) {
        if (det < epsilon)
            return std::nullopt;
    } else if constexpr (FORWARD == -1) {
        if (det > -epsilon)
            return std::nullopt;
    } else {
        if (std::abs(det) < epsilon)
            return std::nullopt;
    }
    const T invDet = T { 1 } / det;

    const std::array<T, 3> tvec { particle.pos[0] - v1[0], particle.pos[1] - v1[1], particle.pos[2] - v1[2] };
    const T u = vectormath::dot(tvec, pvec) * invDet;
    if (u < T { 0 } || u > T { 1 })
        return std::nullopt;

    const auto qvec = vectormath::cross(tvec, v1v2);
    const T v = vectormath::dot(particle.dir, qvec) * invDet;
    if (v < T { 0 } || u + v > T { 1 })
        return std::nullopt;

    return std::optional<T> { vectormath::dot(v1v3, qvec) * invDet };
}

template <Floating T>
std::optional<std::array<T, 2>> intersectAABB(const Particle<T>& p, const std::array<T, 6>& aabb)
{

    std::array<T, 2> t {
        std::numeric_limits<T>::min(),
        std::numeric_limits<T>::max()
    };
    for (std::size_t i = 0; i < 3; i++) {
        if (std::abs(p.dir[i]) > std::numeric_limits<T>::epsilon()) {
            const auto d_inv = T { 1 } / p.dir[i];
            const auto t1 = (aabb[i] - p.pos[i]) * d_inv;
            const auto t2 = (aabb[i + 3] - p.pos[i]) * d_inv;
            const auto t_min_cand = std::min(t1, t2);
            const auto t_max_cand = std::max(t1, t2);
            t[0] = std::max(t[0], t_min_cand);
            t[1] = std::min(t[1], t_max_cand);
        }
    }
    return t[0] > t[1] ? std::nullopt : t;
}

template <Floating T>
class KDTree {
public:
    KDTree(const std::vector<std::array<T, 3>>& vertices, std::vector<std::array<std::size_t, 3>> faces, const unsigned int k = 0)
    {
        m_D = k % 3;
        // sort by center of each triangle
        std::sort(faces.begin(), faces.end(), [&](const auto& lf, const auto& rh) {
            T lf_val = 0;
            T rh_val = 0;
            for (std::size_t i = 0; i < 3; i++) {
                lf_val += vertices[lf[i]][m_D];
                rh_val += vertices[rh[i]][m_D];
            }
            return lf_val < rh_val;
        });
        const auto plane_ind = faces.size() / 2;
        std::array<T, 3> plane_cand {
            vertices[faces[plane_ind][0]][m_D],
            vertices[faces[plane_ind][1]][m_D],
            vertices[faces[plane_ind][2]][m_D]
        };
        std::array<int, 3> fom_cand {
            figureOfMerit(vertices, faces, plane_cand[0]),
            figureOfMerit(vertices, faces, plane_cand[1]),
            figureOfMerit(vertices, faces, plane_cand[2])
        };

        const auto fom = fom_cand[argmin3(fom_cand)];
        if (fom == faces.size()) {
            m_faceIdx = faces;
        } else {
            m_plane = plane_cand[argmin3(plane_cand)];
            std::vector<std::array<std::size_t, 3>> left;
            std::vector<std::array<std::size_t, 3>> right;
            for (const auto& face : faces) {
                const auto side = planeSide(vertices, face, m_plane);
                if (side <= 0)
                    left.push_back(face);
                if (side >= 0)
                    right.push_back(face);
            }
            m_children = std::make_optional<std::array<KDTree<T>, 2>>({ KDTree<T>(vertices, left, k + 1),
                KDTree<T>(vertices, right, k + 1) });
        }
    }

    std::optional<T> intersect(const Particle<T>& particle, const std::vector<std::array<T, 3>>& vertices, const std::array<T, 6>& aabb) const
    {
        const auto& inter = intersectAABB(particle, vertices, aabb);
        return inter ? intersect(particle, vertices, *inter) : std::nullopt;
    }
    std::optional<T> intersect(const Particle<T>& particle, const std::vector<std::array<T, 3>>& vertices, const std::array<T, 2>& tbox) const
    {
        if (!m_children) {
            // intersect triangles between tbox and return;
            T t = std::numeric_limits<T>::max();
            for (const auto& face : m_faceIdx) {

                const auto& v1 = vertices[face[0]];
                const auto& v2 = vertices[face[0]];
                const auto& v3 = vertices[face[0]];
                const auto t_cand = intersectTriangle(particle, v1, v2, v3);
                t = t_cand ? std::min(t, t_cand) : t;
            }
            return t == std::numeric_limits<T>::max() ? std::nullopt : t;
        }

        // må ha test for parallell
        const T t = (particle.pos[m_D] - m_plane) / particle.dir[m_D];

        if (t > tbox[1]) {
            // left only
            if (std::signbit(particle.dir[m_D])) {
                m_children.value()[1].intersect(particle, tbox);
            } else {
                m_children.value()[0].intersect(particle, tbox);
            }

        } else if (t > tbox[0]) {
            // right only
            if (std::signbit(particle.dir[m_D])) {
                m_children.value()[0].intersect(particle, tbox);
            } else {
                m_children.value()[1].intersect(particle, tbox);
            }
        } else {
            // both directions (start with left)
            if (std::signbit(particle.dir[m_D])) {
                const std::array<T, 2> t_first { t, tbox[1] };
                const auto hit = m_children.value()[1].intersect(particle, t_first);
                if (hit) {
                    if (*hit <= t) {
                        return hit;
                    }
                }
                const std::array<T, 2> t_sec { tbox[0], t };
                return m_children.value()[0].intersect(particle, t_sec);

            } else {
                const std::array<T, 2> t_first { t, tbox[1] };
                const auto hit = m_children.value()[0].intersect(particle, t_first);
                if (hit) {
                    if (*hit <= t) {
                        return hit;
                    }
                }
                const std::array<T, 2> t_sec { tbox[0], t };
                return m_children.value()[1].intersect(particle, t_sec);
            }
        }

        /* if (std::signbit(dir[m_D])) {
            // swap left right
        }*/
    };

protected:
    static int argmin3(const std::array<int, 3>& v)
    {
        const auto it = std::min_element(v.cbegin(), v.cend());
        return it - v.cbegin();
    }
    int figureOfMerit(const std::vector<std::array<T, 3>>& vertices, const std::vector<std::array<std::size_t, 3>>& faceIdx, const T planesep) const
    {
        int fom = 0;
        int shared = 0;
        for (const auto& face : faceIdx) {
            const auto side = planeside(vertices, face, planesep);
            fom += side;
            if (side == 0)
                shared++;
        }
        return std::abs(fom) + shared;
    }
    static int planeSide(const std::vector<std::array<T, 3>>& vertices, const std::array<std::size_t, 3>& face, const T plane)
    {
        T max = std::numeric_limits<T>::min();
        T min = std::numeric_limits<T>::max();

        for (const std::size_t& ind : face) {
            min = std::min(vertices[ind][m_D], min);
            max = std::max(vertices[ind][m_D], max);
        }
        if (max <= plane)
            return -1;
        if (min >= plane)
            return 1;
        return 0;
    }

private:
    unsigned int m_D = 0;
    T m_plane = 0;
    std::vector<std::array<std::size_t, 3>> m_faceIdx;

    KDTree<T>* m_child;
    std::optional<std::array<KDTree<T>, 2>> m_children;
    
};

}