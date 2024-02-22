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
#include "dxmc/world/basicshapes/aabb.hpp"

#include <array>
#include <limits>
#include <optional>

namespace dxmc {
namespace visualizationprops {
    class VisualizationProp {
    public:
        virtual std::optional<std::pair<double, std::array<double, 3>>> intersect(const Particle& p) const noexcept = 0;
    };

    class Line : public VisualizationProp {
    public:
        Line(const std::array<double, 3> start, const std::array<double, 3> dir, double length = -1, double radii = 1)
            : VisualizationProp()
            , m_pos(start)
            , m_dir(dir)
            , m_radii(std::abs(radii))
        {
            vectormath::normalize(m_dir);
            if (length <= 0)
                m_length = std::numeric_limits<double>::max();
            else
                m_length = length;
        }

        std::optional<std::pair<double, std::array<double, 3>>> intersect(const Particle& p) const noexcept override
        {
            const auto aabb = AABB();
            const auto inter = basicshape::AABB::intersect(p, aabb);
            if (!inter.valid())
                return std::nullopt;

            const auto n = vectormath::cross(p.dir, m_dir);
            const auto n_length = vectormath::length(n);
            if (std::abs(n_length) < GEOMETRIC_ERROR()) { // lines are parallell
                return std::nullopt;
            }
            const auto n_invlength = 1 / n_length;
            const auto distance = std::abs(vectormath::dot(n, vectormath::subtract(m_pos, p.pos)) * n_invlength);

            if (distance > m_radii) { // intersection is larger than radii
                return std::nullopt;
            }

            // shiftling line aling normal to both lines to find intersection point at radi
            const auto s_pos = vectormath::add(m_pos, vectormath::scale(n, m_radii * n_invlength));

            const auto n2 = vectormath::cross(m_dir, n);
            const auto t1 = vectormath::dot(vectormath::subtract(s_pos, p.pos), n2) / vectormath::dot(n2, p.dir);

            //  intersection pos
            const auto p_int = vectormath::add(p.pos, vectormath::scale(p.dir, t1));
            // distance from intpos to start
            const auto int_distance = vectormath::dot(vectormath::subtract(p_int, m_pos), m_dir);
            if (int_distance < 0 || int_distance > m_length)
                return std::nullopt;

            auto d = vectormath::subtract(p_int, p.pos);
            d = vectormath::subtract(d, vectormath::scale(m_dir, vectormath::dot(d, m_dir)));
            vectormath::normalize(d);
            return std::make_optional(std::make_pair(t1, d));
        }

        std::array<double, 6> AABB() const
        {
            const auto end = vectormath::add(m_pos, vectormath::scale(m_dir, m_length));
            std::array<double, 6> aabb;
            for (std::size_t i = 0; i < 3; ++i) {
                aabb[i] = std::min(m_pos[i], end[i]) - GEOMETRIC_ERROR();
                aabb[i + 3] = std::max(m_pos[i], end[i]) + GEOMETRIC_ERROR();
            }
            return aabb;
        }

    private:
        std::array<double, 3> m_pos = { 0, 0, 0 };
        std::array<double, 3> m_dir = { 0, 0, 1 };
        double m_radii = 1;
        double m_length = 1;
    };

}

}