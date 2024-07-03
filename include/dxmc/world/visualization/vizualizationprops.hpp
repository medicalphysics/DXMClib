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
            const auto v1v2 = vectormath::dot(p.dir, m_dir);
            const auto p1v1 = vectormath::dot(p.pos, p.dir);
            const auto p2v1 = vectormath::dot(m_pos, p.dir);
            const auto p1v2 = vectormath::dot(p.pos, m_dir);
            const auto p2v2 = vectormath::dot(m_pos, m_dir);

            const auto t1 = ((p2v2 + p1v2) * v1v2 - p2v1 + p1v1) / (v1v2 * v1v2 - 1);

            if (t1 < 0.0)
                return std::nullopt;

            const auto t2 = v1v2 * t1 - p2v2 - p1v2;

            if ((t2 < 0.0) || (t2 > m_length))
                return std::nullopt;

            const auto c1 = vectormath::add(p.pos, vectormath::scale(p.dir, t1));
            const auto c2 = vectormath::add(m_pos, vectormath::scale(m_dir, t2));
            const auto d = vectormath::subtract(c1, c2);

            const auto l = vectormath::length(d);
            if (l < m_radii)
                return std::nullopt;
            const auto n = vectormath::scale(d, 1 / l);
            return std::make_optional(std::make_pair(t1, n));
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