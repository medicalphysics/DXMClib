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
#include "dxmc/world/worlditems/worlditembase.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

namespace dxmc {

template <Floating T>
class FluenceScore final : public WorldItemBase<T> {
public:
    FluenceScore(T radius = T { 16 }, const std::array<T, 3>& center = { 0, 0, 0 }, const std::array<T, 3>& normal = { 0, 0, 1 })
        : WorldItemBase<T>()
        , m_radius(radius)
        , m_center(center)
        , m_normal(normal)
    {
        vectormath::normalize(m_normal);
        setEnergyStep(1);
    }

    void setEnergyStep(T step)
    {
        m_energy_step = std::max(step, T { 0.1 });
        const auto N = static_cast<std::uint64_t>(MAX_ENERGY<T>() / m_energy_step) + 1;
        m_intensity.resize(N);
        std::fill(m_intensity.begin(), m_intensity.end(), T { 0 });
    }

    void translate(const std::array<T, 3>& dist)
    {
        m_center = vectormath::add(m_center, dist);
    }
    std::array<T, 3> center() const { return m_center; }
    std::array<T, 6> AABB() const
    {
        const std::array<std::array<T, 3>, 3> span = {
            vectormath::scale(vectormath::cross(m_normal, { 1, 0, 0 }), m_radius),
            vectormath::scale(vectormath::cross(m_normal, { 0, 1, 0 }), m_radius),
            vectormath::scale(vectormath::cross(m_normal, { 0, 0, 1 }), m_radius)
        };
        std::array<T, 6> aabb = {
            std::numeric_limits<T>::max(),
            std::numeric_limits<T>::max(),
            std::numeric_limits<T>::max(),
            std::numeric_limits<T>::lowest(),
            std::numeric_limits<T>::lowest(),
            std::numeric_limits<T>::lowest()
        };

        for (const auto& s : span)
            for (std::size_t i = 0; i < 3; ++i) {
                const auto [mi, ma] = std::minmax(m_center[i] - p[i], m_center[i] + p[i]);
                aabb[i] = std::min(aabb[i], mi);
                aabb[i + 3] = std::max(aabb[i + 3], ma);
            }
        return aabb;
    }
    virtual WorldIntersectionResult<T> intersect(const Particle<T>& p) const = 0;
    virtual VisualizationIntersectionResult<T, WorldItemBase<T>> intersectVisualization(const Particle<T>& p) const = 0;
    virtual const DoseScore<T>& dose(std::size_t index = 0) const = 0;
    virtual void clearDose() = 0;
    virtual void transport(Particle<T>& p, RandomState& state) = 0;

private:
    std::array<T, 3> m_center = { 0, 0, 0 };
    std::array<T, 3> m_normal = { 0, 0, 1 };
    T m_radius = 16;
    T m_energy_step = 1;
    std::vector<std::uint64_t> m_intensity;
    DoseScore<T> m_dummyDose;
}
}