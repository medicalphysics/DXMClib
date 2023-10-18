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
#include "dxmc/world/basicshapes/aabb.hpp"
#include "dxmc/world/worlditems/worlditembase.hpp"

#include <algorithm>
#include <array>
#include <atomic>
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
        , m_center(center)
        , m_radius(radius)
    {
        setPlaneNormal(normal);
        setEnergyStep(1);
    }

    void setEnergyStep(T step)
    {
        m_energy_step = std::max(step, T { 0.1 });
        const auto N = static_cast<std::uint64_t>(MAX_ENERGY<T>() / m_energy_step) + 1;
        m_intensity.resize(N);
        std::fill(m_intensity.begin(), m_intensity.end(), 0);
    }

    void translate(const std::array<T, 3>& dist) override
    {
        for (std::size_t i = 0; i < 3; ++i) {
            m_center[i] += dist[i];
            m_aabb[i] += dist[i];
            m_aabb[i + 3] += dist[i];
        }
    }

    std::array<T, 3> center() const override
    {
        return m_center;
    }
    std::array<T, 6> AABB() const override
    {
        return m_aabb;
    }

    void setPlaneNormal(const std::array<T, 3>& normal)
    {
        m_normal = normal;
        vectormath::normalize(m_normal);
        calculateAABB();
    }

    std::vector<std::pair<T, std::uint64_t>> getSpecter() const
    {
        std::vector<std::pair<T, std::uint64_t>> spec(m_intensity.size());
        for (std::size_t i = 0; i < m_intensity.size(); ++i) {
            spec[i] = std::make_pair(m_energy_step * i, m_intensity[i]);
        }
        return spec;
    }

    WorldIntersectionResult<T> intersect(const Particle<T>& p) const override
    {
        const auto aabb_inter = basicshape::AABB::intersectForwardInterval(p, m_aabb);
        return aabb_inter ? intersectDisc(p) : WorldIntersectionResult<T> {};
    }

    VisualizationIntersectionResult<T, WorldItemBase<T>> intersectVisualization(const Particle<T>& p) const override
    {
        const auto res = intersect(p);
        VisualizationIntersectionResult<T, WorldItemBase<T>> w;
        if (res.valid()) {
            w.rayOriginIsInsideItem = false;
            w.intersection = res.intersection;
            w.intersectionValid = true;
            w.item = this;
            w.normal = vectormath::dot(p.dir, m_normal) <= T { 0 } ? m_normal : vectormath::scale(m_normal, T { -1 });
        }
        return w;
    }

    const EnergyScore<T>& energyScored(std::size_t index = 0) const override
    {
        return m_energyScored;
    }

    void clearEnergyScored() override
    {
        m_energyScored.clear();
        std::fill(m_intensity.begin(), m_intensity.end(), std::uint64_t { 0 });
    }

    void addEnergyScoredToDoseScore(T calibration_factor = 1) override
    {
        // not defined for fleunce counter
        return;
    }

    const DoseScore<T>& doseScored(std::size_t index = 0) const override
    {
        return m_dummyDose;
    }

    void clearDoseScored() override
    {
        return;
    }

    void transport(Particle<T>& particle, RandomState& state) override
    {
        // Assuming particle is on the disc
        const auto eIdx = static_cast<std::size_t>(particle.energy / m_energy_step);
        auto counter = std::atomic_ref(m_intensity[eIdx]);
        counter++;
        m_energyScored.scoreEnergy(particle.energy);
        particle.border_translate(T { 0 });
    }

protected:
    static inline std::pair<T, T> minmax(T v1, T v2)
    {
        // use own minmax instead of std::minamx due to bug or weird feature of MSVC compiler with /O2
        return v1 <= v2 ? std::make_pair(v1, v2) : std::make_pair(v2, v1);
    }

    void calculateAABB()
    {
        const std::array<std::array<T, 3>, 3> span = {
            vectormath::scale(vectormath::cross(m_normal, { 1, 0, 0 }), m_radius),
            vectormath::scale(vectormath::cross(m_normal, { 0, 1, 0 }), m_radius),
            vectormath::scale(vectormath::cross(m_normal, { 0, 0, 1 }), m_radius)
        };
        m_aabb = {
            std::numeric_limits<T>::max(),
            std::numeric_limits<T>::max(),
            std::numeric_limits<T>::max(),
            std::numeric_limits<T>::lowest(),
            std::numeric_limits<T>::lowest(),
            std::numeric_limits<T>::lowest()
        };

        for (const auto& s : span)
            for (std::size_t i = 0; i < 3; ++i) {
                const auto [mi, ma] = minmax(m_center[i] - s[i], m_center[i] + s[i]);
                m_aabb[i] = std::min(m_aabb[i], mi);
                m_aabb[i + 3] = std::max(m_aabb[i + 3], ma);
            }
        // ensure a min size;
        for (std::size_t i = 0; i < 3; ++i) {
            if (m_aabb[i + 3] - m_aabb[i] < 2 * GEOMETRIC_ERROR<T>()) {
                m_aabb[i] -= GEOMETRIC_ERROR<T>();
                m_aabb[i + 3] += GEOMETRIC_ERROR<T>();
            }
        }
    }

    WorldIntersectionResult<T> intersectDisc(const Particle<T>& p) const
    {
        WorldIntersectionResult<T> res;
        const auto D = vectormath::dot(p.dir, m_normal);
        constexpr T minOrt = 1E-6;
        if (D < minOrt && D > -minOrt)
            return res; // dir and normal is orthogonal, we exits
        res.intersection = vectormath::dot(vectormath::subtract(m_center, p.pos), m_normal) / D;
        if (res.intersection <= T { 0 })
            return res;

        // intersection point
        const auto p_int = vectormath::add(p.pos, vectormath::scale(p.dir, res.intersection));

        // distance from center
        const auto c_dist = vectormath::subtract(m_center, p_int);
        // check if distance from center is less than radius
        if (vectormath::dot(c_dist, c_dist) <= m_radius * m_radius) {
            res.intersectionValid = true;
        }
        return res;
    }

private:
    std::array<T, 3> m_center = { 0, 0, 0 };
    std::array<T, 3> m_normal = { 0, 0, 1 };
    T m_radius = 16;
    T m_energy_step = 1;
    std::array<T, 6> m_aabb;
    std::vector<std::uint64_t> m_intensity;
    EnergyScore<T> m_energyScored;
    DoseScore<T> m_dummyDose;
};
}