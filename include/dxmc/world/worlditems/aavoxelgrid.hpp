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

#include "dxmc/dxmcrandom.hpp"
#include "dxmc/floating.hpp"
#include "dxmc/interactions.hpp"
#include "dxmc/material/material.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/basicshapes/aabb.hpp"
#include "dxmc/world/worldintersectionresult.hpp"
#include "dxmc/world/worlditems/worlditembase.hpp"

#include <array>

namespace dxmc {

template <Floating T, int NMaterialShells = 5>
class AAVoxelGrid final : public WorldItemBase<T> {
public:
    bool setData(const std::array<std::size_t, 3>& dim, const std::vector<T>& density, const std::vector<uint8_t>& materialIdx, const std::vector<Material2<T, NMaterialShells>>& materials)
    {
        const auto size = std::reduce(dim.cbegin(), dim.cend(), std::size_t { 1 }, std::multiplies<>());
        if (density.size() != size || materialIdx.size() != size) {
            return false;
        }
        // finding max material idx;
        auto maxMatIdx = std::max_element(std::execution::par_unseq, materialIdx.cbegin(), materialIdx.cend());
        if (*maxMatIdx >= materials.size()) {
            return false;
        }
        m_dim = dim;
        m_data.resize(size);
        DoseScore<T> dummy_dose;
        std::transform(std::execution::par_unseq, density.cbegin(), density.cend(), materialIdx.cbegin(), m_data.begin(), [=](const auto d, const auto mIdx) -> DataElement {
            return { .dose = dummy_dose, .density = d, .materialIndex = mIdx };
        });

        updateAABB();
        return true;
    }
    void setSpacing(const std::array<T, 3>& spacing)
    {
        for (std::size_t i = 0; i < 3; ++i) {
            m_invSpacing[i] = 1 / std::abs(spacing[i]);
        }
        updateAABB();
    }
    std::array<T, 3> spacing() const
    {
        std::array<T, 3> s;
        for (std::size_t i = 0; i < 3; ++i) {
            s[i] = 1 / m_invSpacing[i];
        }
        return s;
    }
    void translate(const std::array<T, 3>& dist) final
    {
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] += dist[i];
            m_aabb[i + 3] += dist[i];
        }
    }

    std::array<T, 3> center() const final
    {
        std::array<T, 3> center;
        for (std::size_t i = 0; i < 3; ++i) {
            center[i] = (m_aabb[i] + m_aabb[i + 3]) / 2;
        }
        return center;
    }

    std::array<T, 6> AABB() const final
    {
        return m_aabb;
    }

    inline std::size_t flatIndex(const std::array<std::size_t, 3>& index) const
    {
        return flatIndex(index[0], index[1], index[2]);
    }

    inline std::size_t flatIndex(const std::size_t x, const std::size_t y, const std::size_t z) const
    {
        return x + m_dim[0] * (y + m_dim[1] * z);
    }

    template <bool BOUNDSCHECK = true>
    inline std::size_t flatIndex(const std::array<T, 3>& pos) const
    {
        return flatIndex(index<BOUNDSCHECK>(pos));
    }

    template <bool BOUNDSCHECK = true>
    inline std::array<std::size_t, 3> index(const std::array<T, 3>& pos) const
    {
        if constexpr (BOUNDSCHECK) {
            std::array<std::size_t, 3> idx;
            for (std::size_t i = 0; i < 3; ++i) {
                const T bpos = std::clamp(pos, m_aabb[i], m_aabb[i + 3]);
                idx[i] = static_cast<std::size_t>(bpos * m_invSpacing[i]);
            }
            return idx;
        } else {
            const std::array<std::size_t, 3> idx = {
                static_cast<std::size_t>((pos[0] - m_aabb[0]) * m_invSpacing[0]),
                static_cast<std::size_t>((pos[1] - m_aabb[1]) * m_invSpacing[1]),
                static_cast<std::size_t>((pos[2] - m_aabb[2]) * m_invSpacing[2])
            };
            return idx;
        }
    }

    inline std::array<std::size_t, 3> index(const std::size_t flatIndex) const
    {
        const auto z = flatIndex / (m_dim[0] * m_dim[1]);
        const auto y = (flatIndex - z * m_dim[0] * m_dim[1]) / m_dim[0];
        const auto x = flatIndex - z * m_dim[0] * m_dim[1] - y * m_dim[0];
        std::array<std::size_t, 3> arr = { x, y, z };
        return arr;
    }

    WorldIntersectionResult<T> intersect(const Particle<T>& p) const final
    {
        auto intersect = basicshape::AABB::intersect(p, m_aabb);
        if (intersect.valid()) {
            voxelIntersect(p, intersect);
        }
        return intersect;
    }

    const DoseScore<T>& dose(std::size_t flatIndex = 0) const final
    {
        return m_data.at(flatIndex).dose;
    }

    void clearDose() final
    {
        std::for_each(std::execution::par_unseq, m_data.begin(), m_data.end(), [](auto& d) { d.dose.clear(); });
    }

    void transport(Particle<T>& p, RandomState& state) final { }

protected:
    void updateAABB()
    {
        const auto c = center();
        for (std::size_t i = 0; i < 3; ++i) {
            const T half_dist = (m_dim[i] * T { 0.5 }) / m_invSpacing[i];
            m_aabb[i] = -half_dist;
            m_aabb[i + 3] = half_dist;
        }
        translate(c);
    }

    static inline std::uint_fast8_t argmin3(const std::array<T, 3>& a)
    {
        return a[0] < a[1] ? a[0] < a[2] ? 0 : 2 : a[1] < a[2] ? 1
                                                               : 2;
    }

    void voxelIntersect(const Particle<T>& p, WorldIntersectionResult<T>& intersection) const
    {
        constexpr std::uint8_t ignoreIdx = 0;
        std::array<T, 3> pos;
        if (intersection.rayOriginIsInsideItem) {
            pos = p.pos;
        } else {
            for (std::size_t i = 0; i < 3; ++i) {
                const auto dir = p.dir[i] < 0 ? std::numeric_limits<T>::lowest() : std::numeric_limits<T>::max();
                if constexpr (sizeof(T) == 4)
                    pos[i] = std::nextafter(p.pos[i] + (intersection.intersection + 1E-4) * p.dir[i], dir);
                else
                    pos[i] = std::nextafter(p.pos[i] + (intersection.intersection + 1E-6) * p.dir[i], dir);
            }
        }
        auto xyz = index<false>(pos);
        auto index_flat = flatIndex(xyz);
        // first voxel is valid
        if (m_data[index_flat].materialIndex != ignoreIdx) {
            return;
        }

        const std::array<int, 3> xyz_step = {
            p.dir[0] < 0 ? -1 : 1,
            p.dir[1] < 0 ? -1 : 1,
            p.dir[2] < 0 ? -1 : 1
        };

        const std::array<int, 3> xyz_step_flat = {
            p.dir[0] < 0 ? -1 : 1,
            p.dir[1] < 0 ? -static_cast<int>(m_dim[0]) : static_cast<int>(m_dim[0]),
            p.dir[2] < 0 ? -static_cast<int>(m_dim[0] * m_dim[1]) : static_cast<int>(m_dim[0] * m_dim[1])
        };

        const std::array<T, 3> delta = {
            1 / (std::abs(p.dir[0]) * m_invSpacing[0]),
            1 / (std::abs(p.dir[1]) * m_invSpacing[1]),
            1 / (std::abs(p.dir[2]) * m_invSpacing[2])
        };
        std::array<T, 3> tMax = delta;

        auto dIdx = argmin3(tMax);
        xyz[dIdx] += xyz_step[dIdx];
        bool still_inside = xyz[dIdx] < m_dim[dIdx];
        index_flat += xyz_step_flat[dIdx];
        bool ignore_hit = still_inside ? m_data[index_flat].materialIndex == ignoreIdx : false;
        while (still_inside && ignore_hit) {
            dIdx = argmin3(tMax);
            xyz[dIdx] += xyz_step[dIdx];
            still_inside = xyz[dIdx] < m_dim[dIdx];
            index_flat += xyz_step_flat[dIdx];
            tMax[dIdx] += delta[dIdx];

            if (still_inside)
                ignore_hit = m_data[index_flat].materialIndex == ignoreIdx;
        }

        if (still_inside) {
            intersection.intersection += tMax[dIdx];
            intersection.rayOriginIsInsideItem = false;
        } else {
            intersection.intersectionValid = false;
        }
    }

    struct DataElement {
        DoseScore<T> dose;
        T density = 0;
        std::uint8_t materialIndex = 0;
    };

private:
    std::array<std::size_t, 3> m_dim = { 1, 1, 1 };
    std::array<T, 3> m_invSpacing = { 1, 1, 1 };
    std::array<T, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };
    std::vector<DataElement> m_data;
    std::vector<Material2<T, NMaterialShells>> m_materials;
};
}