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

template <Floating T, std::size_t NMaterialShells = 5, int LOWENERGYCORRECTION = 2>
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
        m_materials = materials;

        updateAABB();
        return true;
    }
    void setSpacing(const std::array<T, 3>& spacing)
    {
        for (std::size_t i = 0; i < 3; ++i) {
            m_spacing[i] = std::abs(spacing[i]);
            m_invSpacing[i] = 1 / m_spacing[i];
        }
        updateAABB();
    }

    const std::array<T, 3>& spacing() const
    {
        return m_spacing;
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
                const T dt = m_spacing[i] * T { 0.5 };
                const T bpos = std::clamp(pos[i], m_aabb[i] + dt, m_aabb[i + 3] - dt);
                idx[i] = static_cast<std::size_t>((bpos - m_aabb[i]) * m_invSpacing[i]);
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

    void transport(Particle<T>& p, RandomState& state) final
    {
        constexpr std::uint_fast8_t ignoreMaterialIdx = 0;
        voxelTransport<ignoreMaterialIdx>(p, state);
    }

protected:
    void updateAABB()
    {
        const auto c = center();
        for (std::size_t i = 0; i < 3; ++i) {
            const T half_dist = (m_dim[i] * T { 0.5 }) * m_spacing[i];
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
    template <std::uint_fast8_t IGNOREIDX = 0>
    void voxelIntersect(const Particle<T>& p, WorldIntersectionResult<T>& intersection) const
    {
        std::array<std::size_t, 3> xyz;
        if (intersection.rayOriginIsInsideItem) {
            xyz = index<false>(p.pos);
        } else {
            // make sure we are well inside a voxel
            std::array<T, 3> pos;
            for (std::size_t i = 0; i < 3; ++i) {
                pos[i] = p.pos[i] + p.dir[i] * (intersection.intersection + m_spacing[i] * T { 0.5 });
            }
            xyz = index<false>(pos);
        }

        auto index_flat = flatIndex(xyz);
        if (m_data[index_flat].materialIndex != IGNOREIDX) {
            return;
        }

        const std::array<int, 3> xyz_step = {
            p.dir[0] < 0 ? -1 : 1,
            p.dir[1] < 0 ? -1 : 1,
            p.dir[2] < 0 ? -1 : 1
        };

        // by ising int, max x*y size is about 2e9
        const std::array<int, 3> xyz_step_flat = {
            p.dir[0] < 0 ? -1 : 1,
            p.dir[1] < 0 ? -static_cast<int>(m_dim[0]) : static_cast<int>(m_dim[0]),
            p.dir[2] < 0 ? -static_cast<int>(m_dim[0] * m_dim[1]) : static_cast<int>(m_dim[0] * m_dim[1])
        };

        const std::array<T, 3> delta = {
            m_spacing[0] / std::abs(p.dir[0]),
            m_spacing[1] / std::abs(p.dir[1]),
            m_spacing[2] / std::abs(p.dir[2])
        };

        std::array<T, 3> tMax = {
            // abs function in case of delta is infinity
            std::abs(std::nextafter(delta[0], std::numeric_limits<T>::max()) - delta[0]),
            std::abs(std::nextafter(delta[1], std::numeric_limits<T>::max()) - delta[1]),
            std::abs(std::nextafter(delta[2], std::numeric_limits<T>::max()) - delta[2])
        };

        bool still_inside = true;
        std::uint_fast8_t dIdx;
        while (still_inside && m_data[index_flat].materialIndex == IGNOREIDX) {
            dIdx = argmin3(tMax);
            xyz[dIdx] += xyz_step[dIdx];
            still_inside = xyz[dIdx] < m_dim[dIdx];
            index_flat += xyz_step_flat[dIdx];
            tMax[dIdx] += delta[dIdx];
        }

        if (still_inside) {
            intersection.intersection += tMax[dIdx];
            intersection.rayOriginIsInsideItem = false;
        } else {
            intersection.intersectionValid = false;
        }
    }
    template <std::uint_fast8_t IGNOREIDX = 0>
    void voxelTransport(Particle<T>& p, RandomState& state)
    {

        bool still_inside;
        do {
            std::array<std::size_t, 3> xyz = index<false>(p.pos);
            auto index_flat = flatIndex(xyz);

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
                m_spacing[0] / std::abs(p.dir[0]),
                m_spacing[1] / std::abs(p.dir[1]),
                m_spacing[2] / std::abs(p.dir[2])
            };

            std::array<T, 3> tMax = {
                // abs function in case of delta is infinity
                std::abs(std::nextafter(delta[0], std::numeric_limits<T>::max()) - delta[0]),
                std::abs(std::nextafter(delta[1], std::numeric_limits<T>::max()) - delta[1]),
                std::abs(std::nextafter(delta[2], std::numeric_limits<T>::max()) - delta[2])
            };

            // preventing zero threshold
            const T interaction_thres = std::max(state.randomUniform<T>(), std::numeric_limits<T>::min());

            T interaction_accum = 1;
            bool cont = true;
            std::uint_fast8_t dIdx;
            do {
                // current material
                const auto matInd = m_data[index_flat].materialIndex;
                const auto density = m_data[index_flat].density;
                auto& dose = m_data[index_flat].dose;

                // updating stepping
                dIdx = argmin3(tMax);
                xyz[dIdx] += xyz_step[dIdx];
                index_flat += xyz_step_flat[dIdx];
                tMax[dIdx] += delta[dIdx];

                // we are still inside volume and not in an ignored voxel
                still_inside = xyz[dIdx] < m_dim[dIdx] && m_data[index_flat].materialIndex != IGNOREIDX;

                // updating radiological path
                const auto att = m_materials[matInd].attenuationValues(p.energy);
                const auto att_tot = att.sum() * density;
                const auto radio_step = std::exp(-att_tot * delta[dIdx]);
                interaction_accum *= radio_step;
                if (interaction_accum < interaction_thres) {
                    // interaction happends
                    const auto step_correction = std::log(interaction_accum / interaction_thres) / att_tot;
                    tMax[dIdx] += step_correction;
                    const auto intRes = interactions::template interact<T, NMaterialShells, LOWENERGYCORRECTION>(att, p, m_materials[matInd], state);
                    dose.scoreEnergy(intRes.energyImparted);
                    still_inside = still_inside && intRes.particleAlive;
                    // particle energy or direction has changed, we restart stepping
                    cont = false;
                } else {
                    cont = still_inside;
                }
            } while (cont);
            // advancing particle to reached step
            p.translate(tMax[dIdx]);
        } while (still_inside);
    }
    void woodcockTransport(Particle<T>& p, RandomState& state)
    {
    }

    struct DataElement {
        DoseScore<T> dose;
        T density = 0;
        std::uint8_t materialIndex = 0;
    };

private:
    std::array<std::size_t, 3> m_dim = { 1, 1, 1 };
    std::array<T, 3> m_invSpacing = { 1, 1, 1 };
    std::array<T, 3> m_spacing = { 1, 1, 1 };
    std::array<T, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };
    std::vector<DataElement> m_data;
    std::vector<Material2<T, NMaterialShells>> m_materials;
};
}