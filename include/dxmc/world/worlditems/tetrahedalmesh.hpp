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
#include "dxmc/interactions.hpp"
#include "dxmc/material/material.hpp"
#include "dxmc/material/nistmaterials.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/worlditems/tetrahedalmesh/tetrahedalmeshkdtree.hpp"
#include "dxmc/world/worlditems/tetrahedalmesh/tetrahedron.hpp"
#include "dxmc/world/worlditems/worlditembase.hpp"

#include <array>
#include <string>
#include <string_view>
#include <vector>

namespace dxmc {

template <Floating T, int NMaterialShells = 5, int LOWENERGYCORRECTION = 2>
class TetrahedalMesh final : public WorldItemBase<T> {
public:
    TetrahedalMesh()
        : WorldItemBase<T>()
    {
    }

    void translate(const std::array<T, 3>& dist) override
    {
        m_kdtree.translate(dist);
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] += dist[i];
            m_aabb[i + 3] += dist[i];
        }
    }

    std::array<T, 3> center() const override
    {
        std::array<T, 3> c;
        for (std::size_t i = 0; i < 3; ++i) {
            c[i] = (m_aabb[i] + m_aabb[i + 3]) / 2;
        }
        return c;
    }

    std::array<T, 6> AABB() const override
    {
        return m_aabb;
    }
    WorldIntersectionResult<T> intersect(const Particle<T>& p) const override
    {
        const auto res = m_kdtree.intersect(p, m_aabb);
        WorldIntersectionResult<T> w;
        if (res.valid()) {
            w.intersection = res.intersection;
            w.intersectionValid = true;
            w.rayOriginIsInsideItem = res.rayOriginIsInsideItem;
        }
        return w;
    }

    VisualizationIntersectionResult<T, WorldItemBase<T>> intersectVisualization(const Particle<T>& p) const override
    {
        const auto res = m_kdtree.intersect<6100>(p, m_aabb);
        VisualizationIntersectionResult<T, WorldItemBase<T>> w;
        if (res.valid()) {
            w.intersection = res.intersection;
            w.rayOriginIsInsideItem = res.rayOriginIsInsideItem;
            w.intersectionValid = true;
            w.item = this;
            w.normal = res.normal;
        }
        return w;
    }

    const DoseScore<T>& dose(std::size_t index = 0) const override
    {
        return DoseScore<T> {};
    }

    void clearDose() override { }

    void transport(Particle<T>& p, RandomState& state) override
    {
    }

protected:
    void expandAABB(T extra = 0.0001f)
    {
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] -= extra;
            m_aabb[i + 3] += extra;
        }
    }

private:
    struct Collection {
        DoseScore<T> dose;
        T density = 0;
        std::uint16_t materialIdx = 0;
        std::uint16_t organIdx = 0;
        Collection(T dens)
            : density(dens)
        {
        }
    };

    std::array<T, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };
    TetrahedalMeshKDTree<T> m_kdtree;
    std::vector<Collection> m_collections;
    std::vector<Material2<T, NMaterialShells>> m_materials;
    std::vector<std::string> m_organNames;
};
}
