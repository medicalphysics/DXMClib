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
#include "dxmc/world/worlditems/tetrahedalmesh/tetrahedalmeshreader.hpp"
#include "dxmc/world/worlditems/tetrahedalmesh/tetrahedron.hpp"
#include "dxmc/world/worlditems/worlditembase.hpp"

#include <array>
#include <string>
#include <string_view>
#include <vector>

namespace dxmc {

template <Floating T, int NMaterialShells = 5, int LOWENERGYCORRECTION = 2>
class ThetrahedalMesh final : public WorldItemBase<T> {
public:
    ThetrahedalMesh()
        : WorldItemBase<T>()
    {
    }

    void readICRP145Phantom(const std::string& nodeFile, const std::string& elementsFile)
    {
        TetrahedalmeshReader<T> reader;
        m_kdtree.setData(reader.readICRP145Phantom(nodeFile, elementsFile));
        m_aabb = m_kdtree.AABB();
        return;
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
        WorldIntersectionResult<T> res;
        return res;
    }

    VisualizationIntersectionResult<T, WorldItemBase<T>> intersectVisualization(const Particle<T>& p) const override
    {
        VisualizationIntersectionResult<T, WorldItemBase<T>> res;
        return res;
    }

    const DoseScore<T>& dose(std::size_t index = 0) const override
    {
        return DoseScore<T> {};
    }

    void clearDose() override { }

    void transport(Particle<T>& p, RandomState& state) override { }

protected:
private:
    std::array<T, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };
    TetrahedalMeshKDTree<T> m_kdtree;
};
}
