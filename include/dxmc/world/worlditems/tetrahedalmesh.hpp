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
class TetrahedalMesh final : public WorldItemBase<T> {
public:
    TetrahedalMesh()
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

    void setData(const std::vector<std::array<std::size_t, 4>>& nodes, const std::vector<std::array<T, 3>>& vertices)
    {
        const auto max_ind = std::transform_reduce(
            std::execution::par_unseq, nodes.cbegin(), nodes.cend(), std::size_t { 0 }, [](auto lh, auto rh) { return std::max(lh, rh); }, [](const auto& t) {
            std::size_t max = t[3];
            for (std::size_t i = 0; i < 3; ++i)
                max = std::max(max, t[i]);
            return max; });
        if (max_ind < vertices.size()) {
            std::vector<Tetrahedron<T>> tets(nodes.size());
            std::transform(
                std::execution::par_unseq, nodes.cbegin(), nodes.cend(), tets.begin(), [&vertices](const auto& n) {
                    std::array<std::array<T, 3>, 4> verts;
                    for (std::size_t i = 0; i < 4; ++i) {
                        verts[i] = vertices[n[i]];
                    }
                    return Tetrahedron { verts };
                });
            m_kdtree.setData(std::move(tets));
            m_aabb = m_kdtree.AABB();
        }
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
        return m_kdtree.intersect<WorldIntersectionResult<T>>(p, m_aabb);
    }

    VisualizationIntersectionResult<T, WorldItemBase<T>> intersectVisualization(const Particle<T>& p) const override
    {
        auto res = m_kdtree.intersect<VisualizationIntersectionResult<T, WorldItemBase<T>>>(p, m_aabb);
        res.item = this;
        return res;
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
private:
    std::array<T, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };
    TetrahedalMeshKDTree<T> m_kdtree;
};
}
