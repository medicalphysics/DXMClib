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
#include "dxmc/world/triangulatedmesh/meshkdtree.hpp"
#include "dxmc/world/triangulatedmesh/stlreader.hpp"
#include "dxmc/world/triangulatedmesh/triangle.hpp"
#include "dxmc/world/worlditembase.hpp"

#include <array>
#include <vector>

namespace dxmc {

template <Floating T>
class TriangulatedMesh final : public WorldItemBase<T> {
public:
    TriangulatedMesh()
        : WorldItemBase<T>()
    {
    }
    TriangulatedMesh(const std::vector<Triangle<T>>& triangles, const std::size_t max_tree_dept = 8)
        : WorldItemBase<T>()
    {
        setData(triangles, max_tree_dept);
    }
    TriangulatedMesh(const std::string& path, const std::size_t max_tree_dept = 8)
        : WorldItemBase<T>()
    {
        STLReader<T> reader;
        auto triangles = reader(path);
        setData(triangles, max_tree_dept);
    }

    const MeshKDTree<T, Triangle<T>>& kdtree() const
    {
        return m_kdtree;
    }

    std::vector<Triangle<T>> getTriangles() const
    {
        return m_kdtree.items();
    }
    void translate(const std::array<T, 3>& dist) override
    {
        m_kdtree.translate(dist);
        m_aabb = m_kdtree.AABB();
    }

    void setData(std::vector<Triangle<T>> triangles, const std::size_t max_tree_dept)
    {
        m_kdtree.setData(triangles, max_tree_dept);
        m_aabb = m_kdtree.AABB();
    }
    std::array<T, 3> center() const override
    {
        std::array<T, 3> center { 0, 0, 0 };
        const auto triangles = getTriangles();
        std::for_each(std::execution::unseq, triangles.cbegin(), triangles.cend(), [&](const auto& tri) {
            const auto c = tri.center();
            for (std::size_t i = 0; i < 3; ++i) {
                center[i] += c[i];
            }
        });
        for (std::size_t i = 0; i < 3; ++i) {
            center[i] /= triangles.size();
        }
        return center;
    }
    std::optional<T> intersect(const Particle<T>& p) const override
    {
        return m_kdtree.intersect(p, m_aabb);
    }

    T transport(Particle<T>& p, RandomState& state) override
    {
        return 0;
    }

    std::array<T, 6> AABB() const override
    {
        return m_aabb;
    }

private:
    std::array<T, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };
    MeshKDTree<T, Triangle<T>> m_kdtree;
};
}