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
#include "dxmc/world/worlditems/triangulatedmesh/triangle.hpp"
#include "dxmc/world/worlditems/triangulatedmesh/triangulatedmeshkdtree.hpp"
#include "dxmc/world/worlditems/triangulatedmesh/triangulatedmeshstlreader.hpp"
#include "dxmc/world/worlditems/worlditembase.hpp"

#include <array>
#include <execution>
#include <numeric>
#include <vector>

namespace dxmc {

template <Floating T, int NMaterialShells = 5, int LOWENERGYCORRECTION = 2>
class TriangulatedOpenSurface final : public WorldItemBase<T> {
public:
    TriangulatedOpenSurface()
        : WorldItemBase<T>()
        , m_materialDensity(NISTMaterials<T>::density("Air, Dry (near sea level)"))
        , m_material(Material<T, NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
    }

    TriangulatedOpenSurface(const std::vector<Triangle<T>>& triangles, T surfaceThickness = 0.035, const std::size_t max_tree_dept = 8)
        : WorldItemBase<T>()
        , m_materialDensity(NISTMaterials<T>::density("Air, Dry (near sea level)"))
        , m_material(Material<T, NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        setData(triangles, surfaceThickness, max_tree_dept);
    }

    TriangulatedOpenSurface(const std::string& path, T surfaceThickness = 0.035, const std::size_t max_tree_dept = 8)
        : WorldItemBase<T>()
        , m_materialDensity(NISTMaterials<T>::density("Air, Dry (near sea level)"))
        , m_material(Material<T, NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        STLReader<T> reader;
        auto triangles = reader(path);
        setData(triangles, surfaceThickness, max_tree_dept);
    }

    T surfaceAttenuationThickness() const
    {
        return m_thickness;
    }

    void setMaterial(const Material<T, NMaterialShells>& mat)
    {
        m_material = mat;
    }

    void setMaterialDensity(T dens)
    {
        m_materialDensity = std::abs(dens);
    }

    void setMaterial(const Material<T, NMaterialShells>& mat, T dens)
    {
        m_material = mat;
        setMaterialDensity(dens);
    }

    bool setNistMaterial(const std::string& nist_name)
    {
        const auto mat = Material<T, NMaterialShells>::byNistName(nist_name);
        if (mat) {
            m_material = mat.value();
            m_materialDensity = NISTMaterials<T>::density(nist_name);
            return true;
        }
        return false;
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
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] += dist[i];
            m_aabb[i + 3] += dist[i];
        }
    }

    void scale(T s)
    {
        m_kdtree.scale(s);
        for (auto& e : m_aabb)
            e *= s;
    }

    void rotate(T radians, const std::array<T, 3>& axis)
    {
        const auto axisn = vectormath::normalized(axis);
        const auto c = center();
        const auto cm = vectormath::scale(c, T { -1 });
        m_kdtree.translate(cm);
        m_kdtree.rotate(radians, axisn);
        m_kdtree.translate(c);
        m_aabb = expandAABB(m_kdtree.AABB());
    }

    void setData(const std::vector<Triangle<T>>& triangles, T surfaceThickness = 0.035, const std::size_t max_tree_dept = 8)
    {
        m_thickness = std::abs(surfaceThickness);
        m_kdtree.setData(triangles, max_tree_dept);
        m_aabb = expandAABB(m_kdtree.AABB());
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

    WorldIntersectionResult<T> intersect(const Particle<T>& p) const override
    {
        const auto res = m_kdtree.intersect(p, m_aabb);
        return WorldIntersectionResult<T> { .intersection = res.intersection, .rayOriginIsInsideItem = false, .intersectionValid = res.item != nullptr };
    }

    VisualizationIntersectionResult<T, WorldItemBase<T>> intersectVisualization(const Particle<T>& p) const noexcept override
    {
        const auto res = m_kdtree.intersect(p, m_aabb);
        VisualizationIntersectionResult<T, WorldItemBase<T>> res_int;
        if (res.valid()) {
            res_int.normal = vectormath::normalized(res.item->planeVector());
            if (res.rayOriginIsInsideItem) // fix normal vector
                res_int.normal = vectormath::scale(res_int.normal, T { -1 });
            res_int.intersection = res.intersection;
            res_int.rayOriginIsInsideItem = false;
            res_int.intersectionValid = true;
            res_int.value = m_dose.dose();
        }
        return res_int;
    }

    void transport(Particle<T>& p, RandomState& state) override
    {
        const dxmc::AttenuationValues<T> att = m_material.attenuationValues(p.energy);
        const auto attSumInv = 1 / (att.sum() * m_materialDensity);
        const auto stepLen = -std::log(state.randomUniform<T>()) * attSumInv;

        const auto intersection = m_kdtree.intersect(p, m_aabb);
        const auto normal = intersection.item->planeVector();
        const auto length_scale = std::abs(vectormath::dot(p.dir, normal));
        const auto length = m_thickness / length_scale;

        if (stepLen < length) {
            const auto intRes = interactions::template interact<T, NMaterialShells, LOWENERGYCORRECTION>(att, p, m_material, state);
            if (intRes.particleEnergyChanged) {
                m_energyScored.scoreEnergy(intRes.energyImparted);
            }
        }
        p.border_translate(length * T { 0.5 });
    }

    const EnergyScore<T>& energyScored(std::size_t index = 0) const override
    {
        return m_energyScored;
    }

    void clearEnergyScored() override
    {
        m_energyScored.clear();
    }

    void addEnergyScoredToDoseScore(T calibration_factor = 1) override
    {
        const auto triangles = getTriangles();
        const auto volume = calculateVolume(triangles, m_thickness);
        m_dose.addScoredEnergy(m_energyScored, volume, m_materialDensity, calibration_factor);
    }

    const DoseScore<T>& doseScored(std::size_t index = 0) const override
    {
        return m_dose;
    }

    void clearDoseScored() override
    {
        m_dose.clear();
    }

    std::array<T, 6> AABB() const override
    {
        return m_aabb;
    }

protected:
    static std::array<T, 6> expandAABB(std::array<T, 6> aabb)
    {
        // note we take copy of aabb
        for (std::size_t i = 0; i < 3; ++i) {
            aabb[i] -= GEOMETRIC_ERROR<T>();
            aabb[i + 3] += GEOMETRIC_ERROR<T>();
        }
        return aabb;
    }

    static T calculateVolume(const std::vector<Triangle<T>>& triangles, T thickness)
    {
        const auto volume = std::transform_reduce(std::execution::par_unseq, triangles.cbegin(), triangles.cend(), T { 0 }, std::plus<>(), [thickness](const auto& tri) {
            return tri.area() * thickness;
        });
        return volume;
    }

private:
    T m_materialDensity = 1;
    T m_thickness = 0.035; // cm
    std::array<T, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };
    EnergyScore<T> m_energyScored;
    DoseScore<T> m_dose;
    MeshKDTree<T, Triangle<T>> m_kdtree;
    Material<T, NMaterialShells> m_material;
};
}