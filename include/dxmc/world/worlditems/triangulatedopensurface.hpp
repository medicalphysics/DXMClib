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

#include "dxmc/interactions.hpp"
#include "dxmc/material/material.hpp"
#include "dxmc/material/nistmaterials.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/worlditems/triangulatedmesh/triangle.hpp"
#include "dxmc/world/worlditems/triangulatedmesh/triangulatedmeshkdtree.hpp"
#include "dxmc/world/worlditems/triangulatedmesh/triangulatedmeshstlreader.hpp"

#include <array>
#include <execution>
#include <numeric>
#include <vector>

namespace dxmc {

template <int NMaterialShells = 5, int LOWENERGYCORRECTION = 2>
class TriangulatedOpenSurface {
public:
    TriangulatedOpenSurface()
        : m_materialDensity(NISTMaterials::density("Air, Dry (near sea level)"))
        , m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
    }

    TriangulatedOpenSurface(const std::vector<Triangle>& triangles, double surfaceThickness = 0.035, const std::size_t max_tree_dept = 8)
        : m_materialDensity(NISTMaterials::density("Air, Dry (near sea level)"))
        , m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        setData(triangles, surfaceThickness, max_tree_dept);
    }

    TriangulatedOpenSurface(const std::string& path, double surfaceThickness = 0.035, const std::size_t max_tree_dept = 8)
        : m_materialDensity(NISTMaterials::density("Air, Dry (near sea level)"))
        , m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        STLReader reader;
        auto triangles = reader(path);
        setData(triangles, surfaceThickness, max_tree_dept);
    }

    double surfaceThickness() const
    {
        return m_thickness;
    }
    void setSurfaceThickness(double cm)
    {
        m_thickness = std::abs(cm);
    }

    void setMaterial(const Material<NMaterialShells>& mat)
    {
        m_material = mat;
    }

    void setMaterialDensity(double dens)
    {
        m_materialDensity = std::abs(dens);
    }

    void setMaterial(const Material<NMaterialShells>& mat, double dens)
    {
        m_material = mat;
        setMaterialDensity(dens);
    }

    bool setNistMaterial(const std::string& nist_name)
    {
        const auto mat = Material<NMaterialShells>::byNistName(nist_name);
        if (mat) {
            m_material = mat.value();
            m_materialDensity = NISTMaterials::density(nist_name);
            return true;
        }
        return false;
    }

    const MeshKDTree<Triangle>& kdtree() const
    {
        return m_kdtree;
    }

    const std::vector<Triangle>& getTriangles() const
    {
        return m_triangles;
    }

    void translate(const std::array<double, 3>& dist)
    {
        m_kdtree.translate(dist);
        for (std::size_t i = 0; i < 3; ++i) {
            m_aabb[i] += dist[i];
            m_aabb[i + 3] += dist[i];
        }
    }

    void scale(double s)
    {
        m_kdtree.scale(s);
        for (auto& e : m_aabb)
            e *= s;
    }

    void mirror(const std::array<double, 3>& point)
    {
        m_kdtree.mirror(point);
        calculateAABB();
    }

    void mirror(const double value, const std::uint_fast32_t dim)
    {
        m_kdtree.mirror(value, dim);
        calculateAABB();
    }

    void setData(const std::vector<Triangle>& triangles, double surfaceThickness = 0.035, const std::size_t max_tree_dept = 8)
    {

        m_thickness = std::abs(surfaceThickness);
        m_triangles = triangles;
        m_kdtree.setData(m_triangles, max_tree_dept);
        calculateAABB();
    }

    std::array<double, 3> center() const
    {
        std::array<double, 3> center { 0, 0, 0 };
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

    WorldIntersectionResult intersect(const ParticleType auto& p) const
    {
        const auto res = m_kdtree.intersect(p, m_triangles, m_aabb);
        return WorldIntersectionResult { .intersection = res.intersection, .rayOriginIsInsideItem = false, .intersectionValid = res.item != nullptr };
    }

    template <typename U>
    VisualizationIntersectionResult<U> intersectVisualization(const ParticleType auto& p) const noexcept
    {
        const auto res = m_kdtree.intersect(p, m_triangles, m_aabb);
        VisualizationIntersectionResult<U> res_int;
        if (res.valid()) {
            res_int.normal = vectormath::normalized(res.item->planeVector());
            if (res.rayOriginIsInsideItem) // fix normal vector
                res_int.normal = vectormath::scale(res_int.normal, -1.0);
            res_int.intersection = res.intersection;
            res_int.rayOriginIsInsideItem = false;
            res_int.intersectionValid = true;
            res_int.value = m_dose.dose();
        }
        return res_int;
    }

    void transport(ParticleType auto& p, RandomState& state)
    {
        const auto att = m_material.attenuationValues(p.energy);
        const auto attSumInv = 1 / (att.sum() * m_materialDensity);
        const auto stepLen = -std::log(state.randomUniform()) * attSumInv;

        // The particle is transported right behind this plane, we correct this
        p.translate(-2 * GEOMETRIC_ERROR());

        const auto intersection = m_kdtree.intersect(p, m_triangles, m_aabb);
        const auto normal = intersection.item->planeVector();
        const auto length_scale = std::abs(vectormath::dot(p.dir, normal));
        const auto length = m_thickness / length_scale;

        if (stepLen < length) {
            const auto intRes = interactions::template interact<NMaterialShells, LOWENERGYCORRECTION>(att, p, m_material, state);
            if (intRes.particleEnergyChanged) {
                m_energyScored.scoreEnergy(intRes.energyImparted);
            }
        }
        p.border_translate(length * 0.5);
    }

    const EnergyScore& energyScored(std::size_t index = 0) const
    {
        return m_energyScored;
    }

    void clearEnergyScored()
    {
        m_energyScored.clear();
    }

    void addEnergyScoredToDoseScore(double calibration_factor = 1)
    {
        const auto volume = calculateVolume(m_triangles, m_thickness);
        m_dose.addScoredEnergy(m_energyScored, volume, m_materialDensity, calibration_factor);
    }

    const DoseScore& doseScored(std::size_t index = 0) const
    {
        return m_dose;
    }

    void clearDoseScored()
    {
        m_dose.clear();
    }

    const std::array<double, 6>& AABB() const
    {
        return m_aabb;
    }

protected:
    void calculateAABB()
    {
        std::array<double, 6> aabb {
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
            std::numeric_limits<double>::lowest(),
        };

        for (const auto& tri : m_triangles) {
            const auto aabb_tri = tri.AABB();

            for (std::size_t i = 0; i < 3; ++i) {
                aabb[i] = std::min(aabb[i], aabb_tri[i]);
            }
            for (std::size_t i = 3; i < 6; ++i) {
                aabb[i] = std::max(aabb[i], aabb_tri[i]);
            }
        }

        // Extending AABB
        for (std::size_t i = 0; i < 3; ++i) {
            aabb[i] -= GEOMETRIC_ERROR();
            aabb[i + 3] += GEOMETRIC_ERROR();
        }
        m_aabb = aabb;
    }

    static double calculateVolume(const std::vector<Triangle>& triangles, double thickness)
    {
        const auto volume = std::transform_reduce(std::execution::par_unseq, triangles.cbegin(), triangles.cend(), 0.0, std::plus<>(), [thickness](const auto& tri) {
            return tri.area() * thickness;
        });
        return volume;
    }

private:
    double m_materialDensity = 1;
    double m_thickness = 0.035; // cm
    std::array<double, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };
    EnergyScore m_energyScored;
    DoseScore m_dose;
    MeshKDTree<Triangle> m_kdtree;
    std::vector<Triangle> m_triangles;
    Material<NMaterialShells> m_material;
};
}