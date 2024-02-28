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
class TriangulatedMesh {
public:
    TriangulatedMesh()
        : m_materialDensity(NISTMaterials::density("Air, Dry (near sea level)"))
        , m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
    }

    TriangulatedMesh(const std::vector<Triangle>& triangles, const std::size_t max_tree_dept = 8)
        : m_materialDensity(NISTMaterials::density("Air, Dry (near sea level)"))
        , m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        setData(triangles, max_tree_dept);
    }

    TriangulatedMesh(const std::string& path, const std::size_t max_tree_dept = 8)
        : m_materialDensity(NISTMaterials::density("Air, Dry (near sea level)"))
        , m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
    {
        STLReader reader;
        auto triangles = reader(path);
        setData(triangles, max_tree_dept);
    }

    TriangulatedMesh(const std::string& path, double scale, const std::size_t max_tree_dept = 8)
        : m_material(Material<NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
        , m_materialDensity(NISTMaterials::density("Air, Dry (near sea level)"))
    {
        STLReader reader;
        auto triangles = reader(path);
        std::for_each(std::execution::par_unseq, triangles.begin(), triangles.end(), [=](auto& tri) {
            tri.scale(scale);
        });
        setData(triangles, max_tree_dept);
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

    std::vector<Triangle> getTriangles() const
    {
        return m_kdtree.items();
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

    void rotate(double radians, const std::array<double, 3>& axis)
    {
        const auto axisn = vectormath::normalized(axis);
        const auto c = center();
        const auto cm = vectormath::scale(c, -1.0);
        m_kdtree.translate(cm);
        m_kdtree.rotate(radians, axisn);
        m_kdtree.translate(c);
        m_aabb = expandAABB(m_kdtree.AABB());
    }

    void setData(const std::vector<Triangle>& triangles, const std::size_t max_tree_dept = 8)
    {
        m_kdtree.setData(triangles, max_tree_dept);
        m_aabb = expandAABB(m_kdtree.AABB());
        m_volume = calculateVolume(triangles, m_aabb);
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
        const auto res = m_kdtree.intersect(p, m_aabb);
        return WorldIntersectionResult { .intersection = res.intersection, .rayOriginIsInsideItem = res.rayOriginIsInsideItem, .intersectionValid = res.item != nullptr };
    }

    template <typename U>
    VisualizationIntersectionResult<U> intersectVisualization(const ParticleType auto& p) const noexcept
    {
        const auto res = m_kdtree.intersect(p, m_aabb);
        VisualizationIntersectionResult<U> res_int;
        if (res.valid()) {
            res_int.normal = vectormath::normalized(res.item->planeVector());
            res_int.intersection = res.intersection;
            res_int.rayOriginIsInsideItem = res.rayOriginIsInsideItem;
            res_int.intersectionValid = true;
            res_int.value = m_dose.dose();
        }
        return res_int;
    }

    void transport(ParticleType auto& p, RandomState& state)
    {
        bool cont = true;
        bool updateAtt = true;
        AttenuationValues att;
        double attSumInv;
        while (cont) {
            if (updateAtt) {
                att = m_material.attenuationValues(p.energy);
                attSumInv = 1 / (att.sum() * m_materialDensity);
                updateAtt = false;
            }

            const auto intM = intersect(p);
            const auto stepLen = -std::log(state.randomUniform()) * attSumInv;
            if (intM.intersection < stepLen) {
                // transport to edge
                p.border_translate(intM.intersection);
                cont = false;
            } else {
                p.translate(stepLen);
                const auto intRes = interactions::template interact<NMaterialShells, LOWENERGYCORRECTION>(att, p, m_material, state);
                updateAtt = intRes.particleEnergyChanged;
                cont = intRes.particleAlive;
                if (intRes.particleEnergyChanged) {
                    m_energyScored.scoreEnergy(intRes.energyImparted);
                }
            }
        }
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
        // not defined for fleunce counter
        m_dose.addScoredEnergy(m_energyScored, m_volume, m_materialDensity, calibration_factor);
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
    static std::array<double, 6> expandAABB(std::array<double, 6> aabb)
    {
        // note we take copy of aabb
        for (std::size_t i = 0; i < 3; ++i) {
            aabb[i] -= GEOMETRIC_ERROR();
            aabb[i + 3] += GEOMETRIC_ERROR();
        }
        return aabb;
    }

    static double calculateVolume(const std::vector<Triangle>& triangles, const std::array<double, 6>& aabb)
    {
        // using signed volume of a thetrahedron by a point and the three vertices (Gauss theorem of divergence)
        // sign_vol = v1 * (v2 x v3) / 6 for three vectors defining a triangle, i.e triple product
        // volume is then the sum of signed volumes over all triangles

        // finding a neat reference point instead of {0,0,0}
        const std::array center = {
            (aabb[0] + aabb[3]) / 2,
            (aabb[1] + aabb[4]) / 2,
            (aabb[2] + aabb[5]) / 2
        };

        const auto volume = std::transform_reduce(std::execution::par_unseq, triangles.cbegin(), triangles.cend(), 0.0, std::plus<>(), [&center](const auto& tri) {
            const auto& v = tri.vertices();
            const auto v1 = vectormath::subtract(v[0], center);
            const auto v2 = vectormath::subtract(v[1], center);
            const auto v3 = vectormath::subtract(v[2], center);
            return vectormath::tripleProduct(v1, v2, v3);
        });
        return std::abs(volume) / 6;
    }

private:
    double m_materialDensity = 1;
    double m_volume = 0; // cm3
    std::array<double, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };
    EnergyScore m_energyScored;
    DoseScore m_dose;
    MeshKDTree<Triangle> m_kdtree;
    Material<NMaterialShells> m_material;
};
}