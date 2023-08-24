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
class TriangulatedMesh final : public WorldItemBase<T> {
public:
    TriangulatedMesh()
        : WorldItemBase<T>()
        , m_material(Material<T, NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
        , m_materialDensity(NISTMaterials<T>::density("Air, Dry (near sea level)"))
    {
    }

    TriangulatedMesh(const std::vector<Triangle<T>>& triangles, const std::size_t max_tree_dept = 8)
        : WorldItemBase<T>()
        , m_material(Material<T, NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
        , m_materialDensity(NISTMaterials<T>::density("Air, Dry (near sea level)"))
    {
        setData(triangles, max_tree_dept);
    }

    TriangulatedMesh(const std::string& path, const std::size_t max_tree_dept = 8)
        : WorldItemBase<T>()
        , m_material(Material<T, NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
        , m_materialDensity(NISTMaterials<T>::density("Air, Dry (near sea level)"))
    {
        STLReader<T> reader;
        auto triangles = reader(path);
        setData(triangles, max_tree_dept);
    }

    TriangulatedMesh(const std::string& path, T scale, const std::size_t max_tree_dept = 8)
        : WorldItemBase<T>()
        , m_material(Material<T, NMaterialShells>::byNistName("Air, Dry (near sea level)").value())
        , m_materialDensity(NISTMaterials<T>::density("Air, Dry (near sea level)"))
    {
        STLReader<T> reader;
        auto triangles = reader(path);
        std::for_each(std::execution::par_unseq, triangles.begin(), triangles.end(), [=](auto& tri) {
            tri.scale(scale);
        });
        setData(triangles, max_tree_dept);
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

    void setData(const std::vector<Triangle<T>>& triangles, const std::size_t max_tree_dept)
    {
        m_kdtree.setData(triangles, max_tree_dept);
        m_aabb = expandAABB(m_kdtree.AABB());
        m_volume = calculateVolume(triangles, m_aabb);
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
        auto res = m_kdtree.intersect(p, m_aabb);
        return WorldIntersectionResult<T> { .intersection = res.intersection, .rayOriginIsInsideItem = res.rayOriginIsInsideItem, .intersectionValid = res.item != nullptr };
    }

    VisualizationIntersectionResult<T, WorldItemBase<T>> intersectVisualization(const Particle<T>& p) const noexcept override
    {
        auto res = m_kdtree.intersect(p, m_aabb);
        VisualizationIntersectionResult<T, WorldItemBase<T>> res_int;
        if (res.valid()) {
            auto normal = res.item->planeVector();
            vectormath::normalize(normal);
            res_int.normal = normal;
            res_int.intersection = res.intersection;
            res_int.rayOriginIsInsideItem = res.rayOriginIsInsideItem;
            res_int.intersectionValid = true;
        }
        return res_int;
    }

    void transport(Particle<T>& p, RandomState& state) override
    {
        bool cont = true;
        bool updateAtt = true;
        dxmc::AttenuationValues<T> att;
        T attSumInv;
        while (cont) {
            if (updateAtt) {
                att = m_material.attenuationValues(p.energy);
                attSumInv = 1 / (att.sum() * m_materialDensity);
                updateAtt = false;
            }

            const auto intM = intersect(p);
            const auto stepLen = -std::log(state.randomUniform<T>()) * attSumInv;
            if (intM.intersection < stepLen) {
                // transport to edge
                p.border_translate(intM.intersection);
                cont = false;
            } else {
                p.translate(stepLen);
                const auto intRes = interactions::template interact<T, NMaterialShells, LOWENERGYCORRECTION>(att, p, m_material, state);
                updateAtt = intRes.particleEnergyChanged;
                cont = intRes.particleAlive;
                if (intRes.particleEnergyChanged) {
                    m_energyScored.scoreEnergy(intRes.energyImparted);
                }
            }
        }
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
        // not defined for fleunce counter
        m_dose.addScoredEnergy(m_energyScored, m_volume, m_materialDensity, calibration_factor);
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

    static T calculateVolume(const std::vector<Triangle<T>>& triangles, const std::array<T, 6>& aabb)
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

        const auto volume = std::transform_reduce(std::execution::par_unseq, triangles.cbegin(), triangles.cend(), T { 0 }, std::plus<>(), [&center](const auto& tri) {
            const auto& v = tri.vertices();
            const auto v1 = vectormath::subtract(v[0], center);
            const auto v2 = vectormath::subtract(v[1], center);
            const auto v3 = vectormath::subtract(v[2], center);
            return vectormath::tripleProduct(v1, v2, v3);
        });
        return std::abs(volume) / 6;
    }

private:
    T m_materialDensity = 1;
    T m_volume = 0; // cm3
    std::array<T, 6> m_aabb = { 0, 0, 0, 0, 0, 0 };
    EnergyScore<T> m_energyScored;
    DoseScore<T> m_dose;
    MeshKDTree<T, Triangle<T>> m_kdtree;
    Material<T, NMaterialShells> m_material;
};
}