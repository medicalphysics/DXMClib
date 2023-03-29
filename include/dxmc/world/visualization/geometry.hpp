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
#include "dxmc/world/kdtreeintersectionresult.hpp"
#include "dxmc/world/world.hpp"
#include "dxmc/world/worlditems/worlditembase.hpp"

#include <array>
#include <vector>

namespace dxmc {
namespace visualization {

    template <typename U, typename T>
    concept WorldType = requires(U world, Particle<T> p, KDTreeIntersectionResult<T, WorldItemBase<T>> res) {
                            Floating<T>;
                            {
                                world.intersect(p)
                                } -> std::same_as<KDTreeIntersectionResult<T, WorldItemBase<T>>>;
                        };

    template <Floating T>
    std::array<T, 3> projectToPlane(const std::array<T, 3>& planePos, const std::array<T, 3>& planeNormal, const std::array<T, 3>& point)
    {

        const auto np = vectormath::dot(planeNormal, point);
        const auto no = vectormath::dot(planeNormal, planePos);
        const auto d = vectormath::scale(np - no, planeNormal);
        return vectormath::subtract(planePos, d);
    }

    template <Floating T>
    std::pair<T, T> fieldOfView(const std::array<T, 6>& aabb, const std::array<T, 3>& cameraPos, const std::array<T, 3>& cameraCosinex, const std::array<T, 3>& cameraCosiney)
    {
        const auto dir = vectormath::cross(cameraCosinex, cameraCosiney);
        const std::array<T, 3> p1 = {
            (aabb[0] + aabb[3]) / 2,
            (aabb[1] + aabb[4]) / 2,
            (aabb[2] + aabb[5]) / 2
        };
        const auto [c0, c2] = vectormath::splice(aabb);

        const auto p0 = projectToPlane(p1, dir, c0);
        const auto p2 = projectToPlane(p1, dir, c2);

        const auto b0 = vectormath::changeBasis(cameraCosinex, cameraCosiney, dir, p0);
        const auto b1 = vectormath::changeBasis(cameraCosinex, cameraCosiney, dir, p1);
        const auto b2 = vectormath::changeBasis(cameraCosinex, cameraCosiney, dir, p2);

        auto test = vectormath::changeBasis(cameraCosinex, cameraCosiney, dir, dir);

        const auto fov2x = vectormath::lenght(vectormath::subtract(b0, b1));
        const auto fov2y = vectormath::lenght(vectormath::subtract(b2, b1));
        return std::make_pair(fov2x * 2, fov2y * 2);
    }

    template <Floating T, WorldType<T> W>
    std::vector<T> rayTraceGeometryDistance(W& world, const std::array<T, 3>& cameraPos, const std::array<T, 3>& cameraCosinex, const std::array<T, 3>& cameraCosiney, std::size_t resolution = 128)
    {
        const auto aabb = world.AABB();
        const std::array<T, 3> c = {
            (aabb[0] + aabb[3]) / 2,
            (aabb[1] + aabb[4]) / 2,
            (aabb[2] + aabb[5]) / 2
        };

        const auto [fov_x, fov_y] = fieldOfView(world.AABB(), cameraPos, cameraCosinex, cameraCosiney);
        const auto fov = std::max(fov_x, fov_y);
        const auto pc = vectormath::subtract(c, cameraPos);

        const auto dd = std::tan(fov / resolution / vectormath::lenght(pc) );
        const auto dd_start = -dd * resolution / 2;

        const auto dir = vectormath::cross(cameraCosinex, cameraCosiney);

        std::vector<T> data(resolution * resolution, T { 0 });

        for (std::size_t y = 0; y < resolution; ++y) {
            const auto y_flat = y * resolution;
            for (std::size_t x = 0; x < resolution; ++x) {
                const auto flat = y_flat + x;
                Particle<T> p = { .pos = cameraPos, .dir = dir };
                p.dir = vectormath::rotate(p.dir, cameraCosiney, dd_start + dd * x);
                p.dir = vectormath::rotate(p.dir, cameraCosinex, dd_start + dd * y);
                const auto r = world.intersect(p);
                if (r.valid())
                    data[flat] = r.intersection;
            }
        }

        return data;
    }
}

}