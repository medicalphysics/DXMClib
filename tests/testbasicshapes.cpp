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

#include "dxmc/world/basicshapes/aabb.hpp"
#include "dxmc/world/basicshapes/cylinder.hpp"
#include "dxmc/world/basicshapes/sphere.hpp"

#include "dxmc/dxmcrandom.hpp"

#include <iostream>

template <typename T>
bool testSphereIntersection()
{
    bool success = true;
    dxmc::Particle<T> p;

    T radii = 2;
    std::array<T, 3> center = { 0, 0, 0 };

    // Large value (triggers catastrophic cancelations for floats for naive solution)
    p.pos = { 0, 0, -8000 };
    p.dir = { 0, 0, 1 };
    auto t = dxmc::basicshape::sphere::intersect(p, center, radii);
    success = success && t.valid() && t.intersection == 7998;

    p.pos = { 0, 0, 1 };
    p.dir = { 0, 0, 1 };
    t = dxmc::basicshape::sphere::intersect(p, center, radii);
    success = success && t.valid() && t.intersection == 1;

    p.pos = { 0, 0, 1 };
    p.dir = { 0, 0, -1 };
    t = dxmc::basicshape::sphere::intersect(p, center, radii);
    success = success && t.valid() && t.intersection == 3;

    p.pos = { 0, 0, -1 };
    p.dir = { 0, 0, 1 };
    t = dxmc::basicshape::sphere::intersect(p, center, radii);
    success = success && t.valid() && t.intersection == 3;

    p.pos = { 0, 0, -1 };
    p.dir = { 0, 0, -1 };
    t = dxmc::basicshape::sphere::intersect(p, center, radii);
    success = success && t.valid() && t.intersection == 1;

    p.pos = { 0, 0, 8000 };
    p.dir = { 0, 0, -1 };
    t = dxmc::basicshape::sphere::intersect(p, center, radii);
    success = success && t.valid() && t.intersection == 7998;

    p.pos = { 0, 0, 8000 };
    p.dir = { 0, 0, 1 };
    t = dxmc::basicshape::sphere::intersect(p, center, radii);
    success = success && !t.valid();

    if (success)
        std::cout << "SUCCESS";
    else
        std::cout << "FAILURE";
    std::cout << " Test for ray sphere intersections with sizeof(T): " << sizeof(T) << std::endl;

    return success;
}

template <typename T>
bool testAABBIntersection()
{
    bool success = true;
    dxmc::Particle<T> p;

    std::array<T, 6> aabb = { -2, -2, -2, 2, 2, 2 };

    p.pos = { 0, 0, -8000 };
    p.dir = { 0, 0, 1 };
    auto t = dxmc::basicshape::AABB::intersect(p, aabb);
    success = success && t.valid() && t.intersection == 7998 && !t.rayOriginIsInsideItem;

    p.pos = { 0, 0, 1 };
    p.dir = { 0, 0, 1 };
    t = dxmc::basicshape::AABB::intersect(p, aabb);
    success = success && t.valid() && t.intersection == 1 && t.rayOriginIsInsideItem;

    p.pos = { 0, 0, 1 };
    p.dir = { 0, 0, -1 };
    t = dxmc::basicshape::AABB::intersect(p, aabb);
    success = success && t.valid() && t.intersection == 3 && t.rayOriginIsInsideItem;

    p.pos = { 0, 0, -1 };
    p.dir = { 0, 0, 1 };
    t = dxmc::basicshape::AABB::intersect(p, aabb);
    success = success && t.valid() && t.intersection == 3 && t.rayOriginIsInsideItem;

    p.pos = { 0, 0, -1 };
    p.dir = { 0, 0, -1 };
    t = dxmc::basicshape::AABB::intersect(p, aabb);
    success = success && t.valid() && t.intersection == 1 && t.rayOriginIsInsideItem;

    p.pos = { 0, 0, 8000 };
    p.dir = { 0, 0, -1 };
    t = dxmc::basicshape::AABB::intersect(p, aabb);
    success = success && t.valid() && t.intersection == 7998 && !t.rayOriginIsInsideItem;

    p.pos = { 0, 0, 8000 };
    p.dir = { 0, 0, 1 };
    t = dxmc::basicshape::AABB::intersect(p, aabb);
    success = success && !t.valid();

    if (success)
        std::cout << "SUCCESS";
    else
        std::cout << "FAILURE";
    std::cout << " Test for ray AABB intersections with sizeof(T): " << sizeof(T) << std::endl;

    return success;
}

template <typename T>
bool testCylindarIntersection()
{
    bool success = true;
    dxmc::Particle<T> p;

    const T radii = 2;
    const std::array<T, 3> center = { 0, 0, 0 };
    const T half_height = 2;

    // Large value (triggers catastrophic cancelations for floats for naive solution)
    p.pos = { 0, 0, -8000 };
    p.dir = { 0, 0, 1 };
    auto t = dxmc::basicshape::cylinder::intersect(p, center, radii, half_height);
    success = success && t.valid() && t.intersection == 7998;

    p.pos = { -8000, 0, 0 };
    p.dir = { 1, 0, 0 };
    t = dxmc::basicshape::cylinder::intersect(p, center, radii, half_height);
    success = success && t.valid() && t.intersection == 7998;

    p.pos = { 0, 8000, 0 };
    p.dir = { 0, -1, 0 };
    t = dxmc::basicshape::cylinder::intersect(p, center, radii, half_height);
    success = success && t.valid() && t.intersection == 7998;

    p.pos = { -1, 0, center[0] - half_height - 1 };
    p.dir = { 1, 0, 1 };
    dxmc::vectormath::normalize(p.dir);
    t = dxmc::basicshape::cylinder::intersect(p, center, radii, half_height);
    success = success && t.valid() && std::abs(t.intersection - std::sqrt(T { 2 })) <= std::numeric_limits<T>::epsilon() * 10;

    p.pos = { -1, 0, center[0] - half_height + 1 };
    p.dir = { 1, 0, 1 };
    dxmc::vectormath::normalize(p.dir);
    t = dxmc::basicshape::cylinder::intersect(p, center, radii, half_height);
    success = success && t.valid() && std::abs(t.intersection - 3 * std::sqrt(T { 2 })) <= std::numeric_limits<T>::epsilon() * 10;

    p.pos = { -1, 0, center[0] - half_height - 1 };
    p.dir = { -1, 0, -1 };
    dxmc::vectormath::normalize(p.dir);
    t = dxmc::basicshape::cylinder::intersect(p, center, radii, half_height);
    success = success && !t.valid();

    p.pos = { -1, -1, 0 };
    p.dir = { 0, 0, -1 };
    dxmc::vectormath::normalize(p.dir);
    t = dxmc::basicshape::cylinder::intersect(p, center, radii, half_height);
    success = success && t.valid() && t.intersection == half_height;

    p.pos = { -1, -1, -half_height * 2 };
    p.dir = { 0, 0, -1 };
    dxmc::vectormath::normalize(p.dir);
    t = dxmc::basicshape::cylinder::intersect(p, center, radii, half_height);
    success = success && !t.valid();

    if (success)
        std::cout << "SUCCESS";
    else
        std::cout << "FAILURE";
    std::cout << " Test for ray cylinder intersections with sizeof(T): " << sizeof(T) << std::endl;

    return success;
}

template <typename T>
bool testCylinderForwardIntersection()
{

    bool success = true;

    const T radii = 4;
    const std::array<T, 3> center = { 8, 0, 0 };
    const T half_height = 6;

    dxmc::Particle<T> p;

    std::optional<std::array<T, 2>> res;

    std::array dummy = { std::numeric_limits<T>::max(), std::numeric_limits<T>::max() };
    std::array<T, 2> t;
    
    p.pos = { -10, 0, 0 };
    p.dir = { 1, 0, 0 };
    res = dxmc::basicshape::cylinder::intersectForwardInterval(p, center, radii, half_height);
    t = res.value_or(dummy);
    success = success && t[0] == 14;
    success = success && t[1] == 22;

    p.pos = { 6, 0, 0 };
    p.dir = { 1, 0, 0 };
    res = dxmc::basicshape::cylinder::intersectForwardInterval(p, center, radii, half_height);
    t = res.value_or(dummy);
    success = success && t[0] == 0;
    success = success && t[1] == 6;

    p.pos = { 8, 0, 0 };
    p.dir = { 0, 0, 1 };
    res = dxmc::basicshape::cylinder::intersectForwardInterval(p, center, radii, half_height);
    t = res.value_or(dummy);
    success = success && t[0] == 0;
    success = success && t[1] == 6;
    
    dxmc::RandomState state;
    for(std::size_t i =0;i < 1e6;++i) {
        for (std::size_t i = 0; i < 3; ++i) {
            p.pos[i] = state.randomUniform(center[i] - radii * 2, center[i] + radii * 2);
            p.dir[i] = state.randomUniform(T { -1 }, T { 1 });
        }
        dxmc::vectormath::normalize(p.dir);
        auto r = dxmc::basicshape::cylinder::intersectForwardInterval(p, center, radii, half_height);
        if (r) {
            auto t2 = r.value();
            if (t2[0] < T { 0 } || t2[1] < T { 0 } || t2[1] <= t2[0])
                return false;
        }
    }

    return success;
}

int main()
{
    std::cout << "Testing ray intersection on basic shapes\n";

    bool success = true;

    success = success && testCylinderForwardIntersection<double>();
    success = success && testCylinderForwardIntersection<float>();

    success = success && testSphereIntersection<float>();
    success = success && testSphereIntersection<double>();

    success = success && testAABBIntersection<float>();
    success = success && testAABBIntersection<double>();

    success = success && testCylindarIntersection<float>();
    success = success && testCylindarIntersection<double>();

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}