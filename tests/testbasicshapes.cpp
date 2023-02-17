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
    auto t1 = dxmc::basicshape::sphere::intersect(p, center, radii).value();
    success = success && t1[0] == 7998 && t1[1] == 8002;

    p.pos = { 0, 0, 1 };
    p.dir = { 0, 0, 1 };
    auto t2 = dxmc::basicshape::sphere::intersect(p, center, radii).value();
    success = success && t2[0] == -3 && t2[1] == 1;

    p.pos = { 0, 0, 1 };
    p.dir = { 0, 0, -1 };
    auto t3 = dxmc::basicshape::sphere::intersect(p, center, radii).value();
    success = success && t3[0] == -1 && t3[1] == 3;

    p.pos = { 0, 0, -1 };
    p.dir = { 0, 0, 1 };
    auto t4 = dxmc::basicshape::sphere::intersect(p, center, radii).value();
    success = success && t4[0] == -1 && t4[1] == 3;

    p.pos = { 0, 0, -1 };
    p.dir = { 0, 0, -1 };
    auto t5 = dxmc::basicshape::sphere::intersect(p, center, radii).value();
    success = success && t5[0] == -3 && t5[1] == 1;

    p.pos = { 0, 0, 8000 };
    p.dir = { 0, 0, -1 };
    auto t6 = dxmc::basicshape::sphere::intersect(p, center, radii).value();
    success = success && t6[0] == 7998 && t6[1] == 8002;

    p.pos = { 0, 0, 8000 };
    p.dir = { 0, 0, 1 };
    auto t7 = dxmc::basicshape::sphere::intersect(p, center, radii);
    success = success && !t7;

    return success;
}

template <typename T>
bool testSphereIntersectionForward()
{
    bool success = true;
    dxmc::Particle<T> p;

    T radii = 2;
    std::array<T, 3> center = { 0, 0, 0 };

    // Large value (triggers catastrophic cancelations for floats for naive solution)
    p.pos = { 0, 0, -8000 };
    p.dir = { 0, 0, 1 };
    auto t1 = dxmc::basicshape::sphere::intersectForward(p, center, radii);
    success = success && t1.value() == 7998;

    p.pos = { 0, 0, 1 };
    p.dir = { 0, 0, 1 };
    auto t2 = dxmc::basicshape::sphere::intersectForward(p, center, radii);
    success = success && t2.value() == 1;

    p.pos = { 0, 0, 1 };
    p.dir = { 0, 0, -1 };
    auto t3 = dxmc::basicshape::sphere::intersectForward(p, center, radii);
    success = success && t3.value() == 3;

    p.pos = { 0, 0, -1 };
    p.dir = { 0, 0, 1 };
    auto t4 = dxmc::basicshape::sphere::intersectForward(p, center, radii);
    success = success && t4.value() == 3;

    p.pos = { 0, 0, -1 };
    p.dir = { 0, 0, -1 };
    auto t5 = dxmc::basicshape::sphere::intersectForward(p, center, radii);
    success = success && t5.value() == 1;

    p.pos = { 0, 0, 8000 };
    p.dir = { 0, 0, -1 };
    auto t6 = dxmc::basicshape::sphere::intersectForward(p, center, radii);
    success = success && t6.value() == 7998;
    return success;
}
template <typename T>
bool testAABBIntersection()
{
    bool success = true;
    dxmc::Particle<T> p;

    std::array<T, 6> aabb = { -2, -2, -2, 2, 2, 2 };

    // Large value (triggers catastrophic cancelations for floats for naive solution)
    p.pos = { 0, 0, -8000 };
    p.dir = { 0, 0, 1 };
    auto t1 = dxmc::basicshape::AABB::intersect(p, aabb).value();
    success = success && t1[0] == 7998 && t1[1] == 8002;

    p.pos = { 0, 0, 1 };
    p.dir = { 0, 0, 1 };
    auto t2 = dxmc::basicshape::AABB::intersect(p, aabb).value();
    success = success && t2[0] == -3 && t2[1] == 1;

    p.pos = { 0, 0, 1 };
    p.dir = { 0, 0, -1 };
    auto t3 = dxmc::basicshape::AABB::intersect(p, aabb).value();
    success = success && t3[0] == -1 && t3[1] == 3;

    p.pos = { 0, 0, -1 };
    p.dir = { 0, 0, 1 };
    auto t4 = dxmc::basicshape::AABB::intersect(p, aabb).value();
    success = success && t4[0] == -1 && t4[1] == 3;

    p.pos = { 0, 0, -1 };
    p.dir = { 0, 0, -1 };
    auto t5 = dxmc::basicshape::AABB::intersect(p, aabb).value();
    success = success && t5[0] == -3 && t5[1] == 1;

    p.pos = { 0, 0, 8000 };
    p.dir = { 0, 0, -1 };
    auto t6 = dxmc::basicshape::AABB::intersect(p, aabb).value();
    success = success && t6[0] == 7998 && t6[1] == 8002;

    p.pos = { 0, 0, 8000 };
    p.dir = { 0, 0, 1 };
    auto t7 = dxmc::basicshape::AABB::intersect(p, aabb);
    success = success && !t7;
    return success;
}

template <typename T>
bool testAABBIntersectionForward()
{
    bool success = true;
    dxmc::Particle<T> p;

    std::array<T, 6> aabb = { -2, -2, -2, 2, 2, 2 };

    // Large value (triggers catastrophic cancelations for floats for naive solution)
    p.pos = { 0, 0, -8000 };
    p.dir = { 0, 0, 1 };
    auto t1 = dxmc::basicshape::AABB::intersectForward(p, aabb).value();
    success = success && t1 == 7998;

    p.pos = { 0, 0, 1 };
    p.dir = { 0, 0, 1 };
    auto t2 = dxmc::basicshape::AABB::intersectForward(p, aabb).value();
    success = success && t2 == 1;

    p.pos = { 0, 0, 1 };
    p.dir = { 0, 0, -1 };
    auto t3 = dxmc::basicshape::AABB::intersectForward(p, aabb).value();
    success = success && t3 == 3;

    p.pos = { 0, 0, -1 };
    p.dir = { 0, 0, 1 };
    auto t4 = dxmc::basicshape::AABB::intersectForward(p, aabb).value();
    success = success && t4 == 3;

    p.pos = { 0, 0, -1 };
    p.dir = { 0, 0, -1 };
    auto t5 = dxmc::basicshape::AABB::intersectForward(p, aabb).value();
    success = success && t5 == 1;

    p.pos = { 0, 0, 8000 };
    p.dir = { 0, 0, -1 };
    auto t6 = dxmc::basicshape::AABB::intersectForward(p, aabb).value();
    success = success && t6 == 7998;

    p.pos = { 0, 0, 8000 };
    p.dir = { 0, 0, 1 };
    auto t7 = dxmc::basicshape::AABB::intersectForward(p, aabb);
    success = success && !t7;
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
    auto t1 = dxmc::basicshape::cylinder::intersect(p, center, radii, half_height).value();
    success = success && t1[0] == 7998 && t1[1] == 8002;

    p.pos = { -8000, 0, 0 };
    p.dir = { 1, 0, 0 };
    auto t2 = dxmc::basicshape::cylinder::intersect(p, center, radii, half_height).value();
    success = success && t2[0] == 7998 && t2[1] == 8002;

    p.pos = { 0, 8000, 0 };
    p.dir = { 0, -1, 0 };
    auto t3 = dxmc::basicshape::cylinder::intersect(p, center, radii, half_height).value();
    success = success && t3[0] == 7998 && t3[1] == 8002;

    p.pos = { -1, 0, center[0] - half_height - 1 };
    p.dir = { 1, 0, 1 };
    dxmc::vectormath::normalize(p.dir);
    auto t4 = dxmc::basicshape::cylinder::intersect(p, center, radii, half_height).value();
    success = success && std::abs(t4[0] - std::sqrt(T { 2 })) <= std::numeric_limits<T>::epsilon() * 10;
    success = success && std::abs(t4[1] - 3 * std::sqrt(T { 2 })) <= std::numeric_limits<T>::epsilon() * 10;

    p.pos = { -1, 0, center[0] - half_height + 1 };
    p.dir = { 1, 0, 1 };
    dxmc::vectormath::normalize(p.dir);
    auto t5 = dxmc::basicshape::cylinder::intersect(p, center, radii, half_height).value();
    success = success && std::abs(t5[0] + std::sqrt(T { 2 })) <= std::numeric_limits<T>::epsilon() * 10;
    success = success && std::abs(t5[1] - 3 * std::sqrt(T { 2 })) <= std::numeric_limits<T>::epsilon() * 10;

    p.pos = { -1, 0, center[0] - half_height - 1 };
    p.dir = { -1, 0, -1 };
    dxmc::vectormath::normalize(p.dir);
    auto t6 = dxmc::basicshape::cylinder::intersect(p, center, radii, half_height);
    success = success && !t6;

    p.pos = { -1, -1, 0 };
    p.dir = { 0, 0, -1 };
    dxmc::vectormath::normalize(p.dir);
    auto t7 = dxmc::basicshape::cylinder::intersect(p, center, radii, half_height).value();
    success = success && t7[0] == -half_height;
    success = success && t7[1] == half_height;

    p.pos = { -1, -1, -half_height * 2 };
    p.dir = { 0, 0, -1 };
    dxmc::vectormath::normalize(p.dir);
    auto t8 = dxmc::basicshape::cylinder::intersect(p, center, radii, half_height);
    success = success && !t8;

    return success;
}
template <typename T>
bool testCylindarIntersectionForward()
{
    bool success = true;
    dxmc::Particle<T> p;

    const T radii = 2;
    const std::array<T, 3> center = { 0, 0, 0 };
    const T half_height = 2;

    // Large value (triggers catastrophic cancelations for floats for naive solution)
    p.pos = { 0, 0, -8000 };
    p.dir = { 0, 0, 1 };
    auto t1 = dxmc::basicshape::cylinder::intersectForward(p, center, radii, half_height).value();
    success = success && t1 == 7998;

    p.pos = { -8000, 0, 0 };
    p.dir = { 1, 0, 0 };
    auto t2 = dxmc::basicshape::cylinder::intersectForward(p, center, radii, half_height).value();
    success = success && t2 == 7998;

    p.pos = { 0, 8000, 0 };
    p.dir = { 0, -1, 0 };
    auto t3 = dxmc::basicshape::cylinder::intersectForward(p, center, radii, half_height).value();
    success = success && t3 == 7998;

    p.pos = { -1, 0, center[0] - half_height - 1 };
    p.dir = { 1, 0, 1 };
    dxmc::vectormath::normalize(p.dir);
    auto t4 = dxmc::basicshape::cylinder::intersectForward(p, center, radii, half_height).value();
    success = success && std::abs(t4 - std::sqrt(T { 2 })) <= std::numeric_limits<T>::epsilon() * 10;

    p.pos = { -1, 0, center[0] - half_height + 1 };
    p.dir = { 1, 0, 1 };
    dxmc::vectormath::normalize(p.dir);
    auto t5 = dxmc::basicshape::cylinder::intersectForward(p, center, radii, half_height).value();
    success = success && std::abs(t5 - 3 * std::sqrt(T { 2 })) <= std::numeric_limits<T>::epsilon() * 10;

    p.pos = { -1, 0, center[0] - half_height - 1 };
    p.dir = { -1, 0, -1 };
    dxmc::vectormath::normalize(p.dir);
    auto t6 = dxmc::basicshape::cylinder::intersectForward(p, center, radii, half_height);
    success = success && !t6;

    p.pos = { -1, -1, 0 };
    p.dir = { 0, 0, -1 };
    dxmc::vectormath::normalize(p.dir);
    auto t7 = dxmc::basicshape::cylinder::intersectForward(p, center, radii, half_height).value();
    success = success && t7 == half_height;

    p.pos = { -1, -1, -half_height * 2 };
    p.dir = { 0, 0, -1 };
    dxmc::vectormath::normalize(p.dir);
    auto t8 = dxmc::basicshape::cylinder::intersectForward(p, center, radii, half_height);
    success = success && !t8;

    return success;
}

int main()
{
    bool success = true;
    success = success && testSphereIntersection<float>();
    success = success && testSphereIntersection<double>();
    success = success && testSphereIntersectionForward<float>();
    success = success && testSphereIntersectionForward<double>();
    success = success && testAABBIntersection<float>();
    success = success && testAABBIntersection<double>();
    success = success && testAABBIntersectionForward<float>();
    success = success && testAABBIntersectionForward<double>();
    success = success && testCylindarIntersection<float>();
    success = success && testCylindarIntersection<double>();
    success = success && testCylindarIntersectionForward<float>();
    success = success && testCylindarIntersectionForward<double>();
    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}