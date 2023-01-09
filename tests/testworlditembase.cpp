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

#include "dxmc/world/ctdiphantom.hpp"

#include <chrono>
#include <fstream>
#include <iostream>

template <typename T>
bool testIntersectAABB()
{

    std::array<T, 6> aabb { -1, -1, -1, 1, 1, 1 };
    std::array<T, 3> pos { 5, 5, 3 };
    std::array<T, 3> dir = dxmc::vectormath::scale(pos, T { -1 });
    dxmc::vectormath::normalize(dir);
    dxmc::Particle<T> p;
    p.pos = pos;
    p.dir = dir;

    auto inter = dxmc::WorldItemBase<T>::intersectAABB(p, aabb);
    bool success = inter ? true : false;

    p.pos = dxmc::vectormath::scale(p.pos, T { 0 });
    auto inter2 = dxmc::WorldItemBase<T>::intersectAABB<0>(p, aabb);
    success = success && (inter2 ? true : false);

    return success;
}

template <dxmc::Floating T>
bool testCTDIPhantom()
{
    std::array<T, 3> pos { 16, 16, -7.5 };
    std::array<T, 3> dir { 0, 0, 1 };
    dxmc::vectormath::normalize(dir);
    dxmc::CTDIPhantom<T> phantom(8, pos);
    dxmc::Particle<T> p { .pos = pos, .dir = dir, .energy = 0, .weight = 0 };
    p.pos = { 16, 16, -30 };
    auto t = phantom.intersect(p);
    bool valid = true;
    if (t.item)
        valid = valid && std::abs(t.intersection - T { 15 }) <= std::numeric_limits<T>::epsilon();
    p.pos = { 0, 0, T { -7.5 } };
    p.dir = { 1, 1, 0 };
    dxmc::vectormath::normalize(p.dir);
    t = phantom.intersect(p);

    if (t.item)
        valid = valid && std::abs(t.intersection - std::sqrt(T { 16 } * 16 + 16 * 16) + 8) <= std::numeric_limits<T>::epsilon() * 100;
    else
        valid = false;

    dxmc::CTDIPhantom<T> phantom2(9, pos);
    t.item = &phantom2;

    return valid;
}

template <dxmc::Floating T>
bool testItemCollection()
{
    return false;
}

int main(int argc, char* argv[])
{
    auto success = testCTDIPhantom<float>();
    success = success && testCTDIPhantom<double>();

    // success = success && testItemCollection<float>();
    // success = success && testItemCollection<double>();

    success = success && testIntersectAABB<float>();
    success = success && testIntersectAABB<double>();

    if (success)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}
