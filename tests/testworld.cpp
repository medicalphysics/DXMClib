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

#include "dxmc/particle.hpp"
#include "dxmc/world/ctdiphantom.hpp"

#include <chrono>
#include <fstream>
#include <iostream>

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
    bool valid = t.has_value();
    if (t)
        valid = valid && std::abs(t.value() - T { 15 }) <= std::numeric_limits<T>::epsilon();
    p.pos = { 0, 0, T { -7.5 } };
    p.dir = { 1, 1, 0 };
    dxmc::vectormath::normalize(p.dir);
    t = phantom.intersect(p);
    valid = valid && t.has_value();
    if (t)
        valid = valid && std::abs(t.value() - std::sqrt(T { 16 } * 16 + 16 * 16) + 8) <= std::numeric_limits<T>::epsilon() * 100;
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

    success = success && testItemCollection<float>();
    success = success && testItemCollection<double>();

    if (success)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}
