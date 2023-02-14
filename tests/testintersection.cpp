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

#include "dxmc/world/intersection.hpp"

#include <iostream>

template <typename T>
bool testSphereIntersection()
{
    bool success = true;
    dxmc::Particle<T> p;
    p.pos = { 0, 0, -8000 };
    p.dir = { 0, 0, 1 };

    T radii = 1;
    std::array<T, 3> center = { 0, 0, 0 };
    auto t1 = intersectSphere(p, center, radii);
    if (t1)
        success = success && true;

    p.pos = { 0, 0, 0 };
    auto t2 = intersectSphere(p, center, radii);
    if (t2)
        success = success && true;

    p.pos = { 0, 0, -8000 };
    auto t3 = intersectSphereForward(p, center, radii);
    if (t3)
        success = success && true;

    p.pos = { 0, 0, 0 };
    auto t4 = intersectSphereForward(p, center, radii);
    if (t4)
        success = success && true;

    return success;
}

int main()
{
    bool success = true;
    success = success && testSphereIntersection<float>();
    success = success && testSphereIntersection<double>();

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}