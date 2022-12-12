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

#include "dxmc/world/ctdiphantom.hpp"
#include "dxmc/world/statickdtree.hpp"
#include "dxmc/world/triangle.hpp"

template <typename T>
bool teststatickdtree()
{
    bool success = true;

    dxmc::KDTree<T, dxmc::Triangle<T>, dxmc::CTDIPhantom<T>> tree;

    std::array<T, 9> vertices = { 0, 0, 0, 1, 0, 0, 0, 0, 1 };
    dxmc::Triangle<T> tri(vertices.data());

    dxmc::CTDIPhantom<T> ph1;
    dxmc::CTDIPhantom<T> ph2(8, { 20, 20, 20 }, 8);

    // adding objects
    tree.insert(tri);
    tree.insert(ph1);
    tree.insert(ph2);

    // translate all objects
    std::array<T, 3> trans { 1, 1, 1 };
    tree.translate(trans);

    // calculate aabb
    auto aabb = tree.AABB();

    return false;
}

int main()
{

    auto success = teststatickdtree<float>();
    success = success && teststatickdtree<double>();
    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}