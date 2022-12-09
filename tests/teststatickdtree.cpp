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

template <typename T>
bool teststatickdtree()
{

    dxmc::CTDIPhantom<T> ph1;
    dxmc::CTDIPhantom<T> ph2(8, { 20, 20, 20 }, 8);




    auto test = dxmc::Test(T { 1 }, ph1);


    //dxmc::StaticKDTree<T, dxmc::CTDIPhantom<T>> tree(8);
    //tree.insert(ph1);
    //tree.insert(ph2);


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