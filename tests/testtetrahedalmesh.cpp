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

#include "dxmc/world/worlditems/tetrahedalmesh.hpp"

#include <iostream>

template <typename T>
bool testReader()
{
    dxmc::ThetrahedalMesh<T> mesh;
    mesh.readICRP145Phantom("MRCP_AM.node", "MRCP_AM.ele");
    return false;
}

int main()
{
    std::cout << "Testing ray intersection on basic shapes\n";

    bool success = true;

    success = success && testReader<double>();
    success = success && testReader<float>();

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}