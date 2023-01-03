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


#include "dxmc/world/world.hpp"
#include "dxmc/vectormath.hpp"

#include <chrono>
#include <fstream>
#include <iostream>

template <dxmc::Floating T>
bool testWorld()
{


}

int main(int argc, char* argv[])
{
    auto success = true;
    success = success && testWorld<double>();
    success = success && testWorld<float>();

    if (success)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}
