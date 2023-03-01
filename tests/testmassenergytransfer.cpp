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

#include "dxmc/material/massenergytransfer.hpp"
#include "dxmc/material/material.hpp"

#include <iostream>

template <typename T>
bool testZ()
{
    dxmc::MassEnergyTransfer<T> me(6);

    auto m = dxmc::Material2<T>::byZ(6).value();

    auto test = me(T { 60 });
    auto test2 = m.attenuationValues(T { 60 });
    auto test3 = test2.sum();
    return false;
}

int main(int argc, char* argv[])
{
    std::cout << "Mass energy transfer tests\n";
    auto success = true;

    success = success && testZ<float>();
    success = success && testZ<double>();

    if (success)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}
