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

#include "dxmc/material/material.hpp"
#include "dxmc/material/nistmaterials.hpp"

#include <iostream>
#include <string>
#include <vector>

bool testMaterials()
{
    for (std::size_t i = 1; i < 101; ++i) {
        auto mat_opt = dxmc::Material<>::byZ(i);
        if (mat_opt) {
            auto& mat = mat_opt.value();
        } else
            return false;
    }
    return true;
}

bool testNistMaterials()
{
    for (const auto& name : dxmc::NISTMaterials::listNames()) {
        auto mat_opt = dxmc::Material<>::byNistName(name);
        if (mat_opt) {
            auto& mat = mat_opt.value();
        } else
            return false;
    }
    return true;
}

int main(int argc, char* argv[])
{
    std::cout << "Basic tests of materials, please run without --fast_math flags: ";

    auto stxt = [](bool v) -> std::string { return v ? " SUCCSESS " : " FAILED "; };

    bool success = true;
    success = success && testMaterials();
    success = success && testNistMaterials();
    std::cout << stxt(success) << std::endl;

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
