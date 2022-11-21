
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

#include "dxmc/material/atomhandler.hpp"
#include "dxmc/material/atomicelement.hpp"
#include "dxmc/material/atomicshell.hpp"
#include "dxmc/material/atomserializer.hpp"
#include "dxmc/material/material.hpp"

#include "xraylib.h"

#include <format>
#include <iostream>
using namespace dxmc;
void testMaterial()
{
    auto m1 = dxmc::Material2<double>::byChemicalFormula("H2O");
    auto m2 = dxmc::Material2<float>::byChemicalFormula("Ca5(PO4)3");
    auto m3 = dxmc::Material2<double>::byChemicalFormula("Ca5(PO4)3");
}
void testShell()
{
}

void testInterpolator()
{

    for (std::size_t i = 1; i <= 83; i++) {
        auto O = AtomHandler<double>::Atom(i);
        std::vector<double> x(O.photoel.size());
        std::vector<double> y(O.photoel.size());
        std::transform(O.photoel.begin(), O.photoel.end(), x.begin(), [](const auto& val) { return std::log(val.first); });
        std::transform(O.photoel.begin(), O.photoel.end(), y.begin(), [](const auto& val) { return std::log(val.second); });

        auto intp = CubicLSInterpolator<double>(x, y, 20);

        double max_error_int = 0;
        double max_error_lib = 0;
        std::size_t ind = 0;
        std::size_t ind_max = 0;
        std::size_t ind_max_lib = 0;

        std::vector<double> binding_e;
        for (const auto& [sidx, shell] : O.shells) {
            binding_e.push_back(shell.bindingEnergy);
        }

        for (auto [e, a] : O.photoel) {
            auto bindingIdx = std::find(binding_e.cbegin(), binding_e.cend(), e);
            if (bindingIdx == binding_e.cend()) {
                auto ai = std::exp(intp(std::log(e)));
                auto ei = CS_Photo(O.Z, e, nullptr);
                if (std::abs(ai / a - 1) > max_error_int)
                    ind_max = ind;
                if (std::abs(ei / a - 1) > max_error_lib)
                    ind_max_lib = ind;
                max_error_int = std::max(std::abs(ai / a - 1), max_error_int);
                max_error_lib = std::max(std::abs(ei / a - 1), max_error_lib);
            }
            ind++;
        }
        auto s = std::format("{}: Max Error :{}  xlib: {}  at {} and {} of {}", O.Z, max_error_int, max_error_lib, ind_max, ind_max_lib, ind);
        std::cout << s << std::endl;
        // for (auto [e, a] : O.photoel) {
        //    std::cout << e << ", " << a << ", " << std::exp(intp(std::log(e))) << ", " << CS_Photo(O.Z, e, nullptr) << std::endl;
        //}
    }
    return;
}

int main(int argc, char* argv[])
{
    testMaterial();
    testInterpolator();
    auto success = true;
    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
