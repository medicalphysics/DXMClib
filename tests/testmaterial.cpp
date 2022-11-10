
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

#include <iostream>
using namespace dxmc;

void testInterpolator()
{

    auto O = AtomHandler<double>::Atom(82);

    std::vector<double> x(O.photoel.size());
    std::vector<double> y(O.photoel.size());
    std::transform(O.photoel.begin(), O.photoel.end(), x.begin(), [](const auto& val) { return std::log(val.first); });
    std::transform(O.photoel.begin(), O.photoel.end(), y.begin(), [](const auto& val) { return std::log(val.second); });

    auto intp = CubicLSInterpolator<double>(x, y, 15);

    for (auto [e, a] : O.photoel) {
        std::cout << e << ", " << a << ", " << std::exp(intp(std::log(e))) << ", " << CS_Photo(O.Z, e, nullptr) << std::endl;
    }
    return;
}

bool atomHandler()
{
    const auto& atom = AtomHandler<float>::Atom(6);
    auto test = interpolate(atom.photoel, float { 60 });
    return true;
}

bool testParser()
{

    auto test = Material2<float>::parseCompoundStr("Ca5(PO4)3");
    return true;
}

int main(int argc, char* argv[])
{
    testInterpolator();
    testParser();
    auto success = atomHandler();
    success = success && testParser();
    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
