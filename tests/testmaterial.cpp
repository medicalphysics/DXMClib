
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
#include <fstream>
#include <iostream>
#include <numbers>
using namespace dxmc;

template <typename T>
void writeMaterialData(const dxmc::Material2<T>& m, const std::vector<T>& energy)
{
    std::ofstream file("material.txt");
    file << "e,photo,coherent,incoherent,x,formfactor,scatterfactor" << std::endl;
    constexpr auto pi = std::numbers::pi_v<T>;
    for (auto e : energy) {
        auto a = m.attenuationValues(e);
        auto x = m.momentumTransfer(e, pi);
        file << std::format("{},{},{},{},{},{},{}", e, a[0], a[1], a[2], x, m.formFactor(e, pi), m.scatterFactor(e, pi)) << std::endl;
    }

    file.close();
}

void testMaterial()
{
    //auto m0 = dxmc::Material2<double>::byChemicalFormula("Ca5(PO4)3");
    auto m0 = dxmc::Material2<double>::byChemicalFormula("Ca5(PO4)3");
    if (m0) {
        auto m = m0.value();
        std::vector<double> e(150);
        std::iota(e.begin(), e.end(), 1.0);
        writeMaterialData(m, e);
    }

    auto m1 = dxmc::Material2<double>::byChemicalFormula("Ca5(PO4)3");
    if (m1) {
        auto m = m1.value();
        auto att = m.attenuationValues(60.0);
        auto sum = std::reduce(att.cbegin(), att.cend());
        std::cout << sum << std::endl;

        std::vector<double> e(150);
        std::iota(e.begin(), e.end(), 1.0);
        writeMaterialData(m, e);

    }
    std::map<std::uint64_t, double> w;
    w[1] = 0.034000;
    w[6] = 0.155000;
    w[7] = 0.042000;
    w[8] = 0.435000;
    w[11] = 0.001000;
    w[12] = 0.002000;
    w[15] = 0.103000;
    w[16] = 0.003000;
    w[20] = 0.225000;
    auto bone = dxmc::Material2<double>::byWeight(w);
    if (bone) {
        auto m = bone.value();
        auto att1 = m.attenuationValues(4.038);
        auto sum1 = std::reduce(att1.cbegin(), att1.cend());
        auto att2 = m.attenuationValues(4.04);
        auto sum2 = std::reduce(att2.cbegin(), att2.cend());
        std::cout << sum1 << std::endl;
    }

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
