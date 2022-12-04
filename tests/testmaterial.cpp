
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
void writeAtomData(std::size_t Z)
{
    std::ofstream file("material.txt");
    file << "e,att,type,kind" << std ::endl;

    std::vector<T> energy;
    const T step = 0.1;
    std::size_t it = 0;
    do {
        energy.push_back(step * it++ + dxmc::MIN_ENERGY<T>());
    } while (energy.back() < dxmc::MAX_ENERGY<T>());

    auto m = dxmc::Material2<T>::byZ(Z).value();
    auto a = dxmc::AtomHandler<T>::Atom(Z);

    for (auto e : energy) {
        auto att = m.attenuationValues(e);
        file << std::format("{},{},{},{}", e, att.photoelectric, "photoelectric", "dxmc") << std::endl;
        file << std::format("{},{},{},{}", e, att.coherent, "coherent", "dxmc") << std::endl;
        file << std::format("{},{},{},{}", e, att.incoherent, "incoherent", "dxmc") << std::endl;
        auto p = dxmc::interpolate(a.photoel, e);
        file << std::format("{},{},{},{}", e, p, "photoelectric", "lin") << std::endl;
        auto co = dxmc::interpolate(a.coherent, e);
        file << std::format("{},{},{},{}", e, co, "coherent", "lin") << std::endl;
        auto inco = dxmc::interpolate(a.incoherent, e);
        file << std::format("{},{},{},{}", e, inco, "incoherent", "lin") << std::endl;

        file << std::format("{},{},{},{}", e, att.incoherent - inco, "incoherent", "diff") << std::endl;
        file << std::format("{},{},{},{}", e, att.coherent - co, "coherent", "diff") << std::endl;
        file << std::format("{},{},{},{}", e, att.photoelectric - p, "photoelectric", "diff") << std::endl;

        file << std::format("{},{},{},{}", e, (att.incoherent / inco - 1) * 100, "incoherent", "diffp") << std::endl;
        file << std::format("{},{},{},{}", e, (att.coherent / co - 1) * 100, "coherent", "diffp") << std::endl;
        file << std::format("{},{},{},{}", e, (att.photoelectric / p - 1) * 100, "photoelectric", "diffp") << std::endl;

        auto tot = p + inco + co;
        file << std::format("{},{},{},{}", e, tot, "total", "lin") << std::endl;
        file << std::format("{},{},{},{}", e, att.sum(), "total", "dxmc") << std::endl;
        file << std::format("{},{},{},{}", e, att.sum() - tot, "total", "diff") << std::endl;
        file << std::format("{},{},{},{}", e, (att.sum() / tot - 1) * 100, "total", "diffp") << std::endl;
    }
    file.close();
}
template <dxmc::Floating T = double>
bool testAtoms()
{
    const auto emin = dxmc::MIN_ENERGY<T>();
    const auto emax = dxmc::MAX_ENERGY<T>();

    std::vector<T> earr(static_cast<std::size_t>(emax - emin));
    std::iota(earr.begin(), earr.end(), emin);
    bool valid = true;
    constexpr T lim = 1;
    for (std::size_t Z = 11; Z < 85; ++Z) {
        const auto atom = dxmc::AtomHandler<T>::Atom(Z);
        const auto material = dxmc::Material2<T>::byZ(Z).value();
        for (const auto& e : earr) {
            auto att = material.attenuationValues(e);
            auto photo_val = dxmc::interpolate(atom.photoel, e);

            // valid = valid && (std::abs(att.photoelectric - photo_val) / photo_val * 100) < lim;
            auto lne = std::log(e);

            auto coher_val = dxmc::interpolate(atom.coherent, e);
            if (e > T { 8 }) {
                //  valid = valid && (std::abs(att.coherent - coher_val) / coher_val * 100) < lim;
            }

            auto incoher_val = dxmc::interpolate(atom.incoherent, e);
            // valid = valid && (std::abs(att.incoherent - incoher_val) / incoher_val * 100) < lim;

            auto att_sum = photo_val + incoher_val + coher_val;
            auto att_sum_dxmc = att.sum();
            valid = valid && (std::abs(att.sum() - att_sum) / att_sum * 100) < lim;

            if (!valid) {
                writeAtomData<T>(Z);
                return valid;
            }
        }
    }

    return valid;
}

void testMaterial()
{
    auto m0 = dxmc::Material2<double>::byChemicalFormula("He");
    // auto m0 = dxmc::Material2<double>::byZ(2);
    if (m0) {
        auto m = m0.value();
        std::vector<double> e(150);
        std::iota(e.begin(), e.end(), 1.0);
    }

    auto m1 = dxmc::Material2<double>::byChemicalFormula("Ca5(PO4)3");
    if (m1) {
        auto m = m1.value();
        auto att = m.attenuationValues(60.0);
        auto sum = att.coherent + att.incoherent + att.photoelectric;
        std::cout << sum << std::endl;

        std::vector<double> e(150);
        std::iota(e.begin(), e.end(), 1.0);
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
        auto sum1 = att1.sum();
        auto att2 = m.attenuationValues(4.04);
        auto sum2 = att2.sum();
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
    auto valid = testAtoms();
    // testMaterial();
    // testInterpolator();
    auto success = valid;
    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
