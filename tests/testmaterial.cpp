
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

Copyright 2023 Erlend Andersen
*/

#include "dxmc/material/atomhandler.hpp"
#include "dxmc/material/atomicelement.hpp"
#include "dxmc/material/atomicshell.hpp"
#include "dxmc/material/atomserializer.hpp"
#include "dxmc/material/material.hpp"

#include <format>
#include <fstream>
#include <iostream>
#include <numbers>

using namespace dxmc;

template <typename T>
void writeAtomTestData(std::size_t Z)
{
    std::ofstream file("atomTestData.csv");
    file << "e,att,type,kind" << std ::endl;

    std::vector<T> energy;
    const T step = 0.5;
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

        auto x = m.momentumTransferMax(e);
        auto ff = dxmc::interpolate(a.formFactor, x);
        file << std::format("{},{},{},{}", e, ff, "formfactor", "lin") << std::endl;
        auto ff_dx = m.formFactor(x);
        file << std::format("{},{},{},{}", e, ff_dx, "formfactor", "dxmc") << std::endl;
        file << std::format("{},{},{},{}", e, ff_dx - ff, "formfactor", "diff") << std::endl;
        file << std::format("{},{},{},{}", e, (ff_dx / ff - 1) * 100, "formfactor", "diffp") << std::endl;

        auto sf = dxmc::interpolate(a.incoherentSF, x);
        file << std::format("{},{},{},{}", e, sf, "scatterfactor", "lin") << std::endl;
        auto sf_dx = m.scatterFactor(x);
        file << std::format("{},{},{},{}", e, sf_dx, "scatterfactor", "dxmc") << std::endl;
        file << std::format("{},{},{},{}", e, sf_dx - sf, "scatterfactor", "diff") << std::endl;
        file << std::format("{},{},{},{}", e, (sf_dx / sf - 1) * 100, "scatterfactor", "diffp") << std::endl;
    }
    file.close();
}
template <dxmc::Floating T = double>
bool testAtomAttenuation()
{
    const auto emin = dxmc::MIN_ENERGY<T>() + 1;
    const auto emax = dxmc::MAX_ENERGY<T>();

    std::vector<T> earr(static_cast<std::size_t>(emax - emin));
    std::iota(earr.begin(), earr.end(), emin);
    bool valid = true;
    constexpr T lim = 2;
    for (std::size_t Z = 13; Z < 85; ++Z) {
        const auto atom = dxmc::AtomHandler<T>::Atom(Z);
        const auto material = dxmc::Material2<T>::byZ(Z).value();
        for (const auto& e : earr) {
            auto att = material.attenuationValues(e);
            auto photo_val = dxmc::interpolate(atom.photoel, e);

            valid = valid && (std::abs(att.photoelectric - photo_val) / photo_val * 100) < lim;
            auto lne = std::log(e);

            auto coher_val = dxmc::interpolate(atom.coherent, e);
            valid = valid && (std::abs(att.coherent - coher_val) / coher_val * 100) < 100;

            auto incoher_val = dxmc::interpolate(atom.incoherent, e);
            valid = valid && (std::abs(att.incoherent - incoher_val) / incoher_val * 100) < lim;

            auto att_sum = photo_val + incoher_val + coher_val;
            auto att_sum_dxmc = att.sum();
            valid = valid && (std::abs(att.sum() - att_sum) / att_sum * 100) < lim;

            auto x = material.momentumTransferMax(e);
            auto ff_lin = dxmc::interpolate(atom.formFactor, x);
            auto ff_dx = material.formFactor(x);
            valid = valid && (std::abs(ff_dx / ff_lin) - 1) * 100 < 20;

            auto sf_lin = dxmc::interpolate(atom.incoherentSF, x);
            auto sf_dx = material.scatterFactor(x);
            valid = valid && (std::abs(sf_dx / sf_lin) - 1) * 100 < 20;

            if (!valid) {
                writeAtomTestData<T>(Z);
                return valid;
            }
        }
    }
    return valid;
}

int main(int argc, char* argv[])
{
    auto success = true;
    success = success && testAtomAttenuation();

    // writeAtomTestData<double>(13);

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
