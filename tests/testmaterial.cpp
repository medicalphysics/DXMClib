
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
#include "dxmc/material/nistmaterials.hpp"

#include "xraylib.h"

#include <format>
#include <fstream>
#include <iostream>
#include <numbers>

using namespace dxmc;

template <typename T>
void writeAtomTestData(std::size_t Z)
{
    std::ofstream file("atomTestData.csv");
    file << Z << std::endl;
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

        file << std::format("{},{},{},{}", e, CS_Photo(Z, e, nullptr), "photoelectric", "xlib") << std::endl;
        file << std::format("{},{},{},{}", e, CS_Rayl(Z, e, nullptr), "coherent", "xlib") << std::endl;
        file << std::format("{},{},{},{}", e, CS_Compt(Z, e, nullptr), "incoherent", "xlib") << std::endl;

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
        file << std::format("{},{},{},{}", e, CS_Total(Z, e, nullptr), "total", "xlib") << std::endl;
        file << std::format("{},{},{},{}", e, att.sum(), "total", "dxmc") << std::endl;
        file << std::format("{},{},{},{}", e, att.sum(), "total", "dxmc") << std::endl;
        file << std::format("{},{},{},{}", e, att.sum() - tot, "total", "diff") << std::endl;
        file << std::format("{},{},{},{}", e, (att.sum() / tot - 1) * 100, "total", "diffp") << std::endl;

        auto x = m.momentumTransferMax(e);
        auto ff = dxmc::interpolate(a.formFactor, x);
        file << std::format("{},{},{},{}", e, ff, "formfactor", "lin") << std::endl;
        auto ff_dx = m.formFactor(x);
        file << std::format("{},{},{},{}", e, ff_dx, "formfactor", "dxmc") << std::endl;
        T ff_xlib = FF_Rayl(Z, x, nullptr);
        file << std::format("{},{},{},{}", e, ff_xlib, "formfactor", "xlib") << std::endl;
        file << std::format("{},{},{},{}", e, ff_dx - ff, "formfactor", "diff") << std::endl;
        file << std::format("{},{},{},{}", e, (ff_dx / ff - 1) * 100, "formfactor", "diffp") << std::endl;

        auto sf = dxmc::interpolate(a.incoherentSF, x);
        file << std::format("{},{},{},{}", e, sf, "scatterfactor", "lin") << std::endl;
        auto sf_dx = m.scatterFactor(x);
        file << std::format("{},{},{},{}", e, sf_dx, "scatterfactor", "dxmc") << std::endl;
        T sf_xlib = SF_Compt(Z, x, nullptr);
        file << std::format("{},{},{},{}", e, sf_xlib, "scatterfactor", "xlib") << std::endl;
        file << std::format("{},{},{},{}", e, sf_dx - sf, "scatterfactor", "diff") << std::endl;
        file << std::format("{},{},{},{}", e, (sf_dx / sf - 1) * 100, "scatterfactor", "diffp") << std::endl;
    }
    file.close();
}

template <typename T, int N = 12>
void writeCompoundTestData(const std::string& name)
{
    std::ofstream file("compTestData.csv");
    file << name << std::endl;
    file << "e,att,type,kind" << std ::endl;

    std::vector<T> energy;
    const T step = 0.5;
    std::size_t it = 0;
    do {
        energy.push_back(step * it++ + dxmc::MIN_ENERGY<T>());
    } while (energy.back() < dxmc::MAX_ENERGY<T>());

    const auto m = dxmc::Material2<T, N>::byNistName(name).value();

    const auto Z = name.c_str();

    for (auto e : energy) {
        auto att = m.attenuationValues(e);
        file << std::format("{},{},{},{}", e, att.photoelectric, "photoelectric", "dxmc") << std::endl;
        file << std::format("{},{},{},{}", e, att.coherent, "coherent", "dxmc") << std::endl;
        file << std::format("{},{},{},{}", e, att.incoherent, "incoherent", "dxmc") << std::endl;

        file << std::format("{},{},{},{}", e, CS_Photo_CP(Z, e, nullptr), "photoelectric", "xlib") << std::endl;
        file << std::format("{},{},{},{}", e, CS_Rayl_CP(Z, e, nullptr), "coherent", "xlib") << std::endl;
        file << std::format("{},{},{},{}", e, CS_Compt_CP(Z, e, nullptr), "incoherent", "xlib") << std::endl;

        file << std::format("{},{},{},{}", e, CS_Total_CP(Z, e, nullptr), "total", "xlib") << std::endl;
        file << std::format("{},{},{},{}", e, att.sum(), "total", "dxmc") << std::endl;
        file << std::format("{},{},{},{}", e, att.sum(), "total", "dxmc") << std::endl;
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
    for (std::size_t Z = 1; Z < 85; ++Z) {
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

template <typename T>
bool testCompoundAttenuation()
{
    const auto emin = dxmc::MIN_ENERGY<T>() + 1;
    const auto emax = dxmc::MAX_ENERGY<T>();

    std::vector<T> earr(static_cast<std::size_t>(emax - emin));
    std::iota(earr.begin(), earr.end(), emin);
    bool valid = true;
    constexpr T lim = 0.2;

    const auto names = dxmc::Material2<T>::listNistCompoundNames();

    for (const auto& name : names) {

        const auto mat_opt = dxmc::Material2<T>::byNistName(name);
        const auto& comp = NISTMaterials<T>::Composition(name);
        const auto& material = mat_opt.value();
        const auto n_c = name.c_str();

        for (std::size_t i = 0; i < earr.size(); ++i) {
            const auto e = earr[i];
            const auto att = material.attenuationValues(e);
            std::array<T, 3> xlib = {
                static_cast<T>(CS_Photo_CP(n_c, e, nullptr)),
                static_cast<T>(CS_Rayl_CP(n_c, e, nullptr)),
                static_cast<T>(CS_Compt_CP(n_c, e, nullptr)),
            };

            valid = valid && std::abs(xlib[0] / att.photoelectric - 1) < lim;
            // valid = valid && std::abs(xlib[1] / att.coherent - 1) < lim;
            valid = valid && std::abs(xlib[2] / att.incoherent - 1) < lim;

            if (!valid) {
                writeCompoundTestData<T>(name);
                return false;
            }
        }
    }

    return valid;
}

int main(int argc, char* argv[])
{

    auto success = true;

    success = success && testCompoundAttenuation<float>();
    success = success && testCompoundAttenuation<double>();
    success = success && testAtomAttenuation<double>();
    success = success && testAtomAttenuation<float>();

    // writeAtomTestData<double>(13);

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
