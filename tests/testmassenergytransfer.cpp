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

#include "dxmc/constants.hpp"
#include "dxmc/material/massenergytransfer.hpp"
#include "dxmc/material/material.hpp"

#include <iostream>
#include <string>
#include <vector>

template <typename T>
bool testCompoundNIST()
{
    bool success = true;

    struct data_t {
        std::string name;
        T energy = 0;
        T nist = 0;
    };

    std::vector<data_t> data;
    data.push_back({ .name = "Air, Dry (near sea level)", .energy = 10, .nist = 4.742 });
    data.push_back({ .name = "Air, Dry (near sea level)", .energy = 60, .nist = 3.041E-02 });
    data.push_back({ .name = "Air, Dry (near sea level)", .energy = 150, .nist = 2.496E-02 });

    data.push_back({ .name = "Polymethyl Methacralate (Lucite, Perspex)", .energy = 10, .nist = 3.026E+00 });
    data.push_back({ .name = "Polymethyl Methacralate (Lucite, Perspex)", .energy = 60, .nist = 2.530E-02 });
    data.push_back({ .name = "Polymethyl Methacralate (Lucite, Perspex)", .energy = 150, .nist = 2.755E-02 });

    for (const auto& d : data) {
        if (d.energy < dxmc::MAX_ENERGY<T>()) {
            dxmc::MassEnergyTransfer<T> m(d.name);
            T uen = m(d.energy);
            T diff = (uen / d.nist - 1) * 100;
            bool valid = std::abs(diff) < 2 || std::abs(uen - d.nist) < T { 0.001 };
            success = success && valid;

            if (valid)
                std::cout << "SUCCESS ";
            else
                std::cout << "FAILURE ";
            std::cout << "Mass energy absorbtion for " << d.name << ", energy = " << d.energy << " [keV] sizeof(T) = " << sizeof(T) << "\ndxmc: ";
            std::cout << uen;
            std::cout << " NIST: ";
            std::cout << d.nist;
            std::cout << " Difference [%]: " << diff << std::endl;
        }
    }
    return success;
}

template <typename D>
struct zdata_t {
    std::size_t Z = 0;
    D energy = 0;
    D nist = 0;
};

template <typename T>
bool testZNIST()
{
    bool success = true;

    std::vector<zdata_t<T>> data;

    data.push_back({ .Z = 82, .energy = 150, .nist = 1.056E+00 });
    data.push_back({ .Z = 82, .energy = 10, .nist = 1.247E+02 });
    data.push_back({ .Z = 82, .energy = 60, .nist = 4.149 });

    data.push_back({ .Z = 6, .energy = 10, .nist = 2.078 });
    data.push_back({ .Z = 6, .energy = 60, .nist = 2.098E-2 });
    data.push_back({ .Z = 6, .energy = 150, .nist = 2.449E-2 });

    dxmc::MassEnergyTransfer<T> m(82);

    for (const auto& d : data) {
        if (d.energy < dxmc::MAX_ENERGY<T>()) {
            dxmc::MassEnergyTransfer<T> m(d.Z);
            T uen = m(d.energy);
            T diff = (uen / d.nist - 1) * 100;
            bool valid = std::abs(diff) < 2 || std::abs(uen - d.nist) < T { 0.001 };
            success = success && valid;

            if (valid)
                std::cout << "SUCCESS ";
            else
                std::cout << "FAILURE ";
            std::cout << "Mass energy absorbtion for Z = " << d.Z << ", energy = " << d.energy << " [keV] sizeof(T) = " << sizeof(T) << "\ndxmc: ";
            std::cout << uen;
            std::cout << " NIST: ";
            std::cout << d.nist;
            std::cout << " Difference [%]: " << diff << std::endl;
        }
    }
    return success;
}

template <typename T>
bool testInterpolation()
{
    bool success = true;
    // atoms

    std::vector<T> energy;
    T e = dxmc::MIN_ENERGY<T>();
    while (e < dxmc::MAX_ENERGY<T>()) {
        energy.push_back(e);
        e += T { 0.5 };
    }

    for (const auto& [Z, atom] : dxmc::AtomHandler<T>::allAtoms()) {
        int N = 0;
        dxmc::MassEnergyTransfer<T, true> mCub(Z);
        dxmc::MassEnergyTransfer<T, false> mLin(Z);

        std::vector<T> cub, lin;

        for (const auto e : energy) {
            const auto u_lin = mLin(e);
            const auto u_cub = mCub(e);

            cub.push_back(u_cub);
            lin.push_back(u_lin);

            const T diffp = std::abs(u_cub / u_lin - 1) * 100;
            const T diffa = std::abs(u_cub - u_lin);
            bool valid = diffp < 2 || diffa < 1;

            if (!valid) {
                std::cout << Z << ": Energy: " << e << " Diff abs: " << diffa << " Diff per: " << diffp << std::endl;
                int test = 0;
                N++;
            }
        }

        const auto sqr = std::transform_reduce(cub.cbegin(), cub.cend(), lin.cbegin(), T { 0 }, std::plus<>(), [](const auto c, const auto l) { const auto d = c - l; return d*d; });
        const auto rms = std::sqrt(sqr / cub.size());

        success = success && N < 10;
    }

    return success;
}

int main(int argc, char* argv[])
{
    std::cout << "Mass energy transfer tests\n";
    auto success = true;

    // success = success && testInterpolation<float>();

    success = success && testZNIST<double>();
    success = success && testZNIST<float>();
    success = success && testCompoundNIST<float>();
    success = success && testCompoundNIST<double>();

    if (success)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}
