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
bool testCompound()
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

    for (const auto& d : data) {
        if (d.energy < dxmc::MAX_ENERGY<T>()) {
            dxmc::MassEnergyTransfer<T> m(d.name);
            T uen = m(d.energy);
            T diff = 1 - uen / d.nist;

            success = success && std::abs(diff) < T { 0.02 };

            if (std::abs(diff) < T { 0.1 })
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
bool testZ()
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
    std::vector<T> e(150);
    std::iota(e.begin(), e.end(), T { 1 });
    auto u = m(e);
    for (int i = 0; i < 150; ++i) {
        std::cout << e[i] << ", " << u[i] << std::endl;
    }

    for (const auto& d : data) {
        if (d.energy < dxmc::MAX_ENERGY<T>()) {
            dxmc::MassEnergyTransfer<T> m(d.Z);
            T uen = m(d.energy);
            T diff = (1 - uen / d.nist) * 100;

            success = success && std::abs(diff) < T { 2 };

            if (std::abs(diff) < T { 2 })
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

int main(int argc, char* argv[])
{
    std::cout << "Mass energy transfer tests\n";
    auto success = true;

    success = success && testZ<float>();
    success = success && testZ<double>();
    success = success && testCompound<float>();
    success = success && testCompound<double>();

    if (success)
        return EXIT_SUCCESS;
    else
        return EXIT_FAILURE;
}
