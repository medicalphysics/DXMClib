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

#include "epicsparser.hpp"

#include <format>
#include <fstream>
#include <iostream>
#include <numbers>
#include <string>

int main()
{
    const std::string eadl = EADLPATH;
    const std::string epdl = EPDLPATH;

    EPICSparser parser(eadl);
    parser.read(epdl);

    for (const auto& [key, value] : parser.getElements()) {

        std::cout << "Z: " << static_cast<int>(key) << "  ";
        std::cout << "N_photo: " << value.photoelectricData().size() << "  ";
        std::cout << "N_inco: " << value.incoherentData().size() << "  ";
        std::cout << "N_coher: " << value.coherentData().size() << "  ";
        std::cout << std::endl;
    }

    const std::uint8_t Z = 16;
    const double angle = std::numbers::pi;
    // write Form factor data
    std::ofstream f;

    const auto& elements = parser.getElements();
    f.open("formfactor_test.txt");
    for (const auto& [x, v] : elements.at(Z).formFactor()) {
        f << x << ", " << v << std::endl;
    }
    f.close();

    f.open("imagSF_test.txt");
    for (const auto& [x, v] : elements.at(Z).imaginaryAnomalousSF()) {
        f << x << ", " << v << std::endl;
    }
    f.close();
    f.open("realSF_test.txt");
    for (const auto& [x, v] : elements.at(Z).realAnomalousSF()) {
        f << x << ", " << v << std::endl;
    }
    f.close();
    f.open("SF_test.txt");
    for (const auto& [x, v] : elements.at(Z).incoherentSF()) {
        f << x << ", " << v << std::endl;
    }
    f.close();

    return 1;
}
