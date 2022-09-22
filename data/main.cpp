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

#include <string>

#include <iostream>

int main()
{
    const std::string eadl = EADLPATH;
    const std::string epdl = EPDLPATH;

    EPICSparser parser(epdl);

    for (const auto& [key, value] : parser.getElements()) {

        std::cout << "Z: " << static_cast<int>(key) << "  ";

        std::cout << "N_photo: " << value.photoelectricData().size() << "  ";
        std::cout << "N_inco: " << value.incoherentData().size() << "  ";
        std::cout << "N_coher: " << value.coherentData().size() << "  ";
        std::cout << std::endl;
    }

    return 1;
}
