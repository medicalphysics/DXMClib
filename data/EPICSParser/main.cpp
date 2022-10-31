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

#include "atomicshell.hpp"
#include "epicsparser.hpp"

#include <fstream>
#include <iostream>
#include <numbers>
#include <string>
#include <filesystem>

bool test_serializer(EPICSparser& parser)
{
    auto data = parser.serializeElements();
    EPICSparser parser2(data);

    auto& elements1 = parser.getElements();
    auto& elements2 = parser2.getElements();

    if (elements1.size() != elements2.size()) {
        return false;
    }

    for (auto& [key, el1] : elements1) {
        auto& el2 = elements2.at(key);
        auto valid = el1 == el2;
        if (!valid)
            return false;
    }
    return true;
}

int main()
{
    const std::string eadl(EADLPATH);
    const std::string epdl(EPDLPATH);
    const std::string outpath(PHYSICSLISTSPATH);

    auto file_exists = std::filesystem::exists("helloworld.txt");
    if (!file_exists) {
        EPICSparser parser(eadl);
        parser.read(epdl);

        auto data = parser.serializeElements();

        std::ofstream of;
        of.open(outpath, std::ios::binary);
        of.write(data.data(), data.size());
        of.close();
    }
    return EXIT_SUCCESS;
}
