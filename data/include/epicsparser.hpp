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

#pragma once

#include "atomicelement.hpp"

#include <map>
#include <string>

class EPICSparser {
public:
    EPICSparser(const std::string& path);
    void read(const std::string& path);

    std::map<std::uint8_t, AtomicElement>& getElements() { return m_elements; }

    bool writeMaterialHeaderFile(const std::string& filename) const; 

protected:
    constexpr static std::size_t endIdx() { return 71; }

private:
    std::map<std::uint8_t, AtomicElement> m_elements;
};