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

#include "dxmc/floating.hpp"
#include "dxmc/material/atomicelement.hpp"
#include "dxmc/material/atomicshell.hpp"

#include <map>
#include <vector>

namespace dxmc {

template <Floating T>
class AtomHandler {

    static std::vector<char> dataToBinary(const std::map<std::uint64_t, AtomicElement<double>>& elements) {

    }

    protected:
    static fromBuffer(const std::vector<char>& buffer) 
    {
    }



private:
    std::map<std::uint64_t, AtomicElement> m_elements;
};

}