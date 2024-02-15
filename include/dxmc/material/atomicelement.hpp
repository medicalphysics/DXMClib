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

#include "atomicshell.hpp"

#include <concepts>
#include <cstdint>
#include <map>
#include <vector>


namespace dxmc {

template <std::floating_point T = double>
struct AtomicElement {
    std::uint64_t Z = 0;
    T atomicWeight = 0;
    T standardDensity = 0;
    std::vector<std::pair<T, T>> coherent;
    std::vector<std::pair<T, T>> incoherent;
    std::vector<std::pair<T, T>> photoel;
    std::vector<std::pair<T, T>> formFactor;
    std::vector<std::pair<T, T>> incoherentSF;
    std::vector<std::pair<T, T>> incoherentMeanScatterEnergy;
    std::map<std::uint64_t, AtomicShell<T>> shells;
};

}