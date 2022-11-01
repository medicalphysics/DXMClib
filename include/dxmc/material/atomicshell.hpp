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

#include <concepts>
#include <utility>
#include <vector>

namespace dxmc {
template <std::floating_point T>
struct AtomicShell {
    std::uint64_t shell = 0;
    T numberOfElectrons = 0;
    T bindingEnergy = 0;
    T HartreeFockOrbital_0 = 0;
    T numberOfPhotonsPerInitVacancy = 0;
    T energyOfPhotonsPerInitVacancy = 0;
    std::vector<std::pair<T, T>> photoel;
};

}