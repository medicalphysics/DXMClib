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

Copyright 2020 Erlend Andersen
*/

#pragma once

#include "dxmc/floating.h"
#include <numbers>

namespace dxmc {

template <Floating T>
constexpr T KEV_TO_ANGSTROM()
{
    return T { 12.398520 };
}

template <Floating T>
constexpr T ELECTRON_REST_MASS()
{
    return T { 510.9989461 };
}

template <Floating T>
constexpr T PI_VAL()
{
    return std::numbers::pi_v<T>;
}

template <Floating T>
constexpr T DEG_TO_RAD()
{
    return PI_VAL<T>() / T { 180 };
}
template <Floating T>
constexpr T RAD_TO_DEG()
{
    return T { 180 } / PI_VAL<T>();
}

template <Floating T>
constexpr T KEV_TO_MJ()
{
    return T { 1.6021773e-13 }; // milli Joules}
}
template <Floating T>
constexpr T MJ_TO_KEV()
{
    return T { 1 } / KEV_TO_MJ<T>();
}
}