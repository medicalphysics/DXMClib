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

Copyright 2019 Erlend Andersen
*/

#pragma once

#include "dxmc/floating.hpp"

#include <array>

namespace dxmc {
/**
 * @brief Simple struct to describe a photon
 */
template <Floating T = double>
struct Particle {
    /**
     * @brief Position vector in three dimensions.
     */
    std::array<T, 3> pos;
    /**
     * @brief Direction vector in three dimension. This vector is threated as a normal vector.
     */
    std::array<T, 3> dir;
    /**
     * @brief Photon energy in keV.
     */
    T energy;
    /**
     * @brief Photon relative weight.
     */
    T weight;

    inline void translate(const T dist)
    {
        pos[0] += dir[0] * dist;
        pos[1] += dir[1] * dist;
        pos[2] += dir[2] * dist;
    }

    inline static constexpr T border_translate_minimum()
    {
        if constexpr (sizeof(T) == 4)
            return 1e-4;
        else
            return 1e-6;
    }
    inline void border_translate(const T dist)
    {
        // Make sure we translate particle beyond any border we want to translate to
        // We simply add 100 nm to the distance, works for float and double
        translate(dist + border_translate_minimum());
    }
};
}