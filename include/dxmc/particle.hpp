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

#include "dxmc/constants.hpp"

#include <array>

namespace dxmc {

template <std::uint_fast32_t N = 5>
class ParticleTrack {
public:
    void registerPosition(const std::array<double, 3>& pos)
    {
        // incrementing and circle around
        m_index++;
        m_history[m_index % N] = pos;
    }

    std::array<std::array<double, 3>, N> getHistory() const
    {
        std::array<std::array<double, 3>, N> r;
        if (m_index < N) {
            for (std::size_t i = 0; i < m_index; ++i)
                r[i] = m_history[i];            
        } else {
            const std::uint_fast32_t first = (m_index + std::uint_fast32_t { 1 }) % N;
            std::uint_fast32_t teller = 0;            
            while (teller < N) {
                r[teller] = m_history[(first + teller) % N];
                ++teller;
            }
        }
        return r;
    }

    std::uint_fast32_t getSize() const
    {
        return std::min(N, m_index);
    }

private:
    std::array<std::array<double, 3>, N> m_history;
    std::uint_fast32_t m_index = 0;
};

/**
 * @brief Simple struct to describe a photon
 */

struct Particle {
    /**
     * @brief Position vector in three dimensions.
     */
    std::array<double, 3> pos;
    /**
     * @brief Direction vector in three dimension. This vector is threated as a normal vector.
     */
    std::array<double, 3> dir;
    /**
     * @brief Photon energy in keV.
     */
    double energy;
    /**
     * @brief Photon relative weight.
     */
    double weight;

    ParticleTrack<> track;

    inline void translate(const double dist)
    {
        track.registerPosition(pos);
        pos[0] += dir[0] * dist;
        pos[1] += dir[1] * dist;
        pos[2] += dir[2] * dist;
        
    }

    inline static constexpr double border_translate_minimum()
    {
        return GEOMETRIC_ERROR<double>();
    }

    inline void border_translate(const double dist)
    {
        // Make sure we translate particle beyond any border we want to translate to
        // We simply add 100 nm to the distance, works for float and double
        translate(dist + border_translate_minimum());
    }
};
}