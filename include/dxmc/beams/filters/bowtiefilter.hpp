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

Copyright 2023 Erlend Andersen
*/

#pragma once

#include "dxmc/floating.hpp"
#include "dxmc/interpolation.hpp"

#include <array>
#include <vector>

namespace dxmc {

template <Floating T, bool ONESIDED = true>
class BowtieFilter {
public:
    BowtieFilter(const std::vector<T>& angles_r, const std::vector<T>& intensity_r)
    {
        const auto N = std::min(angles_r.size(), intensity_r.size());
        std::vector<std::pair<T, T>> data(N);

        for (std::size_t i = 0; i < N; ++i) {
            if constexpr (ONESIDED)
                data[i].first = std::abs(angles_r[i]);
            else
                data[i].first = angles_r[i];
            data[i].second = std::abs(intensity_r[i]);
        }

        std::sort(data.begin(), data.end(), [](const auto& lh, const auto& rh) { return lh.first < rh.first; });

        m_inter.setup(data);

        const auto start = data.front().first;
        const auto stop = data.back().first;
        const auto area = m_inter.integral(start, stop);
        const auto scale = (stop - start) / area;

        m_inter.scale((stop - start) / area);
    }

    T operator()(T angle) const
    {
        if constexpr (ONESIDED)
            return m_inter(std::abs(angle));
        else
            return m_inter(angle);
    }

private:
    AkimaSplineStatic<T, 6> m_inter;
};
}