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

template <Floating T>
class BowtieFilter {
public:
    BowtieFilter(const std::vector<T>& angles_r, const std::vector<T>& intensity_r)
    {
        const auto N = std::min(angles_r.size(), intensity_r.size());
        std::vector<std::pair<T, T>> data(N);

        for (std::size_t i = 0; i < N; ++i) {
            data[i].first = angles_r[i];
            data[i].second = std::abs(intensity_r[i]);
        }

        std::sort(data.begin(), data.end(), [](const auto& lh, const auto& rh) { return lh.first < rh.first; });

        const auto mean_intensity = trapz(data);
        for (auto& p : data) {
            p.second /= mean_intensity;
        }
        m_inter.setup(data);
    }

    T operator()(T angle) const
    {
        return m_inter(angle);
    }

private:
    CubicSplineInterpolator<T> m_inter;
};
}