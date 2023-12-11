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

namespace dxmc {

template <Floating T, bool ONESIDED = true>
class BowtieFilter {
public:
    BowtieFilter(const std::vector<T>& angles_r, const std::vector<T>& intensity_r)
    {
        setData(angles_r, intensity_r);
    }

    BowtieFilter(const std::vector<std::pair<T, T>>& data)
    {
        setData(data);
    }

    BowtieFilter()
    {
        // generic filter from a Siemens Definition Flash
        static std::vector<std::pair<T, T>> data = {
            { 0.166511074f, 3.53208f },
            { 0.000000000f, 13.9167f },
            { 0.041992107f, 12.5868f },
            { 0.083836642f, 9.41943f },
            { 0.246954945f, 1.96665f },
            { 0.324269441f, 1.27605f },
            { 0.390607044f, 0.947716f }
        };
        setData(data);
    }

    T operator()(T angle) const
    {
        if constexpr (ONESIDED)
            return m_inter(std::abs(angle));
        else
            return m_inter(angle);
    }

    void setData(std::vector<std::pair<T, T>> data)
    {
        for (auto& d : data) {
            if constexpr (ONESIDED)
                d.first = std::abs(d.first);
            d.second = std::abs(d.second);
        }
        std::sort(data.begin(), data.end(), [](const auto& lh, const auto& rh) { return lh.first < rh.first; });

        m_inter.setup(data);

        const auto start = data.front().first;
        const auto stop = data.back().first;
        const auto area = m_inter.integral(start, stop);
        const auto scale = (stop - start) / area;

        m_inter.scale((stop - start) / area);
    }

    void setData(const std::vector<T>& angles_r, const std::vector<T>& intensity_r)
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

        setData(data);
    }

private:
    AkimaSplineStatic<T, 6> m_inter;
};
}