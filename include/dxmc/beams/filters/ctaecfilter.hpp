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
#include <pair>
#include <vector>

namespace dxmc {

template <Floating T>
class CTAECFilter {
public:
    CTAECFilter(const std::array<T, 3>& start, const std::array<T, 3>& stop, const std::vector<std::pair<T, T>>& values)
    {
    }

    void setStart(const std::array<T, 3>& pos)
    {
        m_start = pos;
    }

    void setStop(const std::array<T, 3>& pos)
    {
        m_stop = pos;
    }

    void setData(const std::vector<std::pair<T, T>>& values)
    {
        m_data = data;
        // Sorting based on pos
        std::sort(m_data.begin(), m_data.end(), [](const auto& lh, const auto& rh) { return lh.first < rh.first; });

        // normalize to mean of 1
        const auto E = expectationValue(m_data);
        std::for_each(std::execution::par_unseq, m_data.begin(), m_data.end(), [=](auto& d) { d.second /= E; });
    }
    void setData(const std::array<T, 3>& start, const std::array<T, 3>& stop, const std::vector<std::pair<T, T>>& values)
    {
        setStart(start);
        setStop(stop);
        setData(data);
    }

    T operator()(const std::array<T, 3>& pos) const
    {
        if constexpr (ONESIDED)
            return m_inter(std::abs(angle));
        else
            return m_inter(angle);
    }

protected:
    static T expectationValue(std::vector<std::pair<T, T>> data)
    {
        const auto N = data.size();
        if (N < 2)
            return T { 1 };

        T in = 0;
        T ine = 0;
        for (std::size_t i = 1; i < N; ++i) {
            const auto& p0 = data[i - 1];
            const auto& p1 = data[i];
            const auto dx = p1.first - p0.first;
            fortsett her
        }
        return in;
    }

private:
    std::array<T, 3> m_start = { 0, 0, 0 };
    std::array<T, 3> m_stop = { 0, 0, 0 };
    std::vector<std::pair<T, T>> m_data;
};
}