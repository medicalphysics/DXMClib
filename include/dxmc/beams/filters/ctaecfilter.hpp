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
#include "dxmc/vectormath.hpp"

#include <array>
#include <vector>

namespace dxmc {

template <Floating T>
class CTAECFilter {
public:
    CTAECFilter()
    {
        m_data.resize(2);
        std::fill(m_data.begin(), m_data.end(), T { 1 });
    }

    CTAECFilter(const std::array<T, 3>& start, const std::array<T, 3>& stop, const std::vector<T>& values)
    {
        setData(start, stop, values);
    }

    void setData(const std::array<T, 3>& start, const std::array<T, 3>& stop, const std::vector<T>& data)
    {
        m_start = start;
        const auto dir = vectormath::subtract(stop, start);
        m_dir = vectormath::normalized(dir);
        m_length = vectormath::lenght(dir);

        if (data.size() < 2) {
            m_data.resize(2);
            std::fill(m_data.begin(), m_data.end(), T { 1 });
            m_step = m_length;
            return;
        }

        m_step = m_length / (data.size() - 1);
        m_data = data;
        // normalize
        const auto mean = std::reduce(m_data.cbegin(), m_data.cend()) / m_data.size();
        for (auto& d : m_data)
            d /= mean;
    }

    T operator()(const std::array<T, 3>& pos) const
    {
        const auto dist = vectormath::subtract(m_start, pos);
        const auto proj = vectormath::dot(dist, m_dir);
        if (proj > m_length)
            return m_data.back();
        else if (proj < 0)
            return m_data.front();

        const auto idx0 = std::clamp(static_cast<std::size_t>(proj / m_step), std::size_t { 0 }, m_data.size() - 2);
        const auto idx1 = idx0 + 1;
        return interp(m_step * idx0, m_step * idx1, m_data[idx0], m_data[idx1], proj);
    }

protected:
private:
    std::array<T, 3> m_start = { 0, 0, 0 };
    std::array<T, 3> m_dir = { 0, 0, 1 };
    T m_length = 1;
    T m_step = 1;
    std::vector<T> m_data;
};
}