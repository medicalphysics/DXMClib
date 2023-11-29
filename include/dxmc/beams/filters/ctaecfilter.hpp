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

    std::size_t size() const
    {
        return m_data.size();
    }

    const std::vector<T>& weights() const
    {
        return m_data;
    }

    const std::array<T, 3>& start() const
    {
        return m_start;
    }

    std::array<T, 3> stop() const
    {
        dxmc::vectormath::add(m_start, vectormath::scale(m_dir, m_lenght));
    }

    T length() const
    {
        return m_lenght;
    }

    void setData(const std::array<T, 3>& start, const std::array<T, 3>& stop, const std::vector<T>& data)
    {
        m_start = start;
        const auto dir = vectormath::subtract(stop, start);
        m_dir = vectormath::normalized(dir);
        m_lenght = vectormath::lenght(dir);

        if (data.size() < 2) {
            m_data.resize(2);
            std::fill(m_data.begin(), m_data.end(), T { 1 });
            m_step = m_lenght;
            return;
        }

        m_step = m_lenght / (data.size() - 1);
        m_data = data;
        // normalize
        normalize();
    }

    T operator()(const std::array<T, 3>& pos) const
    {
        const auto dist = vectormath::subtract(m_start, pos);
        const auto proj = vectormath::dot(dist, m_dir);
        return this->operator()(proj);
    }

    T operator()(T d) const
    {
        const auto dc = std::clamp(d, T { 0 }, m_lenght);
        const auto idx0 = std::clamp(static_cast<std::size_t>(d / m_step), std::size_t { 0 }, m_data.size() - 2);
        const auto idx1 = idx0 + 1;
        return interp(m_step * idx0, m_step * idx1, m_data[idx0], m_data[idx1], dc);
    }

    T integrate() const
    {
        T s = 0;
        for (std::size_t i = 0; i < m_data.size() - 1; ++i) {
            const auto d0 = m_data[i];
            const auto d1 = m_data[i + 1];
            s += (d0 + d1) * m_step * T { 0.5 };
        }
        return s;
    }

    T integrate(T start_r, T stop_r) const
    {
        const auto start = std::clamp(std::min(start_r, stop_r), T { 0 }, m_lenght);
        const auto stop = std::clamp(std::max(start_r, stop_r), T { 0 }, m_lenght);
        const auto idx_start = std::clamp(static_cast<std::size_t>(start / m_step), std::size_t { 0 }, m_data.size() - 2);
        const auto idx_stop = std::clamp(static_cast<std::size_t>(stop / m_step) + 1, std::size_t { 0 }, m_data.size() - 2);

        T s = 0;
        for (std::size_t i = idx_start; i < idx_stop; ++i) {
            const auto d0 = m_data[i];
            const auto d1 = m_data[i + 1];
            s += (d0 + d1) * m_step * T { 0.5 };
        }

        const auto begin_part = (m_step * idx_start - start) * (this->operator()(start) + m_data[idx_start + 1]) / 2;
        const auto end_part = (m_step * idx_stop - stop) * (this->operator()(stop) + m_data[idx_stop]) / 2;

        return s - begin_part - end_part;
    }

protected:
    void normalize()
    {
        const auto area = integrate();
        // we want the total area equal to m_lenght * 1 for an expected value of 1.0;
        const auto k = m_lenght / area;
        for (auto& d : m_data)
            d *= k;
    }

    void normalize(T start, T stop)
    {
        const auto area = integrate(start, stop);
        // we want the total area equal to m_lenght * 1 for an expected value of 1.0;
        const auto k = m_lenght / area;
        for (auto& d : m_data)
            d *= k;
    }

private:
    std::array<T, 3> m_start = { 0, 0, 0 };
    std::array<T, 3> m_dir = { 0, 0, 1 };
    T m_lenght = 1;
    T m_step = 1;
    std::vector<T> m_data;
};
}