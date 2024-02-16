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

class CTAECFilter {
public:
    CTAECFilter()
    {
        m_data.resize(2);
        std::fill(m_data.begin(), m_data.end(), 1.0);
    }

    CTAECFilter(const std::array<double, 3>& start, const std::array<double, 3>& stop, const std::vector<double>& values)
    {
        setData(start, stop, values);
    }

    void normalizeBetween(const std::array<double, 3>& start, const std::array<double, 3>& stop)
    {
        const auto dist_start = vectormath::subtract(start, m_start);
        const auto proj_start = vectormath::dot(dist_start, m_dir);
        const auto dist_stop = vectormath::subtract(stop, m_start);
        const auto proj_stop = vectormath::dot(dist_stop, m_dir);
        normalize(proj_start, proj_stop);
    }

    std::size_t size() const
    {
        return m_data.size();
    }

    bool isEmpty() const
    {
        return m_data.size() == 2;
    }

    const std::vector<double>& weights() const
    {
        return m_data;
    }

    const std::array<double, 3>& start() const
    {
        return m_start;
    }

    std::array<double, 3> stop() const
    {
        return dxmc::vectormath::add(m_start, vectormath::scale(m_dir, m_length));
    }

    double length() const
    {
        return m_length;
    }

    void setData(const std::array<double, 3>& start, const std::array<double, 3>& stop, const std::vector<double>& data)
    {
        m_start = start;
        const auto dir = vectormath::subtract(stop, start);
        m_dir = vectormath::normalized(dir);
        m_length = vectormath::length(dir);

        if (data.size() < 2) {
            m_data.resize(2);
            std::fill(m_data.begin(), m_data.end(), 1.0);
            m_step = m_length;
            return;
        }

        m_step = m_length / (data.size() - 1);
        m_data = data;
        // normalize
        normalize();
    }

    double operator()(const std::array<double, 3>& pos) const
    {
        const auto dist = vectormath::subtract(pos, m_start);
        const auto proj = vectormath::dot(dist, m_dir);
        return this->operator()(proj);
    }

    double operator()(double d) const
    {
        const auto dc = std::clamp(d, 0.0, m_length);
        const auto idx0 = std::clamp(static_cast<std::size_t>(dc / m_step), std::size_t { 0 }, m_data.size() - 2);
        const auto idx1 = idx0 + 1;
        return interp(m_step * idx0, m_step * idx1, m_data[idx0], m_data[idx1], dc);
    }

    double integrate() const
    {
        double s = 0;
        for (std::size_t i = 0; i < m_data.size() - 1; ++i) {
            const auto d0 = m_data[i];
            const auto d1 = m_data[i + 1];
            s += (d0 + d1) * m_step * 0.5;
        }
        return s;
    }

    double integrate(double start_r, double stop_r) const
    {
        const auto start = std::clamp(std::min(start_r, stop_r), 0.0, m_length);
        const auto stop = std::clamp(std::max(start_r, stop_r), 0.0, m_length);
        const auto idx_start = std::clamp(static_cast<std::size_t>(start / m_step), std::size_t { 0 }, m_data.size() - 2);
        const auto idx_stop = std::clamp(static_cast<std::size_t>(stop / m_step) + 1, std::size_t { 0 }, m_data.size() - 2);

        double s = 0;
        for (std::size_t i = idx_start; i < idx_stop; ++i) {
            const auto d0 = m_data[i];
            const auto d1 = m_data[i + 1];
            s += (d0 + d1) * m_step * 0.5;
        }

        const auto begin_part = (m_step * idx_start - start) * (this->operator()(start) + m_data[idx_start + 1]) / 2;
        const auto end_part = (m_step * idx_stop - stop) * (this->operator()(stop) + m_data[idx_stop]) / 2;

        return s - begin_part - end_part;
    }

    double integrate(const std::array<double, 3>& start, const std::array<double, 3>& stop) const
    {
        const auto dist_start = vectormath::subtract(start, m_start);
        const auto proj_start = vectormath::dot(dist_start, m_dir);
        const auto dist_stop = vectormath::subtract(stop, m_start);
        const auto proj_stop = vectormath::dot(dist_stop, m_dir);
        return integrate(proj_start, proj_stop);
    }

protected:
    void normalize()
    {
        const auto area = integrate();
        // we want the total area equal to m_length * 1 for an expected value of 1.0;
        const auto k = m_length / area;
        for (auto& d : m_data)
            d *= k;
    }

    void normalize(double start, double stop)
    {
        const auto area = integrate(start, stop);
        // we want the total area equal to m_length * 1 for an expected value of 1.0;
        const auto k = m_length / area;
        for (auto& d : m_data)
            d *= k;
    }

private:
    std::array<double, 3> m_start = { 0, 0, 0 };
    std::array<double, 3> m_dir = { 0, 0, 1 };
    double m_length = 1;
    double m_step = 1;
    std::vector<double> m_data;
};
}