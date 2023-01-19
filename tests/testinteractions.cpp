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

#include "dxmc/interactions.hpp"
#include "dxmc/material/material.hpp"

#include "xraylib.h"

#include <iostream>

template <typename T>
class Histogram {
public:
    Histogram(T start, T stop, std::size_t Nbins)
        : m_start(start)
        , m_stop(stop)
    {
        intensity.resize(Nbins);
        std::fill(intensity.begin(), intensity.end(), 0);
        edges.resize(Nbins + 1);
        m_step = (stop - start) / Nbins;
        for (std::size_t i = 0; i <= Nbins; ++i) {
            edges[i] = m_start + m_step * i;
        }
    }
    inline std::size_t idx(T val) const
    {
        const T frac = (val - m_start) / (m_stop - m_start);
        const auto idx_f = std::nextafter(frac * intensity.size(), T { 0 });
        const auto idx = static_cast<std::size_t>(idx_f);
        return idx;
    }

    void operator()(T value)
    {
        const auto i = idx(value);
        intensity[i] += 1;
    }
    std::pair<std::vector<T>, std::vector<T>> plotData(bool normalize = true)
    {
        std::vector<T> x(intensity.size(), 0);
        std::vector<T> y(intensity.size(), 0);

        const auto sum = normalize ? std::reduce(intensity.cbegin(), intensity.cend(), T { 0 }) : T { 1 };
        for (std::size_t i = 0; i < intensity.size(); ++i) {
            y[i] = static_cast<T>(intensity[i]) / sum;
            x[i] = (edges[i] + edges[i + 1]) / 2;
        }
        return std::make_pair(x, y);
    }

private:
    std::vector<std::size_t> intensity;
    std::vector<T> edges;
    T m_start;
    T m_stop;
    T m_step = 0;
};

template <typename T, int type = 0>
bool testCoherent(std::size_t Z = 13, T energy = 50)
{
    bool success = true;
    constexpr std::size_t N = 1E7;
    Histogram hist(T { -1 }, T { 1 }, 150);

    auto material_opt = dxmc::Material2<T>::byZ(Z);
    auto material = material_opt.value();

    constexpr std::array<T, 3> dir = { 0, 0, 1 };
    dxmc::Particle<T> particle;
    particle.pos = { 0, 0, 0 };
    particle.dir = dir;
    particle.energy = energy;
    particle.weight = 1;

    dxmc::RandomState state;
    for (std::size_t i = 0; i < N; ++i) {
        particle.dir = dir;
        dxmc::interactions::rayleightScatter<T, type>(particle, material, state);
        const T costheta = dxmc::vectormath::dot(particle.dir, dir);
        hist(costheta);
    }

    const auto [x, y] = hist.plotData();
    std::vector<T> y_ana(x.size(), 0);
    std::vector<T> y_xraylib(x.size(), 0);
    // analytical
    if constexpr (type == 0) {
        std::transform(x.cbegin(), x.cend(), y_ana.begin(), [](const auto xv) {
            return (1 + xv * xv);
        });
        std::transform(x.cbegin(), x.cend(), y_xraylib.begin(), [](const auto xv) {
            const auto angle = std::acos(xv);
            return static_cast<T>(DCS_Thoms(angle, nullptr));
        });
    } else {
        std::transform(x.cbegin(), x.cend(), y_ana.begin(), [energy, &material](const auto xv) {
            const auto angle = std::acos(xv);
            const auto q = material.momentumTransfer(energy, angle);
            const auto formFactor = material.formFactor(q);
            return (1 + xv * xv) * formFactor * formFactor / 2;
        });
        std::transform(x.cbegin(), x.cend(), y_xraylib.begin(), [Z, energy](const auto xv) {
            const auto angle = std::acos(xv);
            return static_cast<T>(DCS_Rayl(Z, energy, angle, nullptr));
        });
    }
    const auto y_ana_sum = std::reduce(y_ana.cbegin(), y_ana.cend(), T { 0 });
    std::for_each(y_ana.begin(), y_ana.end(), [y_ana_sum](auto& yv) { yv /= y_ana_sum; });
    const auto y_xraylib_sum = std::reduce(y_xraylib.cbegin(), y_xraylib.cend(), T { 0 });
    std::for_each(y_xraylib.begin(), y_xraylib.end(), [y_xraylib_sum](auto& yv) { yv /= y_xraylib_sum; });

    std::cout << "Rayleight\ncostheta, samples, analytical, xraylib\n";
    for (std::size_t i = 0; i < y.size(); ++i) {
        std::cout << x[i] << ", " << y[i] << ", " << y_ana[i] << ", " << y_xraylib[i] << std::endl;
    }

    return success;
}

int main()
{
    std::cout << "Testing interactions" << std::endl;
    bool success = true;
    success = success && testCoherent<double, 1>(80, 5.0);
    success = success && testCoherent<float, 1>();
    success = success && testCoherent<double, 0>();
    success = success && testCoherent<float, 0>();

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
