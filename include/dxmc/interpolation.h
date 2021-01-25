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

#pragma once // include guard
#include "dxmc/floating.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <execution>
#include <iterator>
#include <numeric>
#include <type_traits>
#include <vector>

namespace dxmc {

template <Floating T>
inline T interp(T x0, T x1, T y0, T y1, T x)
{
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}
template <Floating T, Floating U>
inline U interp(T x[2], T y[2], U xi)
{
    return y[0] + (y[1] - y[0]) * (xi - x[0]) / (x[1] - x[0]);
}

template <Floating T>
inline T logloginterp(T x0, T x1, T y0, T y1, T x)
{
    // we do not test for negative values, instead we rely on std::log10 to return NaN values
    const double value = std::log10(y0) + (std::log10(y1 / y0) / std::log10(x1 / x0) * std::log10(x / x0)); // std::log 10 always promotes to double
    return std::pow(10., value); // std::pow always promotes to doubles
}

template <typename It, Floating T>
requires std::is_same_v<typename std::iterator_traits<It>::value_type, T>
    T interpolate(It xbegin, It xend, It ybegin, It yend, T xvalue)
{
    auto upper = std::upper_bound(xbegin, xend, xvalue);
    if (upper == xbegin)
        return *ybegin;
    if (upper == xend)
        return *(yend - 1);
    auto lower = upper;
    std::advance(lower, -1);

    auto lowery = ybegin;
    std::advance(lowery, std::distance(xbegin, lower));
    auto uppery = ybegin;
    std::advance(uppery, std::distance(xbegin, upper));

    return interp(*lower, *upper, *lowery, *uppery, xvalue);
}

template <Floating T>
std::vector<T> trapz(const std::vector<T>& f,
    const std::vector<T>& x)
{
    std::vector<T> integ(f.size(), 0);
    integ[0] = T { 0.0 };
    for (std::size_t i = 1; i < f.size(); ++i) {
        integ[i] = integ[i - 1] + (f[i - 1] + f[i]) * T { 0.5 } * (x[i] - x[i - 1]);
    }
    return integ;
}

template <Floating T>
constexpr T gaussIntegration(const T start, const T stop, std::array<T, 20> gaussPoints)
{
    constexpr std::array<T, 20> weights = {
        1.5275338713072585E-01,
        1.4917298647260375E-01,
        1.4209610931838205E-01,
        1.3168863844917663E-01,
        1.1819453196151842E-01,
        1.0193011981724044E-01,
        8.3276741576704749E-02,
        6.2672048334109064E-02,
        4.0601429800386941E-02,
        1.7614007139152118E-02,
        1.5275338713072585E-01,
        1.4917298647260375E-01,
        1.4209610931838205E-01,
        1.3168863844917663E-01,
        1.1819453196151842E-01,
        1.0193011981724044E-01,
        8.3276741576704749E-02,
        6.2672048334109064E-02,
        4.0601429800386941E-02,
        1.7614007139152118E-02
    };
    const T interval_half = (stop - start) * T { 0.5 };
    const T value = std::transform_reduce(std::execution::unseq, weights.cbegin(), weights.cend(), gaussPoints.cbegin(), T { 0 }, std::plus<>(), std::multiplies<>());
    return value * interval_half;
}
template <Floating T>
constexpr std::array<T, 20> gaussIntegrationPoints(const T start, const T stop)
{
    constexpr std::array<T, 20> x_val = {
        7.6526521133497334E-02,
        2.2778585114164508E-01,
        3.7370608871541956E-01,
        5.1086700195082710E-01,
        6.3605368072651503E-01,
        7.4633190646015079E-01,
        8.3911697182221882E-01,
        9.1223442825132591E-01,
        9.6397192727791379E-01,
        9.9312859918509492E-01,
        -7.6526521133497334E-02,
        -2.2778585114164508E-01,
        -3.7370608871541956E-01,
        -5.1086700195082710E-01,
        -6.3605368072651503E-01,
        -7.4633190646015079E-01,
        -8.3911697182221882E-01,
        -9.1223442825132591E-01,
        -9.6397192727791379E-01,
        -9.9312859918509492E-01
    };
    std::array<T, 20> function_points;
    const T xi = (stop - start) * T { 0.5 };
    const T xii = (stop + start) * T { 0.5 };
    std::transform(std::execution::unseq, x_val.cbegin(), x_val.cend(), function_points.begin(), [=](const auto x) { return x * xi + xii; });
    return function_points;
};

template <Floating T, std::regular_invocable<T> F>
requires std::is_same<std::invoke_result_t<F, T>, T>::value constexpr T gaussIntegration(const T start, const T stop, const F function)
{
    // F is a class or function that can be called with a Floating type argument and will return a Floating type value
    auto function_points = gaussIntegrationPoints(start, stop);
    std::array<T, 20> function_values;
    std::transform(std::execution::unseq, function_points.cbegin(), function_points.cend(), function_values.begin(), [&](const auto x) { return function(x); });
    return gaussIntegration(start, stop, function_values);
}

//class for numerical inverse transform of analytical probability density functions (pdfs)
template <Floating T, int N = 20>
class RITA {
private:
    std::array<T, N> m_x;
    std::array<T, N> m_e;
    std::array<T, N> m_a;
    std::array<T, N> m_b;

protected:
    template <typename F>
    static T simpson_integral(const T start, const T stop, F pdf)
    {
        const T h = (stop - start) / 50;
        T result = pdf(start) + pdf(stop);
        for (std::size_t i = 1; i < 50; ++i) {
            const T prod = 1 % 2 == 0 ? 2 : 4; // prod is 2 when i is even else 4
            result += prod * pdf(start + h * i);
        }
        return h * result / 3;
    }
    static void insert_elements(std::array<T, N>& array, std::size_t index, T value)
    {
        for (std::size_t i = N - 1; i > index; --i) {
            array[i] = array[i - 1];
        }
        array[index + 1] = value;
    }
    T integral_p_bar(std::size_t index) const
    {
        const T ai = m_a[index];
        const T bi = m_b[index];
        const T xi = m_x[index];
        const T xii = m_x[index + 1];

        auto p = [=](const T x) -> T {
            T n;
            if (x == xi) {
                n = 0;
            } else if (x == xii) {
                n = 1;
            } else {
                const T t = (x - xi) / (xii - xi);
                const T f = (1 + ai + bi - ai * t) / (2 * bi * t);
                const T nom = 1 + ai + bi - ai * t;
                const T l = 1 - std::sqrt(1 - (4 * bi * t * t) / (nom * nom));
                n = f * l;
            }
            const T upper = (1 + ai * n * bi * n * n);
            const T lower = (1 + ai + bi) * (1 - bi * n * n);
            const T res = upper * (m_e[index + 1] - m_e[index]) / (lower * (xii - xi));
            return res;
        };
        return simpson_integral(xi, xii, p);
    }

public:
    template <std::regular_invocable<T> F>
    requires std::is_same<std::invoke_result_t<F, T>, T>::value
    RITA(const T min, const T max, F pdf)
    {
        m_x[0] = min;
        while (pdf(m_x[0]) <= 0) {
            m_x[0] += (max - min) / (N*10);
        }
        m_x[1] = (max - min) / 2;
        m_x[2] = max;
        m_e.fill(1);
        m_e[0] = 0;

        std::array<T, N> error;
        error.fill(0);

        std::size_t current_point = 2;
        while (current_point < N ) {

            for (std::size_t i = 1; i <= current_point; ++i) {
                m_e[i] = simpson_integral(m_x[i - 1], m_x[i], pdf) + m_e[i - 1];
            }

            for (std::size_t i = 0; i < current_point; ++i) {
                const T frac = (m_e[i + 1] - m_e[i]) / (m_x[i + 1] - m_x[i]);

                m_b[i] = 1 - frac * frac / (pdf(m_x[i + 1]) * pdf(m_x[i]));
                m_a[i] = frac / pdf(m_x[i]) - m_b[i] - 1;

                error[i] = std::abs(simpson_integral(m_x[i], m_x[i + 1], pdf) - integral_p_bar(i));
            }
            const auto it_max = std::max_element(error.cbegin(), error.cbegin() + current_point);
            const std::size_t error_index = std::distance(error.cbegin(), it_max);
            T new_value=(m_x[error_index] + m_x[error_index + 1]) / 2;

            insert_elements(m_x, error_index, new_value);
            current_point++;
        }
        
        std::transform(m_e.cbegin(), m_e.cend(), m_e.begin(), [=](const auto e) { return e/ m_e[N - 1]; });
    }

    T operator()(const T random) const
    {
        auto upper_bound = std::upper_bound(m_e.cbegin(), m_e.cend(), random);
        const std::size_t index = std::distance(m_e.cbegin(), upper_bound);
        if (index == 0)
            return m_x[0];

        const auto v = random - m_e[index - 1];
        const auto d = m_e[index] - m_e[index - 1];
        return m_x[index - 1] + (1 + m_a[index - 1] + m_b[index - 1]) * d * v / (d * d + m_a[index - 1] * d * v + m_b[index - 1] * v * v) * (m_x[index + 1] - m_x[index]);
    }
};

}