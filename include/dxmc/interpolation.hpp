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
#include "dxmc/floating.hpp"

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

template <Floating T, Floating U>
inline T logloginterp(T x[2], T y[2], U xi)
{
    const T x0 = x[0];
    const T x1 = x[1];
    const T y0 = y[0];
    const T y1 = y[1];
    // we do not test for negative values, instead we rely on std::log10 to return NaN values
    const double value = std::log10(y0) + (std::log10(y1 / y0) / std::log10(x1 / x0) * std::log10(xi / x0)); // std::log 10 always promotes to double
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

template <Floating T, int N = 30>
class CubicSplineInterpolator {
private:
    std::array<T, (N - 1) * 4> m_coefficients;
    std::array<T, N> m_x;
    T m_step = 0;
    T m_start = 0;
    T m_stop = 0;

protected:
    static std::vector<T> thomasPenSplineElimination(const std::vector<T>& h_p, const std::vector<T>& H_p, const std::vector<T>& d_p)
    {
        // Thomas algorithm for gaussian elimination for a trigonal system of equations
        /*
        |b0 c0  0 0  ..  | x0 |   |d0|
        |a1 b1 c1 0  ... | x1 | = |d1|
        |0  a2 b2 c2  ...| x2 |   |d2|
        |0  0  a3 b3 c3  | .. |   |..|

        */
        std::vector<T> h = h_p;
        std::vector<T> H = H_p;
        std::vector<T> d = d_p;

        std::vector<T> b(d.size(), T { 2 });
        for (std::size_t i = 1; i < d.size(); ++i) {
            const T w = h[i - 1] / H[i - 1];
            H[i] -= w * h[i - 1];
            d[i] -= w * d[i - 1];
        }
        std::vector<T> x(d.size());
        x[d.size() - 1] = d[d.size() - 1] / H[d.size() - 1];
        for (int i = d.size() - 2; i >= 0; --i) {
            x[i] = (d[i] - h[i] * x[i + 1]) / H[i];
        }
        return x;
    }

public:
    template <std::regular_invocable<T> F>
    requires std::is_same<std::invoke_result_t<F, T>, T>::value
    CubicSplineInterpolator(const T start, const T stop, F function)
    {
        std::vector<T> y(N);
        m_start = start;
        m_step = (stop - start) / (N - 2);
        for (std::size_t i = 0; i < N; ++i) {
            m_x[i] = m_start + m_step * i;
            y[i] = function(m_x[i]);
        }
        m_stop = m_x.back();

        std::vector<T> h(m_x.size());
        for (std::size_t i = 0; i < m_x.size() - 1; ++i) {
            h[i] = m_x[i + 1] - m_x[i];
        }

        std::vector<T> delta(m_x.size(), T { 1 });
        for (std::size_t i = 0; i < m_x.size() - 1; ++i) {
            delta[i] = (y[i + 1] - y[i]) / h[i];
        }
        std::vector<T> H(m_x.size());
        for (std::size_t i = 1; i < m_x.size(); ++i) {
            H[i] = 2 * (h[i - 1] + h[i]);
        }
        H[0] = H[1];

        std::vector<T> d(m_x.size());
        for (std::size_t i = 1; i < m_x.size(); ++i) {
            d[i] = 6 * (delta[i] - delta[i - 1]);
        }
        d[m_x.size() - 1] = 0;
        d[0] = 0;
        auto s = thomasPenSplineElimination(h, H, d);

        for (std::size_t i = 0; i < m_x.size() - 1; ++i) {
            const auto offset = i * 4;

            m_coefficients[offset + 0] = (s[i] * m_x[i + 1] * m_x[i + 1] * m_x[i + 1] - s[i + 1] * m_x[i] * m_x[i] * m_x[i] + 6 * (y[i] * m_x[i + 1] - y[i + 1] * m_x[i])) / (6 * h[i]);
            m_coefficients[offset + 0] += h[i] * (s[i + 1] * m_x[i] - s[i] * m_x[i + 1]) / 6;
            m_coefficients[offset + 1] = (s[i + 1] * m_x[i] * m_x[i] - s[i] * m_x[i + 1] * m_x[i + 1] + 2 * (y[i + 1] - y[i])) / (2 * h[i]) + h[i] * (s[i] - s[i + 1]) / 6;
            m_coefficients[offset + 2] = (s[i] * m_x[i + 1] - s[i + 1] * m_x[i]) / (2 * h[i]);
            m_coefficients[offset + 3] = (s[i + 1] - s[i]) / (6 * h[i]);
        }
    }
    T operator()(const T x_val) const
    {
        const T x = std::clamp(x_val, m_start, m_stop);

        const std::size_t index = x > m_start ? static_cast<std::size_t>((x - m_start) / m_step) : 0;
        const std::size_t offset = index < N - 1 ? index * 4 : (N - 2) * 4;
        return m_coefficients[offset] + m_coefficients[offset + 1] * x + m_coefficients[offset + 2] * x * x + m_coefficients[offset + 3] * x * x * x;
    }
};

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
}

template <Floating T, std::regular_invocable<T> F>
requires std::is_same<std::invoke_result_t<F, T>, T>::value constexpr T gaussIntegration(const T start, const T stop, const F function)
{
    // F is a class or function that can be called with a Floating type argument and will return a Floating type value
    auto function_points = gaussIntegrationPoints(start, stop);
    std::array<T, 20> function_values;
    std::transform(std::execution::unseq, function_points.cbegin(), function_points.cend(), function_values.begin(), [&](const auto x) { return function(x); });
    return gaussIntegration(start, stop, function_values);
}

template <Floating T>
class Matrix {
public:
    /*
    Construct R x C matrix; example
    3 x 2 =
    |c11 c21|
    |c12 c22|
    |c13 c23|
    */
    Matrix(std::size_t R, std::size_t C, const std::vector<T>& el)
        : m_R(R)
        , m_C(C)
    {
        if (m_R * m_C == el.size()) {
            m_elements = el;
        } else {
            m_elements.resize(m_R * m_C);
            std::fill(m_elements.begin(), m_elements.end(), T { 0 });
        }
    }
    Matrix(std::size_t R, std::size_t C)
        : m_R(R)
        , m_C(C)
    {
        m_elements.resize(m_R * m_C, T { 0 });
    }

    void setElement(std::size_t i, std::size_t j, T val)
    {
        const auto ind = index(i, j);
        m_elements[ind] = val;
    }
    T getElement(std::size_t i, std::size_t j) const
    {
        const auto ind = index(i, j);
        return m_elements[ind];
    }
    Matrix<T> transpose() const
    {
        Matrix<T> other(m_C, m_R);
        for (std::size_t j = 0; j < m_C; ++j)
            for (std::size_t i = 0; i < m_R; ++i) {
                other.setElement(j, i, getElement(i, j));
            }
        return other;
    }

    Matrix<T> reducedRowEchelonForm() const
    {
        Matrix<T> other(m_R, m_C, m_elements);
        std::size_t lead = 0;

        for (std::size_t row = 0; row < m_R; ++row) {
            if (lead >= m_C)
                return other;
            std::size_t i = row;
            while (std::abs(other.getElement(i, lead)) < std::numeric_limits<T>::epsilon()) { // essentially other[i, lead]==0
                ++i;
                if (i >= m_R) {
                    i = row;
                    ++lead;
                    if (lead >= m_C)
                        return other;
                }
            }
            other.swapRows(i, row);
            other.divideRow(row, other.getElement(row, lead));
            for (std::size_t ii = 0; ii < m_R; ++ii) {
                if (ii != row) {
                    other.addMultipleRow(ii, row, -other.getElement(ii, lead));
                }
            }
        }
        return other;
    }

    std::vector<T> solve(const std::vector<T>& b)
    {
        // for solving Mx=b
        Matrix<T> ext(m_R, m_C + 1);
        for (std::size_t i = 0; i < m_R; ++i) {
            for (std::size_t j = 0; j < m_C; ++j) {
                ext.setElement(i, j, getElement(i, j));
            }
        }
        for (std::size_t i = 0; i < m_R; ++i) {
            ext.setElement(i, m_C, b[i]);
        }
        auto ansM = ext.reducedRowEchelonForm();
        std::vector<T> ans(b.size());
        for (std::size_t i = 0; i < m_R; ++i) {
            ans[i] = ext.getElement(i, m_C);
        }
        return ans;
    }

    std::size_t nrows() const { return m_R; }
    std::size_t ncolumns() const { return m_C; }

protected:
    std::size_t index(std::size_t i, std::size_t j) const
    {
        return i * m_C + j;
    }

    void swapRows(std::size_t from, std::size_t to)
    {
        if (from != to) {
            for (std::size_t j = 0; j < m_C; ++j) {
                const auto r = getElement(from, j);
                setElement(from, j, getElement(to, j));
                setElement(to, j, r);
            }
        }
    }

    void divideRow(const std::size_t r, const T value)
    {
        const T factor = T { 1 } / value;
        for (std::size_t c = 0; c < m_C; ++c) {
            setElement(r, c, getElement(r, c) * factor);
        }
    }

    void addMultipleRow(const std::size_t to, const std::size_t from, const T value)
    {
        for (std::size_t c = 0; c < m_C; ++c) {
            const auto v = getElement(to, c) + getElement(from, c) * value;
            setElement(to, c, v);
        }
    }

private:
    std::vector<T> m_elements;
    std::size_t m_R = 0;
    std::size_t m_C = 0;
};

template <Floating T, int N>
class CubicLSInterpolator {
public:
    CubicLSInterpolator(const std::vector<T>& x, const std::vector<T>& y, const std::vector<T>& t)
    {
    }
};

}