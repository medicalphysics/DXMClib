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
class CubicSplineInterpolator {
private:
    std::vector<T> m_coefficients;
    std::vector<T> m_x;

protected:
    static std::vector<T> gaussSplineElimination(const std::vector<T>& h, const std::vector<T>& D)
    {
        const std::size_t m = h.size() - 1; // rows
        const std::size_t n = m + 1; // columns

        auto idx = [=](int i, int j) -> int { return j + i * m; }; //(row, column)
        //construct matrix
        std::vector<T> A(n * m, 0);
        for (int j = 0; j < m; ++j) {
            for (int i = 0; i < n; ++i) {
                const auto index = idx(i, j);
                if (i == j)
                    A[index] = 2 * (h[i] + h[i + 1]);
                else if (j - i == 1)
                    A[index] = h[i];
                else if (j - i == -1)
                    A[index] = h[i];
                else if (j == n)
                    A[index] = D[j + 1];
            }
        }

        for (int i = 0; i < m - 1; i++) {
            //Partial Pivoting
            for (int k = i + 1; k < m; k++) {
                //If diagonal element(absolute vallue) is smaller than any of the terms below it
                if (std::abs(A[idx(i, i)]) < std::abs(A[idx(k, i)])) {
                    //Swap the rows
                    for (int j = 0; j < n; j++) {
                        const auto temp = A[idx(i, j)];
                        A[idx(i, j)] = A[idx(k, j)];
                        A[idx(k, j)] = temp;
                    }
                }
            }
            //Begin Gauss Elimination
            for (int k = i + 1; k < m; k++) {
                const auto term = A[idx(k, i)] / A[idx(i, i)];
                for (int j = 0; j < n; j++) {
                    A[idx(k, j)] = A[idx(k, j)] - term * A[idx(i, j)];
                }
            }
        }
        std::vector<T> sigma(h.size()+1, 0);

        //Begin Back-substitution
        for (int i = m - 1; i >= 0; i--) {
            sigma[i + 1] = A[idx(i, n - 1)];
            for (int j = i + 1; j < n - 1; j++) {
                sigma[i + 1] = sigma[i + 1] - A[idx(i, j)] * sigma[j + 1];
            }
            sigma[i + 1] = sigma[i + 1] / A[idx(i, i)];
        }
        return sigma;
    }

public:
    CubicSplineInterpolator(const std::vector<T>& x, const std::vector<T>& y)
    {
        if (x.size() != y.size())
            return;

        if (x.size() < 3)
            return;

        m_coefficients.resize(4 * (x.size() - 1));
        m_x = x; // copy array

        std::vector<T> h(m_x.size() - 1);
        std::vector<T> rho(m_x.size() - 1);
        std::vector<T> delta(m_x.size() - 1);
        for (std::size_t i = 0; i < m_x.size() - 1; ++i) {
            h[i] = m_x[i + 1] - m_x[i];
            delta[i] = (y[i + 1] - y[i]) / h[i];
        }

        const T sigma0 = 1;
        const T sigmaN = 1;
        std::vector<T> D(m_x.size() - 1);
        for (std::size_t i = 1; i < D.size(); ++i) {
            D[i] = 6*(delta[i] - delta[i - 1]);
        }
        D[0] = 0;
        D[1] += -h[1] * sigma0;
        D[m_x.size() - 2] += -h[m_x.size() - 2] * sigmaN;

        auto sigma1N = gaussSplineElimination(h, D);
        sigma1N[0] = sigma0;
        sigma1N[m_x.size() - 1] = sigmaN;

        for (std::size_t i = 0; i < m_x.size() - 1; ++i) {
            m_coefficients[i * 4 + 0] = (sigma1N[i] * x[i + 1] * x[i + 1] * x[i + 1] - sigma1N[i + 1] * x[i] * x[i] * x[i] + 6 * (y[i] *x[i + 1] - y[i + 1] * x[i])) / (h[i] * 6);
            m_coefficients[i * 4 + 0] += h[i] * (sigma1N[i + 1] * x[i] - sigma1N[i] * x[i + 1]) / 6;

            m_coefficients[i * 4 + 1] = (sigma1N[i + 1] * x[i] * x[i] - sigma1N[i] * x[i + 1] * x[i + 1] + 2 * (y[i + 1] - y[i])) / (h[i] * 2) + h[i] * (sigma1N[i] - sigma1N[i + 1]) / 6;
            m_coefficients[i * 4 + 2] = (sigma1N[i] * x[i + 1] - sigma1N[i + 1] * x[i]) / (2 * h[i]);
            m_coefficients[i * 4 + 3] = (sigma1N[i + 1] - sigma1N[i]) / (6 * h[i]);
        }
    }
    T operator()(const T x) const
    {
        auto upper_iter = std::upper_bound(m_x.cbegin(), m_x.cend(), x);
        if (upper_iter == m_x.cend() || upper_iter == m_x.cbegin())
            return 0;
        const std::size_t i = std::distance(m_x.cbegin(), upper_iter) - 1;
        return m_coefficients[i * 4] + m_coefficients[i * 4 + 1] * x + m_coefficients[i * 4 + 2] * x * x + m_coefficients[i * 4 + 3] * x * x * x;
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

}