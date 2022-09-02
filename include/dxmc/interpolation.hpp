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

//#define USE_EIGEN

#ifdef USE_EIGEN
#include "Eigen/Dense"
#endif // USE_EIGEN

#include <iostream>

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
    Matrix() {};
    Matrix(std::size_t R, std::size_t C)
        : m_r(R)
        , m_c(C)
    {
        m_data.resize(m_r * m_c, T { 0 });
    }
    Matrix(std::size_t R, std::size_t C, const std::vector<T>& d)
        : m_r(R)
        , m_c(C)
    {
        m_data.resize(m_r * m_c);
        std::copy(d.cbegin(), d.cend(), m_data.begin());
    }

    T& operator()(std::size_t r, std::size_t c)
    {
        return m_data[index(r, c)];
    }
    const T& operator()(std::size_t r, std::size_t c) const
    {

        return m_data[index(r, c)];
    }
    std::vector<T> operator*(const std::vector<T>& x)
    {
        std::vector<T> y(x.size(), 0);
        for (std::size_t r = 0; r < m_r; ++r) {
            for (std::size_t c = 0; c < m_c; ++c) {
                y[r] += m_data[index(r, c)] * x[c];
            }
        }
        return y;
    }
    std::vector<T> solve(const std::vector<T>& bvalues) const
    {
        // gaussian ellim

        auto b = bvalues;
        Matrix<T> A(m_r, m_c, m_data);

        for (std::size_t k = 0; k < m_r; ++k) {
            auto imax = A.find_row_argmax(k);
            if (std::abs(A(imax, k)) <= std::numeric_limits<T>::epsilon()) {
                // we have a singular matrix return nans
                std::fill(b.begin(), b.end(), std::numeric_limits<T>::quiet_NaN());
                return b;
            }
            A.swapRows(k, imax);
            std::swap(b[k], b[imax]);
            for (std::size_t i = k + 1; i < m_r; ++i) {
                const auto c = A(i, k) / A(k, k);
                A(i, k) = T { 0 };
                for (std::size_t j = k + 1; j < m_r; ++j) {
                    A(i, j) -= c * A(k, j);
                }
                b[i] -= c * b[k];
            }
        }

        // backsubstitution before return
        std::vector<T> x(b.size());
        for (std::size_t r = m_r - 1; r < m_r; --r) { // note we rely on unsigned wraparound
            T s = 0;
            for (std::size_t c = r + 1; c < m_r; ++c) {
                s += A(r, c) * b[c];
            }
            b[r] = (b[r] - s) / A(r, r);
        }
        return b;
    }

    void fill(const T value)
    {
        std::fill(m_data.begin(), m_data.end(), value);
    }
    void setZero()
    {
        fill(T { 0 });
    }
    void setZero(std::size_t r, std::size_t c)
    {
        m_r = r;
        m_c = c;
        m_data.resize(m_r * m_c);
        fill(T { 0 });
    }

    std::size_t rows() const { return m_r; }
    std::size_t cols() const { return m_c; }

    Matrix<T> transpose() const
    {
        Matrix<T> t(m_c, m_r);
        for (std::size_t i = 0; i < m_r; ++i) {
            for (std::size_t j = 0; j < m_c; ++j) {
                t(j, i) = m_data[index(i, j)];
            }
        }
        return t;
    }

protected:
    void swapRows(std::size_t from, std::size_t to)
    {
        if (from == to)
            return;
        for (std::size_t i = 0; i < m_c; ++i) {
            std::swap(m_data[index(from, i)], m_data[index(to, i)]);
        }
    }

    std::size_t index(std::size_t r, std::size_t c) const
    {
        return r * m_c + c;
    }

    std::size_t find_row_argmax(std::size_t r) const
    {
        auto imax = r;
        T max_pivot = std::abs(m_data[index(r, r)]);
        for (std::size_t i = r + 1; i < m_r; ++i) {
            const T a = std::abs(m_data[index(i, r)]);
            if (a > max_pivot) {
                max_pivot = a;
                imax = i;
            }
        }
        return imax;
    }

private:
    std::vector<T> m_data;
    std::size_t m_r = 0;
    std::size_t m_c = 0;
};

template <Floating T>
class CubicLSInterpolator {
public:
    CubicLSInterpolator(const std::vector<T>& x, const std::vector<T>& y, const std::vector<T>& t, bool maybe_dicont = true)
    {
        if (maybe_dicont) {
            // all this to handle dicontinous functions designated by a equal x value
            std::size_t x_start = 0;
            for (std::size_t i = 0; i < x.size() - 1; ++i) {
                constexpr auto epsilon = std::numeric_limits<T>::epsilon();
                if (x[i + 1] - x[i] < epsilon) {
                    // ups, we have a discontinous function
                    auto x_beg = x_start;
                    x_start = i + 1;

                    std::vector<T> x_red, y_red, t_red;
                    for (std::size_t j = x_beg; j < x_start; ++j) {
                        x_red.push_back(x[j]);
                        y_red.push_back(y[j]);
                    }
                    t_red.push_back(x[x_beg]);
                    for (std::size_t j = 0; j < t.size(); ++j) {
                        if ((t[j] < x[i]) && (t[j] > x[x_beg])) {
                            t_red.push_back(t[j]);
                        }
                    }
                    t_red.push_back(x[i]);
                    addLSSplinePart(x_red, y_red, t_red);
                }
            }
            // did we se any discontinoutis?
            if (x_start == 0) {
                addLSSplinePart(x, y, t);
            } else {
                std::vector<T> x_red, y_red, t_red;
                for (std::size_t j = x_start; j < x.size(); ++j) {
                    x_red.push_back(x[j]);
                    y_red.push_back(y[j]);
                }
                t_red.push_back(x[x_start]);
                for (std::size_t j = 0; j < t.size(); ++j) {
                    if (t[j] > x[x_start]) {
                        t_red.push_back(t[j]);
                    }
                }

                addLSSplinePart(x_red, y_red, t_red);
            }
        } else {
            addLSSplinePart(x, y, t);
        }

        m_t.shrink_to_fit();
        m_z.shrink_to_fit();
        m_zp.shrink_to_fit();
    }

    T operator()(T x) const
    {
        auto it = std::upper_bound(m_t.cbegin() + 1, m_t.cend() - 1, x);
        const auto h = *it - *(--it);
        const auto u = (x - *it) / h;
        const auto v = u - 1;
        const auto uu = u * u;
        const auto vv = v * v;

        const auto ind = std::distance(m_t.cbegin(), it);

        const auto f1 = (2 * u + 1) * vv * m_z[ind];
        const auto f2 = uu * (1 - 2 * v) * m_z[ind + 1];
        const auto f3 = h * u * vv * m_zp[ind];
        const auto f4 = h * uu * v * m_zp[ind + 1];
        return f1 + f2 + f3 + f4;
    }

protected:
    void addLSSplinePart(const std::vector<T>& x, const std::vector<T>& y, const std::vector<T>& t)
    { // assume x is sorted
        const std::size_t N = t.size() - 1;

        std::vector<T> h(N);
        for (std::size_t i = 0; i < N; ++i)
            h[i] = t[i + 1] - t[i];

        // deal with t intervals without data
        std::vector<std::size_t> pn(N, 0);
        std::vector<std::size_t> p(N, 0);
        std::size_t j = 0;
        for (std::size_t i = 0; i < x.size(); ++i) {
            if (x[i] <= t[j + 1]) {
                pn[j]++;
            } else {
                j++;
                p[j] = i;
                pn[j] = 0;
                while (x[i] > t[j + 1]) {
                    j++;
                    p[j] = i;
                    pn[j] = 0;
                }
                pn[j] = 1;
            }
        }

        // setting up matrices
        std::vector<T> taa(N, T { 0 });
        std::vector<T> tab(N, T { 0 });
        std::vector<T> tag(N, T { 0 });
        std::vector<T> tad(N, T { 0 });
        std::vector<T> tbb(N, T { 0 });
        std::vector<T> tbg(N, T { 0 });
        std::vector<T> tbd(N, T { 0 });
        std::vector<T> tgg(N, T { 0 });
        std::vector<T> tgd(N, T { 0 });
        std::vector<T> tdd(N, T { 0 });

        std::vector<T> D(N + 1, T { 0 });
        std::vector<T> G(N + 1, T { 0 });

        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = p[i]; j < p[i] + pn[i]; ++j) {
                const auto u = (x[j] - t[i]) / h[i];
                const auto v = u - 1;
                const auto alpha = (2 * u + 1) * v * v;
                const auto beta = u * u * (1 - 2 * v);
                const auto gamma = h[i] * u * v * v;
                const auto delta = h[i] * u * u * v;
                taa[i] += alpha * alpha;
                tab[i] += alpha * beta;
                tag[i] += alpha * gamma;
                tad[i] += alpha * delta;
                tbb[i] += beta * beta;
                tbg[i] += beta * gamma;
                tbd[i] += beta * delta;
                tgg[i] += gamma * gamma;
                tgd[i] += gamma * delta;
                tdd[i] += delta * delta;
                D[i] += 2 * y[j] * alpha;
                D[i + 1] += 2 * y[j] * beta;
                G[i] += 2 * y[j] * gamma;
                G[i + 1] += 2 * y[j] * delta;
            }
        }
#ifdef USE_EIGEN
        using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
#else
        using Matrix = Matrix<T>;
#endif // USE_EIGEN

        Matrix A, B, E, C, F, Z;
        A.setZero(N + 1, N + 1);
        B.setZero(N + 1, N + 1);
        E.setZero(N + 1, N + 1);
        C.setZero(N - 1, N + 1);
        F.setZero(N - 1, N + 1);
        Z.setZero(N - 1, N + 1);

        for (std::size_t i = 0; i < N; ++i) {
            A(i, i) += 2 * taa[i];
            A(i + 1, i + 1) = 2 * tbb[i];
            A(i, i + 1) = 2 * tab[i];
            A(i + 1, i) = A(i, i + 1);

            B(i, i) += 2 * tag[i];
            B(i + 1, i + 1) = 2 * tbd[i];
            B(i, i + 1) = 2 * tbg[i];
            B(i + 1, i) = 2 * tad[i];

            E(i, i) += 2 * tgg[i];
            E(i + 1, i + 1) = 2 * tdd[i];
            E(i, i + 1) = 2 * tgd[i];
            E(i + 1, i) = E(i, i + 1);
        }
        for (std::size_t i = 0; i < N - 1; ++i) {
            C(i, i) = 3 * h[i + 1] / h[i];
            C(i, i + 2) = -3 * h[i] / h[i + 1];
            C(i, i + 1) = -C(i, i) - C(i, i + 2);

            F(i, i) = h[i + 1];
            F(i, i + 1) = 2 * (h[i] + h[i + 1]);
            F(i, i + 2) = h[i];
        }

        auto BT = B.transpose();
        auto CT = C.transpose();
        auto FT = F.transpose();

        Matrix sol(A.rows() + BT.rows() + C.rows(), A.cols() + BT.cols() + CT.cols());
        sol.setZero();

        auto nrows = sol.rows();
        auto ncols = sol.cols();

        for (std::size_t r = 0; r < sol.rows(); ++r) {
            for (std::size_t c = 0; c < sol.cols(); ++c) {
                if ((r < N + 1) && (c < N + 1)) {
                    sol(r, c) = A(r, c);
                } else if ((r < N + 1) && (c >= N + 1) && (c < 2 * (N + 1))) {
                    sol(r, c) = BT(r, c - (N + 1));
                } else if ((r < N + 1) && (c >= 2 * (N + 1))) {
                    sol(r, c) = CT(r, c - (2 * (N + 1)));
                } else if ((r >= N + 1) && (r < 2 * (N + 1)) && (c < N + 1)) {
                    sol(r, c) = B(r - (N + 1), c);
                } else if ((r >= N + 1) && (r < 2 * (N + 1)) && (c >= N + 1) && (c < 2 * (N + 1))) {
                    sol(r, c) = E(r - (N + 1), c - (N + 1));
                } else if ((r >= N + 1) && (r < 2 * (N + 1)) && (c >= 2 * (N + 1))) {
                    sol(r, c) = FT(r - (N + 1), c - 2 * (N + 1));
                } else if ((r >= 2 * (N + 1)) && (c < N + 1)) {
                    sol(r, c) = C(r - 2 * (N + 1), c);
                } else if ((r >= 2 * (N + 1)) && (c >= N + 1) && (c < 2 * (N + 1))) {
                    sol(r, c) = F(r - 2 * (N + 1), c - (N + 1));
                } else if ((r >= 2 * (N + 1)) && (c >= 2 * (N + 1))) {
                    sol(r, c) = 0;
                }
            }
        }
#ifdef USE_EIGEN
        Eigen::VectorX<T> bval(sol.rows());
        bval.setZero();
#else
        std::vector<T> bval(sol.rows(), T { 0 });
#endif // USE_EIGEN
        std::copy(D.begin(), D.end(), bval.begin());
        std::copy(G.begin(), G.end(), bval.begin() + N + 1);

#ifdef USE_EIGEN
        // Eigen::VectorX<T> res = sol.colPivHouseholderQr().solve(bval);
        Eigen::VectorX<T> res = sol.lu().solve(bval);
#else
        auto res = sol.solve(bval);
#endif // USE_EIGEN
        std::copy(res.begin(), res.begin() + N + 1, std::back_inserter(m_z));
        std::copy(res.begin() + N + 1, res.begin() + 2 * (N + 1), std::back_inserter(m_zp));
        std::copy(t.begin(), t.end(), std::back_inserter(m_t));
        // m_z = std::vector<T>(res.begin(), res.begin() + N + 1);
        // m_zp = std::vector<T>(res.begin() + N + 1, res.begin() + 2 * (N + 1));
    }

private:
    std::vector<T> m_z;
    std::vector<T> m_zp;
    std::vector<T> m_t;
};
}