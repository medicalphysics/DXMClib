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

#include <assert.h>
#include <concepts>
#include <execution>
#include <limits>
#include <numeric>
#include <random>
#include <utility>
#include <vector>

namespace dxmc {

/**
 * @brief Class for simple generation of random numbers
 * This class aims to provide a small and fast PRNG, but should perhaps be replaced by a STL random generator.
 */
class RandomState {
public:
    /**
     * @brief Initiate a RandomState with a seed from the local machine random device implementation.
     */
    RandomState()
    {
        std::random_device d;
        std::uniform_int_distribution<std::uint64_t> dist(0);
        m_state[0] = static_cast<std::uint64_t>(dist(d));
        m_state[1] = static_cast<std::uint64_t>(dist(d));
    }
    /**
     * @brief Initiate RandomState with a custom seed
     * @param state Pointer to a array with two unsigned 64 bit numbers. Both numbers must not be zero.
     */
    RandomState(std::uint64_t state[2])
    {
        m_state[0] = state[0];
        m_state[1] = state[1];
    }
    RandomState(const RandomState&) = delete; // non construction-copyable
    RandomState& operator=(const RandomState&) = delete; // non copyable

    /**
     * @brief Generate a random floating point number in range [0, 1.0)
     * @tparam T Must satisfy std::is_floating_point<T>::value == True
     * @return Random floating point in range [0, 1)
     */
    template <typename T>
    inline T randomUniform() noexcept
    {
        static_assert(std::is_floating_point<T>::value, "Uniform random number requires floating point precision");
        const auto ui = pcg32();
        constexpr T uiMaxInv = T { 2.32830643653869628906e-010 };
        return ui * uiMaxInv;
    }

    /**
     * @brief Generate random uniform number in interval from 0 to max, exclusive
     * @tparam T Type of number, either integral or floating point
     * @param max Max of range
     * @return Random number in range [0, max).
     */
    template <typename T>
    inline T randomUniform(const T max) noexcept
    {
        if constexpr (std::is_floating_point<T>::value) {
            return randomUniform<T>() * max;
        } else if constexpr (std::is_integral<T>::value) {
            const std::uint32_t threshold = static_cast<std::uint32_t>(-max % max);
            for (;;) {
                const auto r = pcg32();
                if (r >= threshold)
                    return static_cast<T>(r % static_cast<std::uint32_t>(max));
            }
        } else {
            static_assert(std::is_integral<T>::value || std::is_floating_point<T>::value, "Must be integral or floating point value.");
        }
    }

    /**
     * @brief Generate random uniform number in interval from 0 to max, exclusive
     * @tparam T Type of number, either integral or floating point.
     * @param min Min of range.
     * @param max Max of range.
     * @return Random number in range [min, max).
     */
    template <typename T>
    inline T randomUniform(const T min, const T max) noexcept
    {
        if constexpr (std::is_floating_point<T>::value) {
            const T r = randomUniform<T>();
            const T range = max - min;
            return min + r * range;
        } else if constexpr (std::is_integral<T>::value) {
            const T range = max - min;
            const T r = randomUniform<T>(range);
            return min + r;
        } else
            static_assert(std::is_integral<T>::value || std::is_floating_point<T>::value, "Must be integral or floating point value.");
    }

    /**
     * @brief Generate random unsigned integer on interval 0 to max inclusive
     * @param max
     * @return
     */
    template <std::unsigned_integral T>
    inline T randomInteger(const T max) noexcept
    {
        static_assert(sizeof(max) <= 4, "This prng only supports up to 32 bit random integers, for a capped to 32 bit random integer use randomInteger32BitCapped instead");
        const T threshold = (static_cast<T>(-max)) % max;
        for (;;) {
            const auto r = pcg32();
            if (r >= threshold)
                return r % max;
        }
    }

    template <std::unsigned_integral T>
    inline T randomInteger32BitCapped(const T max) noexcept
    {
        static_assert(sizeof(max) > 4, "This function is intended for 64 bit values or greater, use randomInteger method instead");
        const T threshold = (static_cast<T>(-max)) % max;
        for (;;) {
            const auto r = pcg32();
            if (r >= threshold)
                return r % max;
        }
    }

    /**
     * @brief Implementation of the pcg32 random number generator
     * The function below is borrowed from:
     * Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
     * Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)
     * @return A random 32 bit unsigned integer.
     */
    inline std::uint32_t pcg32() noexcept
    {
        const std::uint64_t oldstate = m_state[0];
        // Advance internal state
        m_state[0] = oldstate * 6364136223846793005ULL + (m_state[1] | 1);
        // Calculate output function (XSH RR), uses old state for max ILP
        const std::uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
        const std::uint32_t rot = oldstate >> 59u;
        return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
    }
    std::uint64_t m_state[2];
};

/**
 * @brief Class for sampling of random numbers from a specified descreet distribution.
 * Class for sampling of random numbers from a specified descreet distribution. This implementation use the filling of histogram method.
 */
template <Floating T = double>
class RandomDistribution {
public:
    /**
     * @brief Initialize RandomDistribution class
     * @param weights A vector of probabilities for each bin. Weights should be atleast of size two. The weights vector will be normalized such that sum of weights is unity.
     */
    RandomDistribution(const std::vector<T>& weights)
    {
        m_data = generateTable(weights);
    }
    RandomDistribution()
    {
        m_data.push_back(std::make_pair(T { 1 }, 0));
    }

    /**
     * @brief Sample an index from 0 to size() - 1 according to weights. This function is thread safe.
     * @param state Random state for sampling of an index.
     * @return index with probability according to weights[index].
     */
    std::size_t sampleIndex(RandomState& state) const
    {
        const auto r = state.randomUniform<T>();
        const auto k = state.randomUniform<std::size_t>(m_data.size());
        // return r < m_probs[k] ? k : m_alias[k];
        return r < m_data[k].first ? k : m_data[k].second;
    }

    static std::vector<std::pair<T, std::uint64_t>> generateTable(const std::vector<T>& weights)
    {
        // Squaring the histogram method
        const auto m_size = weights.size();

        std::vector<std::pair<T, std::uint64_t>> data(m_size);

        std::vector<T> norm_probs(m_size);
        std::vector<std::int64_t> large_block(m_size);
        std::vector<std::int64_t> small_block(m_size);

        std::int64_t num_small_block = 0;
        std::int64_t num_large_block = 0;

        const auto sum = std::reduce(std::execution::par_unseq, weights.begin(), weights.end(), T { 0.0 });
        const auto scale_factor = m_size / sum;
        std::transform(std::execution::par_unseq, weights.cbegin(), weights.cend(), norm_probs.begin(), [=](const auto w) -> T { return w * scale_factor; });

        for (std::int64_t i = m_size - 1; i >= 0; i--) {
            if (norm_probs[i] < T { 1.0 }) {
                small_block[num_small_block++] = i;
            } else {
                large_block[num_large_block++] = i;
            }
        }

        while (num_small_block && num_large_block) {
            const auto cur_small_block = small_block[--num_small_block];
            const auto cur_large_block = large_block[--num_large_block];

            data[cur_small_block] = std::make_pair(norm_probs[cur_small_block], cur_large_block);

            norm_probs[cur_large_block] = norm_probs[cur_large_block] + norm_probs[cur_small_block] - 1;
            if (norm_probs[cur_large_block] < 1) {
                small_block[num_small_block++] = cur_large_block;
            } else {
                large_block[num_large_block++] = cur_large_block;
            }
        }

        while (num_large_block) {
            data[large_block[--num_large_block]].first = 1;
        }

        while (num_small_block) {
            data[small_block[--num_small_block]].first = 1;
        }
        return data;
    }

private:
    std::vector<std::pair<T, std::uint64_t>> m_data;
};

/**
 * @brief Class for sampling of a specter of values according to a distribution.
 */
template <Floating T = double>
class SpecterDistribution {
public:
    SpecterDistribution()
    {
    }

    /**
     * @brief Initialize specter distribution
     * @param energy A vector energy in keV for each bin. Must be monotonical increasing
     * @param weighs A vector of weights for each bin. Weights should be at least of size two. The weights vector will be normalized such that sum of weights is unity. Size of weights must be equal to size of energies.
     */
    SpecterDistribution(const std::vector<T>& energy, const std::vector<T>& weights)
    {
        if (energy.size() != weights.size() || energy.size() == 0)
            return;
        const auto table = RandomDistribution<T>::generateTable(weights);

        m_data.resize(table.size());
        std::transform(std::execution::par_unseq, table.cbegin(), table.cend(), energy.cbegin(), m_data.begin(), [](const auto& tab, const T energy) {
            return DataElement { .alias = tab.second, .prob = tab.first, .energy = energy };
        });
    }

    /**
     * @brief Initialize specter distribution
     * @param weights A vector of pairs {energy, probability} for each bin. Weights should be at least of size two. The weights vector will be normalized such that sum of weights is unity.
     */
    SpecterDistribution(const std::vector<std::pair<T, T>>& specter)
    {
        std::vector<T> weights(specter.size());
        std::vector<T> energies(specter.size());

        std::transform(std::execution::par_unseq, specter.cbegin(), specter.cend(), energies.begin(), [](const auto& p) { return p.first; });
        std::transform(std::execution::par_unseq, specter.cbegin(), specter.cend(), weights.begin(), [](const auto& p) { return p.second; });
        const auto table = RandomDistribution<T>::generateTable(weights);

        m_data.resize(table.size());
        std::transform(std::execution::par_unseq, table.cbegin(), table.cend(), energies.cbegin(), m_data.begin(), [](const auto& tab, const T energy) {
            return DataElement { .alias = tab.second, .prob = tab.first, .energy = energy };
        });
    }

    /**
     * @brief Sample an index from 0 to size() - 1 according to weights. This function is thread safe.
     * @param state Random state for sampling of an index.
     * @return index with probability according to weights[index].
     */
    std::size_t sampleIndex(RandomState& state) const
    {
        const auto r = state.randomUniform<T>();
        const auto k = state.randomUniform<std::size_t>(m_data.size());
        return r < m_data[k].prob ? k : m_data[k].alias;
    }

    /**
     * @brief Sample an energy value according to weights probability. This function is thread safe.
     * The sampling is done by first randomly sample an index into the energy vector. A random uniform energy in the interval energies[sampleIndex] and energies[sampleIndex+1] is returned. If sampleIndex is the last index in weights, the last energy value is returned.
     * @return a random energy according to weights probability
     */
    T sampleValue(RandomState& state) const
    {
        const std::size_t ind = sampleIndex(state);
        return ind < m_data.size() - 1 ? state.randomUniform(m_data[ind].energy, m_data[ind + 1].energy) : m_data[ind].energy;
    }

protected:
private:
    struct DataElement {
        std::uint64_t alias = 0;
        T prob = 0;
        T energy = 0;
    };
    std::vector<DataElement> m_data;
};

template <Floating T, std::size_t N = 20>
class CPDFSampling {
public:
    CPDFSampling() = default;

    template <std::regular_invocable<T> F>
        requires std::is_same<std::invoke_result_t<F, T>, T>::value
    CPDFSampling(T xmin, T xmax, const F& pdf)
    {
        std::size_t n_start = N / 2;

        std::vector<T> points(n_start);
        for (std::size_t i = 0; i < n_start; ++i) {
            points[i] = xmin + i * (xmax - xmin) / (n_start - 1);
        }

        while (points.size() != N) {
            auto grid = buildGrid(points, pdf);
            std::vector<T> error(grid.size() - 1);
            for (std::size_t i = 0; i < grid.size() - 1; ++i) {
                const T est = integrateGrid(grid[i], grid[i + 1]);
                const auto inter = simpson_integral(points[i], points[i + 1], pdf);
                error[i] = std::abs(est - inter);
            }
            auto max = std::max_element(error.cbegin(), error.cend());
            const auto idx = std::distance(error.cbegin(), max);
            auto g0 = grid.cbegin() + idx;
            auto g1 = g0 + 1;
            const T new_point = g0->x + (g1->x - g0->x) * T { 0.5 };
            points.insert(points.cbegin() + idx + 1, new_point);
        }

        auto grid = buildGrid(points, pdf);
        std::copy(grid.cbegin(), grid.cend(), m_grid.begin());
    }
    T operator()(RandomState& state) const
    {
        const T r1 = state.randomUniform<T>(m_grid[0].e, m_grid[N - 1].e);
        auto pos1 = std::upper_bound(m_grid.cbegin() + 1, m_grid.cend() - 1, r1, [](const auto num, const auto& element) -> bool { return num < element.e; });
        auto pos0 = pos1 - 1;

        const auto v = r1 - pos0->e;
        const auto delta = pos1->e - pos0->e;

        const auto nom = (1 + pos0->a + pos0->b) * delta * v;
        const auto den = delta * delta + pos0->a * delta * v + pos0->b * v * v;
        const auto res = pos0->x + nom * (pos1->x - pos0->x) / den;
        return res;
    }
    T operator()(T max_x, RandomState& state) const
    {
        const auto max_xpos = std::upper_bound(m_grid.cbegin() + 1, m_grid.cend() - 1, max_x, [](const auto num, const auto& element) -> bool { return num < element.x; });
        T res;
        do {
            const T r1 = state.randomUniform<T>(m_grid[0].e, max_xpos->e);
            auto pos1 = std::upper_bound(m_grid.cbegin() + 1, max_xpos, r1, [](const auto num, const auto& element) -> bool { return num < element.e; });
            auto pos0 = pos1 - 1;

            const auto v = r1 - pos0->e;
            const auto delta = pos1->e - pos0->e;

            const auto nom = (1 + pos0->a + pos0->b) * delta * v;
            const auto den = delta * delta + pos0->a * delta * v + pos0->b * v * v;
            res = pos0->x + nom * (pos1->x - pos0->x) / den;
        } while (res > max_x);
        return res;
    }

protected:
    struct GridEl {
        T x = 0;
        T e = 0;
        T a = 0;
        T b = 0;
    };

    // calcluate p from points
    static T integrateGrid(const GridEl& g0, const GridEl& g1, std::size_t Nsteps = 50)
    {
        T sum = 0;
        const T step = (g1.x - g0.x) / (Nsteps - 1);
        for (std::size_t i = 1; i < Nsteps; ++i) {
            const T x = g0.x + step * i;
            const auto t = (x - g0.x) / (g1.x - g0.x);
            const auto a = g0.a;
            const auto b = g0.b;
            const auto p1 = 1 + a + b - a * t;
            const auto n = p1 / (2 * b * t) * (1 - std::sqrt(1 - (4 * b * t * t) / (p1 * p1)));

            const auto r1 = 1 + a * n + b * n * n;
            const auto r2 = (1 + a + b) * (1 - b * n * n);

            sum += r1 * r1 * (g1.e - g0.e) / (r2 * (g1.x - g0.x));
        }
        return sum * step;
    }
    template <typename F>
    static std::vector<GridEl> buildGrid(const std::vector<T> points, const F& pdf)
    {
        // cumulative prob (unnormalized)
        std::vector<T> e(points.size());
        e[0] = pdf(points[0]);
        for (std::size_t i = 1; i < e.size(); ++i) {
            const auto x0 = points[i - 1];
            const auto x1 = points[i];
            e[i] = e[i - 1] + simpson_integral(x0, x1, pdf);
        }

        std::vector<GridEl> grid;

        for (std::size_t i = 0; i < points.size() - 1; ++i) {

            const T ex = (e[i + 1] - e[i]) / (points[i + 1] - points[i]);
            const T b = T { 1 } - ex * ex / (pdf(points[i]) * pdf(points[i + 1]));
            const T a = ex / pdf(points[i]) - b - T { 1 };
            grid.push_back({ points[i], e[i], a, b });
        }
        grid.push_back({ points.back(), e.back(), T { 0 }, T { 0 } });

        return grid;
    }

    template <typename F, int Nsimp = 51>
    static T simpson_integral(T x0, T x1, const F& pdf)
    {
        const T h = (x1 - x0) / (Nsimp - 1);
        const T p1 = pdf(x0) + pdf(x1);
        T p2 = 0;
        for (std::size_t i = 1; i < Nsimp - 1; i += 2)
            p2 += pdf(x0 + i * h);
        T p3 = 0;
        for (std::size_t i = 2; i < Nsimp - 2; i += 2)
            p3 += pdf(x0 + i * h);
        return h * (p1 + 4 * p2 + 2 * p3) / T { 3 };
    }

private:
    std::array<GridEl, N> m_grid;
};

// class for numerical inverse transform of analytical probability density functions (pdfs)
//[[deprecated("Use CPDFSampling instead")]]
template <Floating T, int N = 20>
class RITA {
public:
    template <std::regular_invocable<T> F>
        requires std::is_same<std::invoke_result_t<F, T>, T>::value
    RITA(const T min, const T max, F pdf)
    {
        struct param {
            T x, e, a, b, error;
        };

        std::size_t n = 10;
        std::vector<param> v(n);
        v.reserve(N);
        for (std::size_t i = 0; i < n; ++i) {
            const T x = min + i * (max - min) / (n - 1);
            v[i].x = x;
            v[i].e = 0;
            v[i].error = -1;
        }
        // finding e
        for (std::size_t j = 1; j < n; ++j) {
            v[j].e = v[j - 1].e + simpson_integral(v[j - 1].x, v[j].x, pdf);
        }
        const T pdf_max_value = v[n - 1].e;

        for (std::size_t j = 1; j < n; ++j) {
            v[j].e = v[j].e / v[n - 1].e;
        }
        while (n != N) {
            for (std::size_t i = 0; i < n - 1; ++i) {
                if (v[i].error < 0) {
                    // finding a and b
                    const T temp = (v[i + 1].e - v[i].e) / (v[i + 1].x - v[i].x);
                    const T px0 = pdf(v[i].x) / pdf_max_value;
                    const T px1 = pdf(v[i + 1].x) / pdf_max_value;
                    if (px0 > 0 && px1 > 0) {
                        v[i].b = 1 - temp * temp / (px0 * px1);
                    } else {
                        v[i].b = 0;
                    }
                    if (px0 > 0) {
                        v[i].a = temp / px0 - v[i].b - 1;
                    } else {
                        v[i].a = 0;
                    }

                    // calculate erroe
                    v[i].error = 0;
                    for (std::size_t j = 1; j < 50; ++j) {
                        const T x = v[i].x + j * (v[i + 1].x - v[i].x) / 50;

                        const auto t = (x - v[i].x) / (v[i + 1].x - v[i].x);
                        const auto a = v[i].a;
                        const auto b = v[i].b;
                        const auto f = (1 + a + b - a * t);
                        const auto nn = f * (1 - std::sqrt(1 - 4 * b * t * t / (f * f))) / (2 * b * t);
                        const auto p1 = (1 + a * nn + b * nn * nn);
                        const auto p = p1 * p1 * (v[i + 1].e - v[i].e) / ((1 + a + b) * (1 - b * nn * nn) * (v[i + 1].x - v[i].x));
                        const auto res = std::abs(p - pdf(x));

                        v[i].error += res * (v[i + 1].x - v[i].x) / 50;
                    }
                }
            }

            auto max_iter = std::max_element(v.begin(), v.end(), [](const auto& lh, const auto& rh) { return lh.error < rh.error; });
            const T x0 = max_iter->x;
            max_iter->error = -1;
            ++max_iter;
            const T x1 = max_iter->x;
            param new_val = { x0 + (x1 - x0) / 2, 0, 0, 0, -1 };
            v.insert(max_iter, new_val);
            ++n;
            // finding e
            for (std::size_t j = 1; j < n; ++j) {
                v[j].e = v[j - 1].e + simpson_integral(v[j - 1].x, v[j].x, pdf);
            }
            for (std::size_t j = 1; j < n; ++j) {
                v[j].e = v[j].e / v[n - 1].e;
            }
        }

        std::transform(v.cbegin(), v.cend(), m_x.begin(), [](const auto& p) { return p.x; });
        std::transform(v.cbegin(), v.cend(), m_e.begin(), [](const auto& p) { return p.e; });
        std::transform(v.cbegin(), v.cend(), m_a.begin(), [](const auto& p) { return p.a; });
        std::transform(v.cbegin(), v.cend(), m_b.begin(), [](const auto& p) { return p.b; });
    }

    T operator()(RandomState& state) const
    {
        const auto r1 = state.randomUniform<T>();
        auto upper_bound = std::upper_bound(m_e.cbegin(), m_e.cend(), r1);
        const std::size_t index = std::distance(m_e.cbegin(), upper_bound);
        if (index == 0)
            return m_x[0];

        const auto v = r1 - m_e[index - 1];
        const auto d = m_e[index] - m_e[index - 1];
        return m_x[index - 1] + (1 + m_a[index - 1] + m_b[index - 1]) * d * v / (d * d + m_a[index - 1] * d * v + m_b[index - 1] * v * v) * (m_x[index] - m_x[index - 1]);
    }

    T operator()(RandomState& state, const T maxValue) const
    {
        // finding xmax
        auto upper_bound_value = std::upper_bound(m_x.cbegin(), m_x.cend(), maxValue);
        const std::size_t max_value_index = std::distance(m_x.cbegin(), upper_bound_value);
        const T modifier = upper_bound_value != m_x.cend() ? m_e[max_value_index] : 1;
        T res;
        do {
            const auto r1 = state.randomUniform<T>(modifier);
            auto upper_bound = std::upper_bound(m_e.cbegin(), m_e.cend(), r1);
            const auto index = std::distance(m_e.cbegin(), upper_bound) - 1;

            const auto v = r1 - m_e[index];
            const auto d = m_e[index + 1] - m_e[index];
            res = m_x[index] + (1 + m_a[index] + m_b[index]) * d * v / (d * d + m_a[index] * d * v + m_b[index] * v * v) * (m_x[index + 1] - m_x[index]);
        } while (res > maxValue);
        return res;
    }

protected:
    template <typename F>
    static T simpson_integral(const T start, const T stop, F pdf)
    {
        const T h = (stop - start) / 50;
        T result = pdf(start) + pdf(stop);
        for (std::size_t i = 1; i < 50; ++i) {
            const T prod = i % 2 == 0 ? 2 : 4; // prod is 2 when i is even else 4
            result += prod * pdf(start + h * i);
        }
        return h * result / 3;
    }

    T integral_p_bar(std::size_t index) const
    {
        const T bi = m_b[index];
        const T ai = m_a[index];
        const T xi = m_x[index];
        const T xii = m_x[index + 1];

        auto p = [=, *this](const T x) -> T {
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

private:
    std::array<T, N> m_x;
    std::array<T, N> m_e;
    std::array<T, N> m_b;
    std::array<T, N> m_a;
};
}