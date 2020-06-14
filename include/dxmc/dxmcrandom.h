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

#include <cstdint>
#include <limits>
#include <utility>
#include <vector>
#include <random>

// The function below is borrowed from:
// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)

class RandomState {
public:
    RandomState()
    {
        std::random_device d;
        m_state[0] = static_cast<std::uint64_t>(d());
        m_state[1] = static_cast<std::uint64_t>(d());
    }
    RandomState(std::uint64_t* state)
    {
        m_state[0] = state[0];
        m_state[1] = state[1];
    }
    RandomState(const RandomState&) = delete; // non construction-copyable
    RandomState& operator=(const RandomState&) = delete; // non copyable
   
    std::uint64_t m_state[2];
};


inline std::uint32_t pcg32_random_r(RandomState &state)
{
    std::uint64_t oldstate = state.m_state[0];
    // Advance internal state
    state.m_state[0] = oldstate * 6364136223846793005ULL + (state.m_state[1] | 1);
    // Calculate output function (XSH RR), uses old state for max ILP
    std::uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    std::uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

template <typename T>
inline T randomUniform(RandomState &state) noexcept
{
    static_assert(std::is_floating_point<T>::value, "Uniform random number requires floating point precision");
    const T r = static_cast<T>(pcg32_random_r(state));
    return r / (std::numeric_limits<std::uint32_t>::max() - 1);
}

/*
see this: http://prng.di.unimi.it/

inline std::uint64_t xoroshiro128plus(RandomState &state) noexcept
{
    std::uint64_t s0 = state.m_state[0];
    std::uint64_t s1 = state.m_state[1];
    std::uint64_t result = s0 + s1;
    s1 ^= s0;
    state.m_state[0] = ((s0 << 55) | (s0 >> 9)) ^ s1 ^ (s1 << 14);
    state.m_state[1] = (s1 << 36) | (s1 >> 28);
    return result;
}

template <typename T>
inline T randomUniform(RandomState &state) noexcept
{
    static_assert(std::is_floating_point<T>::value, "Uniform random number requires floating point precision");
    const T r = static_cast<T>(xoroshiro128plus(state));
    return r / (std::numeric_limits<std::uint64_t>::max() - 1);
}
*/
template <typename T>
inline T randomUniform(RandomState &state, const T max) noexcept
{
    const T r = randomUniform<T>(state);
    return r * max;
}

template <typename T>
inline T randomUniform(RandomState &state, const T min, const T max) noexcept
{
    const T r = randomUniform<T>(state);
    const T range = max - min;
    return min + r * range;
}

void randomSeed(RandomState &state);

class RandomDistribution {
public:
    RandomDistribution(const std::vector<double>& weights);
    std::size_t sampleIndex();
    std::size_t sampleIndex(RandomState &state) const;

protected:
    RandomState m_state;
    void generateTable(const std::vector<double>& weights);

private:
    std::vector<std::uint64_t> m_alias;
    std::vector<double> m_probs;
    std::size_t m_size = 0;
};

class SpecterDistribution : public RandomDistribution {
public:
    SpecterDistribution(const std::vector<double>& weights, const std::vector<double>& energies);
    double sampleValue();
    double sampleValue(RandomState &state) const; // thread safe
private:
    std::vector<double> m_energies;
};
