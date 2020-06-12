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

// The function below is borrowed from:
// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)

inline std::uint32_t pcg32_random_r(std::uint64_t s[2])
{
    std::uint64_t oldstate = s[0];
    // Advance internal state
    s[0] = oldstate * 6364136223846793005ULL + (s[1] | 1);
    // Calculate output function (XSH RR), uses old state for max ILP
    std::uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
    std::uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

template <typename T>
inline T randomUniform(std::uint64_t s[2]) noexcept
{
    static_assert(std::is_floating_point<T>::value, "Uniform random number requires floating point precision");
    const T r = static_cast<T>(pcg32_random_r(s));
    return r / (std::numeric_limits<std::uint32_t>::max()-1);
}

/*
see this: http://prng.di.unimi.it/

inline std::uint64_t xoroshiro128plus(std::uint64_t s[2]) noexcept
{
    std::uint64_t s0 = s[0];
    std::uint64_t s1 = s[1];
    std::uint64_t result = s0 + s1;
    s1 ^= s0;
    s[0] = ((s0 << 55) | (s0 >> 9)) ^ s1 ^ (s1 << 14);
    s[1] = (s1 << 36) | (s1 >> 28);
    return result;
}

template <typename T>
inline T randomUniform(std::uint64_t s[2]) noexcept
{
    static_assert(std::is_floating_point<T>::value, "Uniform random number requires floating point precision");
    const T r = static_cast<T>(xoroshiro128plus(s));
    return r / (std::numeric_limits<std::uint64_t>::max() - 1);
}
*/
template <typename T>
inline T randomUniform(std::uint64_t s[2], const T max) noexcept
{
    const T r = randomUniform<T>(s);
    return r * max;
}

template <typename T>
inline T randomUniform(std::uint64_t s[2], const T min, const T max) noexcept
{
    const T r = randomUniform<T>(s);
    const T range = max - min;
    return min + r * range;
}

void randomSeed(std::uint64_t s[2]);

class RandomDistribution {
public:
    RandomDistribution(const std::vector<double>& weights);
    std::size_t sampleIndex();
    std::size_t sampleIndex(std::uint64_t seed[2]) const;

protected:
    std::uint64_t m_seed[2];
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
    double sampleValue(std::uint64_t seed[2]) const; // thread safe
private:
    std::vector<double> m_energies;
};
