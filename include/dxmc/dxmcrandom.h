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

#include <limits>
#include <random>
#include <vector>

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
        const T r = static_cast<T>(pcg32());
        return r / (std::numeric_limits<std::uint32_t>::max() - 1);
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
        } else
            static_assert(std::is_integral<T>::value || std::is_floating_point<T>::value, "Must be integral or floating point value.");
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
class RandomDistribution {
public:
    /**
     * @brief Initialize RandomDistribution class
     * @param weights A vector of probabilities for each bin. Weights should be atleast of size two. The weights vector will be normalized such that sum of weights is unity.
    */
    RandomDistribution(const std::vector<double>& weights);
    /**
     * @brief Sample an index from 0 to size() - 1 according to weights.
     * @return index with probability according to weights[index].
    */
    std::size_t sampleIndex();
    /**
     * @brief Sample an index from 0 to size() - 1 according to weights. This function is thread safe.
     * @param state Random state for sampling of an index.
     * @return index with probability according to weights[index].
    */
    std::size_t sampleIndex(RandomState& state) const;
    /**
     * @brief Size of the weights vector
     * @return Size of weights vector.
    */
    std::size_t size() const { return m_size; }

protected:
    RandomState m_state;
    void generateTable(const std::vector<double>& weights);

private:
    std::vector<std::uint64_t> m_alias;
    std::vector<double> m_probs;
    std::size_t m_size = 0;
};

/**
 * @brief Class for sampling of a specter of values according to a distribution.
*/
class SpecterDistribution : public RandomDistribution {
public:
    /**
     * @brief Initialize specter distribution
     * @param weights A vector of probabilities for each bin. Weights should be atleast of size two. The weights vector will be normalized such that sum of weights is unity.
     * @param energies A vector of values that are sampled according to weights. Values must be monotomic increasing. Lenght of energies must be equal to lenght of weights.
    */
    SpecterDistribution(const std::vector<double>& weights, const std::vector<double>& energies);
    /**
     * @brief Sample an energy value according to weights probabiliy. 
     * The sampling is done by first randomly sample an index into the energy vector. A random uniform energy in the interval energies[sampleIndex] and energies[sampleIndex+1] is returned. If sampleIndex is the last index in weights, the last energy value is returned.
     * @return a random energy according to weights probability
    */
    double sampleValue();
    /**
     * @brief Sample an energy value according to weights probabiliy. This function is thread safe.
     * The sampling is done by first randomly sample an index into the energy vector. A random uniform energy in the interval energies[sampleIndex] and energies[sampleIndex+1] is returned. If sampleIndex is the last index in weights, the last energy value is returned.
     * @return a random energy according to weights probability
    */
    double sampleValue(RandomState& state) const;

private:
    std::vector<double> m_energies;
};
