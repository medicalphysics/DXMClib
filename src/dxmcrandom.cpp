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

#include "dxmc/dxmcrandom.h"
#include <numeric>

void randomSeed(RandomState &state)
{
    std::random_device d;
    state.m_state[0] = static_cast<std::uint64_t>(d()) + static_cast<std::uint64_t>(d());
    state.m_state[1] = static_cast<std::uint64_t>(d()) + static_cast<std::uint64_t>(d());

    while (state.m_state[0] == 0) {
        state.m_state[0] = d();
    }
    while (state.m_state[1] == 0) {
        state.m_state[1] = d();
    }
}

RandomDistribution::RandomDistribution(const std::vector<double>& weights)
{
    randomSeed(m_state); // initialize prng
    generateTable(weights);
}

void RandomDistribution::generateTable(const std::vector<double>& weights)
{
    // Squaring the histogram method
    m_size = weights.size();
    m_probs.resize(m_size);
    m_alias.resize(m_size);
    std::vector<double> norm_probs(m_size);
    std::vector<std::int64_t> large_block(m_size);
    std::vector<std::int64_t> small_block(m_size);

    std::int64_t num_small_block = 0;
    std::int64_t num_large_block = 0;

    const double sum = std::accumulate(weights.begin(), weights.end(), 0.0);
    const double scale_factor = m_size / sum;
    std::transform(weights.cbegin(), weights.cend(), norm_probs.begin(), [=](const double w) -> double { return w * scale_factor; });

    for (std::int64_t i = m_size - 1; i >= 0; i--) {
        if (norm_probs[i] < 1.0) {
            small_block[num_small_block++] = i;
        } else {
            large_block[num_large_block++] = i;
        }
    }

    while (num_small_block && num_large_block) {
        const auto cur_small_block = small_block[--num_small_block];
        const auto cur_large_block = large_block[--num_large_block];
        m_probs[cur_small_block] = norm_probs[cur_small_block];
        m_alias[cur_small_block] = cur_large_block;
        norm_probs[cur_large_block] = norm_probs[cur_large_block] + norm_probs[cur_small_block] - 1;
        if (norm_probs[cur_large_block] < 1) {
            small_block[num_small_block++] = cur_large_block;
        } else {
            large_block[num_large_block++] = cur_large_block;
        }
    }

    while (num_large_block) {
        m_probs[large_block[--num_large_block]] = 1;
    }

    while (num_small_block) {
        m_probs[small_block[--num_small_block]] = 1;
    }
}

std::size_t RandomDistribution::sampleIndex()
{
    return sampleIndex(m_state);
}

std::size_t RandomDistribution::sampleIndex(RandomState &state) const
{
    /*const double r1 = state.randomUniform<double>();
    const double r2 = state.randomUniform<double>();
    std::size_t k = static_cast<std::size_t>(m_size * r1);
    return r2 < m_probs[k] ? k : m_alias[k];
    */
    const double r = state.randomUniform<double>();
    const double k = state.randomUniform<std::size_t>(m_size);
    return r < m_probs[k] ? k : m_alias[k];
}

SpecterDistribution::SpecterDistribution(const std::vector<double>& weights, const std::vector<double>& energies)
    : RandomDistribution(weights)
{
    m_energies = energies;
}

double SpecterDistribution::sampleValue()
{
    return sampleValue(m_state);
}
double SpecterDistribution::sampleValue(RandomState &state) const
{
    const std::size_t ind = sampleIndex(state);
    return ind < m_energies.size() - 1 ? state.randomUniform(m_energies[ind], m_energies[ind + 1]) : m_energies[ind];
}
