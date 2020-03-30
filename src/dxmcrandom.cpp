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
#include <random>

void randomSeed(std::uint64_t s[2])
{
	std::random_device d;
    s[0] = static_cast<std::uint64_t>(d()) + static_cast<std::uint64_t>(d());
	s[1] = static_cast<std::uint64_t>(d()) + static_cast<std::uint64_t>(d());

	while (s[0] == 0)
	{
		s[0] = d();
	}
	while (s[1] == 0)
	{
		s[1] = d();
	}
}

RandomDistribution::RandomDistribution(const std::vector<double>& weights) 
{
	randomSeed(m_seed); // initialize prng
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
	std::transform(weights.cbegin(), weights.cend(), norm_probs.begin(), [=](const double w)->double {return w * scale_factor; });
	
	for (std::int64_t i = m_size - 1; i >= 0; i--) 
	{
		if (norm_probs[i] < 1.0)
		{
			small_block[num_small_block++] = i;
		}
		else
		{
			large_block[num_large_block++] = i;
		}
	}

	while (num_small_block && num_large_block) 
	{
		const auto cur_small_block = small_block[--num_small_block];
		const auto cur_large_block = large_block[--num_large_block];
		m_probs[cur_small_block] = norm_probs[cur_small_block];
		m_alias[cur_small_block] = cur_large_block;
		norm_probs[cur_large_block] = norm_probs[cur_large_block] + norm_probs[cur_small_block] - 1;
		if (norm_probs[cur_large_block] < 1)
		{
			small_block[num_small_block++] = cur_large_block;
		}
		else
		{
			large_block[num_large_block++] = cur_large_block;
		}
	}

	while (num_large_block)
	{
		m_probs[large_block[--num_large_block]] = 1;
	}

	while (num_small_block)
	{
		m_probs[small_block[--num_small_block]] = 1;
	}

}

std::size_t RandomDistribution::sampleIndex()
{
	return sampleIndex(m_seed);
}

std::size_t RandomDistribution::sampleIndex(std::uint64_t seed[2]) const
{
	const double r1 = randomUniform<double>(seed);
	const double r2 = randomUniform<double>(seed);
	std::size_t k = static_cast<std::size_t>(m_size * r1);
	return r2 < m_probs[k] ? k : m_alias[k];
}


SpecterDistribution::SpecterDistribution(const std::vector<double>& weights, const std::vector<double>& energies)
	: RandomDistribution(weights)
{
	m_energies = energies;
}


double SpecterDistribution::sampleValue()
{
	return sampleValue(m_seed);
}
double SpecterDistribution::sampleValue(std::uint64_t seed[2]) const
{
	const std::size_t ind = sampleIndex(seed);
	return ind < m_energies.size() - 1 ? randomUniform(seed, m_energies[ind], m_energies[ind + 1]) : m_energies[ind];
	
}
