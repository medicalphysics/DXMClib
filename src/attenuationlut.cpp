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

#include "dxmc/attenuationlut.h"

#include <algorithm>
#include <cmath>
#include <string>


constexpr double KEV_TO_A = 12.398520;

void AttenuationLut::generate(const std::vector<Material>& materials, double minEnergy, double maxEnergy)
{
	m_minEnergy = minEnergy > 0.0 ? minEnergy : 0.0;
	m_maxEnergy = maxEnergy > m_minEnergy ? maxEnergy : m_minEnergy + 1.0;

	m_energyResolution = static_cast<std::size_t>(std::ceil((m_maxEnergy - m_minEnergy) / m_energyStep));
	m_nMaterials = materials.size();

	m_attData.resize(m_energyResolution * m_nMaterials * 4 + m_energyResolution);

	//generating energy table
	for (std::size_t i = 0; i < m_energyResolution; ++i)
		m_attData[i] = m_minEnergy + i * m_energyStep;

	for (std::size_t m = 0; m < m_nMaterials; ++m)
	{
		auto offset = m_energyResolution + m * m_energyResolution * 4;
		for (std::size_t i = 0; i < m_energyResolution; ++i)
			m_attData[i + offset] = materials[m].getTotalAttenuation(m_attData[i]);
		offset = 2 * m_energyResolution + m * m_energyResolution * 4;
		for (std::size_t i = 0; i < m_energyResolution; ++i)
			m_attData[i + offset] = materials[m].getPhotoelectricAttenuation(m_attData[i]);
		offset = 3 * m_energyResolution + m * m_energyResolution * 4;
		for (std::size_t i = 0; i < m_energyResolution; ++i)
			m_attData[i + offset] = materials[m].getComptonAttenuation(m_attData[i]);
		offset = 4 * m_energyResolution + m * m_energyResolution * 4;
		for (std::size_t i = 0; i < m_energyResolution; ++i)
			m_attData[i + offset] = materials[m].getRayleightAttenuation(m_attData[i]);
	}
	generateFFdata(materials);
	generateSFdata(materials);

	//maxenergies
	std::vector<double> densArray;
	std::vector<unsigned char> matArray;
	for (std::size_t i = 0; i < materials.size(); ++i)
	{
		densArray.push_back(materials[i].standardDensity());
		matArray.push_back(static_cast<unsigned char>(i));
	}
	generateMaxMassTotalAttenuation(matArray.begin(), matArray.end(), densArray.begin());
}



template<typename T >
inline double interp(const T x1, const T x2, const T y1, const T y2, T xres)
{
	return y1 + (y2 - y1) * (xres - x1) / (x2 - x1);
}

template<typename T >
inline double interp(const T x[2], const T y[2], T xres)
{
	const T val = y[0] + (y[1] - y[0]) * (xres - x[0]) / (x[1] - x[0]);
	return xres < x[0] ? y[0] : xres > x[1] ? y[1] : val;
}

template<typename It, typename T >
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

	return interp( *lower, *upper, *lowery, *uppery, xvalue);
}

double AttenuationLut::totalAttenuation(std::size_t material, double energy) const
{
	const std::size_t offset = m_energyResolution + m_energyResolution * material * 4;
	if (energy <= m_minEnergy)
		return m_attData[offset];
	else if (energy >= m_maxEnergy)
		return m_attData[offset + m_energyResolution - 1];

	const std::size_t idx = static_cast<std::size_t>((energy - m_minEnergy) / m_energyStep);
	return interp(&m_attData[idx], &m_attData[idx + offset], energy);
}
 double AttenuationLut::rayleightAttenuation(std::size_t material, double energy) const
{
	 const std::size_t offset = m_energyResolution + m_energyResolution * material * 4 + m_energyResolution * 3;
	if (energy <= m_minEnergy)
		return m_attData[offset];
	else if (energy >= m_maxEnergy)
		return m_attData[offset + m_energyResolution - 1];

	const std::size_t idx = static_cast<std::size_t>((energy - m_minEnergy) / m_energyStep);
	return interp(&m_attData[idx], &m_attData[idx + offset], energy);
}
double AttenuationLut::photoelectricAttenuation(std::size_t material, double energy) const 
{
	const std::size_t offset = m_energyResolution + m_energyResolution * material * 4 + m_energyResolution;
	if (energy <= m_minEnergy)
		return m_attData[offset];
	else if (energy >= m_maxEnergy)
		return m_attData[offset + m_energyResolution - 1];

	const std::size_t idx = static_cast<std::size_t>((energy - m_minEnergy) / m_energyStep);
	return interp(&m_attData[idx], &m_attData[idx + offset], energy);
}
double AttenuationLut::comptonAttenuation(std::size_t material, double energy) const {
	const std::size_t offset = m_energyResolution + m_energyResolution * material * 4 + m_energyResolution * 2;
	if (energy <= m_minEnergy)
		return m_attData[offset];
	else if (energy >= m_maxEnergy)
		return m_attData[offset + m_energyResolution - 1];

	const std::size_t idx = static_cast<std::size_t>((energy - m_minEnergy) / m_energyStep);
	return interp(&m_attData[idx], &m_attData[idx + offset], energy);
}

std::array<double, 3> AttenuationLut::photoComptRayAttenuation(std::size_t material, double energy) const
{
	const std::size_t offset = m_energyResolution + m_energyResolution * material * 4 + m_energyResolution;
	std::array<double, 3> att;
	if (energy <= m_minEnergy)
	{
		for (std::size_t i = 0; i < 3; ++i)
		{
			att[i] = m_attData[offset + m_energyResolution * i];
		}
	}
	else if(energy >= m_maxEnergy)
	{
		for (std::size_t i = 0; i < 3; ++i)
		{
			att[i] = m_attData[offset + m_energyResolution * (i + 1) - 1];
		}
	}
	else {
		const std::size_t idx = static_cast<std::size_t>((energy - m_minEnergy) / m_energyStep);
		for (std::size_t i = 0; i < 3; ++i)
		{
			att[i] = interp(&m_attData[idx], &m_attData[offset + m_energyResolution * i + idx], energy);
		}
	}
	return att;
}


double AttenuationLut::maxMassTotalAttenuation(double energy) const
{
	if (energy <= m_minEnergy)
		return m_maxMassAtt[0];
	else if (energy >= m_maxEnergy)
		return m_maxMassAtt[m_energyResolution - 1];
	const std::size_t idx = static_cast<std::size_t>((energy - m_minEnergy) / m_energyStep);
	return interp(&m_attData[idx], &m_maxMassAtt[idx], energy);
}
double AttenuationLut::momentumTransfer(std::size_t material, double cumFormFactorSquared) const
{
	const auto offset = m_momtranfResolution + material * m_momtranfResolution;
	const auto begin = m_rayleighFormFactorSqr.cbegin() + offset;
	const auto end = begin + m_momtranfResolution;

	const auto pos = std::lower_bound(begin, end, cumFormFactorSquared);
	if (pos == end)
		return m_rayleighFormFactorSqr[m_momtranfResolution - 1];
	const auto distance = std::distance(begin, pos);
	if (distance <= 0)
		return m_rayleighFormFactorSqr[0];

	return interp(&m_rayleighFormFactorSqr[offset + distance - 1], &m_rayleighFormFactorSqr[distance - 1], cumFormFactorSquared);
}

double AttenuationLut::cumFormFactorSquared(std::size_t material, double momentumTransfer) const
{
	const auto offset = m_momtranfResolution + material * m_momtranfResolution;
	auto begin = m_rayleighFormFactorSqr.cbegin();
	auto end = begin + m_momtranfResolution;

	const auto pos = std::lower_bound(begin, end, momentumTransfer);

	if (pos == end)
		return m_rayleighFormFactorSqr[offset + m_momtranfResolution - 1];
	
	const auto distance = std::distance(begin, pos);
	if (distance <= 0)
		return m_rayleighFormFactorSqr[offset];
	return interp(&m_rayleighFormFactorSqr[distance-1], &m_rayleighFormFactorSqr[offset + distance-1], momentumTransfer);
}

double AttenuationLut::comptonScatterFactor(std::size_t material, double momentumTransfer) const
{
	const std::size_t offset = m_momtranfResolution + m_momtranfResolution * material;
	if (momentumTransfer <= 0.0)
		return m_comptonScatterFactor[offset];
	
	const std::size_t idx = static_cast<std::size_t>(momentumTransfer / m_momtranfStep);
	if (idx >= m_momtranfResolution)
		return m_comptonScatterFactor[offset + m_momtranfResolution - 1];
	return interp(&m_comptonScatterFactor[idx], &m_comptonScatterFactor[idx + offset], momentumTransfer);
}
double AttenuationLut::momentumTransfer(double energy, double angle)
{
	constexpr double k = 1.0 / KEV_TO_A; //constant for momentum transfer
	return energy * std::sin(angle / 2.0) * k;
}
double AttenuationLut::momentumTransferFromCos(double energy, double cosAngle)
{
	constexpr double k = 1.0 / KEV_TO_A; //constant for momentum transfer
	return energy * k * std::sqrt(0.5 - cosAngle * 0.5);
}
double AttenuationLut::momentumTransferMax(double energy)
{
	constexpr double k = 1.0 / KEV_TO_A; //constant for momentum transfer
	return energy * k;
}

double AttenuationLut::cosAngle(double energy, double momentumTransfer)
{
	const double invE = 1.0 / energy;
	return 1.0 - 2.0 * momentumTransfer * momentumTransfer * KEV_TO_A * KEV_TO_A * invE * invE;
}


std::vector<double> trapz(const std::vector<double>& f, const std::vector<double>& x)
{
	std::vector<double> integ(f.size(), 0.0);
	integ[0] = 0.0;
	for (std::size_t i = 1; i < f.size(); ++i)
	{
		integ[i] = integ[i - 1] + (f[i - 1] + f[i]) * 0.5 * (x[i] - x[i - 1]);
	}
	return integ;
}


void AttenuationLut::generateFFdata(const std::vector<Material>& materials)
{
	//se http://rcwww.kek.jp/research/egs/egs5_manual/slac730-150228.pdf

	const double globalMomMax = momentumTransferMax(m_maxEnergy);
	
	std::vector<double> qvalues(m_momtranfResolution);
	for (std::size_t i = 0; i < m_momtranfResolution; ++i)
		qvalues[i] = (i*globalMomMax/(m_momtranfResolution-1))*(i * globalMomMax / (m_momtranfResolution-1))/ globalMomMax;  // nonlinear integration steps


	m_rayleighFormFactorSqr.resize(m_momtranfResolution + materials.size() * m_momtranfResolution);
	std::copy(qvalues.cbegin(), qvalues.cend(), m_rayleighFormFactorSqr.begin());

	for (std::size_t i = 0; i < materials.size(); ++i)
	{
		auto ffmaxsqr = materials[i].getRayleightFormFactorSquared(qvalues);
		auto ffmaxsqrint = trapz(ffmaxsqr, qvalues);
		auto start = m_rayleighFormFactorSqr.begin() + m_momtranfResolution + i * m_momtranfResolution;
		std::copy(ffmaxsqrint.cbegin(), ffmaxsqrint.cend(), start);
	}
}

void AttenuationLut::generateSFdata(const std::vector<Material>& materials)
{
	// Compton scatter corrections

	//se http://rcwww.kek.jp/research/egs/egs5_manual/slac730-150228.pdf

	const double globalMomMax = momentumTransferMax(m_maxEnergy);
	std::vector<double> qvalues(m_momtranfResolution);
	m_momtranfStep = globalMomMax / (m_momtranfResolution - 1);
	for (std::size_t i = 0; i < m_momtranfResolution; ++i)
		qvalues[i] = i * m_momtranfStep;

	m_comptonScatterFactor.resize(m_momtranfResolution + materials.size() * m_momtranfResolution);
	std::copy(qvalues.cbegin(), qvalues.cend(), m_comptonScatterFactor.begin());
	for (std::size_t i = 0; i < materials.size(); ++i)
	{
		auto sfmax = materials[i].getComptonNormalizedScatterFactor(qvalues);
		auto start = m_comptonScatterFactor.begin() + m_momtranfResolution + i * m_momtranfResolution;
		std::copy(sfmax.cbegin(), sfmax.cend(), start);
	}
}
