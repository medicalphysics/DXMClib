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

double AttenuationLut::momentumTransferSquared(std::size_t material, double cumFormFactorSquared) const
{
	auto begin = m_coherData.begin();
	std::advance(begin, m_energyResolution + m_energyResolution * material);
	auto end = m_coherData.begin();
	std::advance(end, m_energyResolution + m_energyResolution * (material+1));

	auto beginy = m_coherData.begin();
	auto endy = m_coherData.begin();
	std::advance(endy, m_energyResolution);
	return interpolate(begin, end, beginy, endy, cumFormFactorSquared);
}

double AttenuationLut::cumFormFactorSquared(std::size_t material, double momentumTransferSquared) const
{
	const std::size_t idx = static_cast<std::size_t>(momentumTransferSquared / m_momtMaxSqr);
	if (idx >= m_energyResolution-1)
		return m_coherData[m_energyResolution + m_energyResolution * (material + 1) - 1];
	return interp(&m_coherData[idx], &m_coherData[idx + m_energyResolution + m_energyResolution * material], momentumTransferSquared);
}
double AttenuationLut::momentumTransfer(double energy, double angle)
{
	constexpr double k = 0.0806554; //constant for momentum transfer
	return energy * std::sin(angle / 2.0) * k;
}
double AttenuationLut::momentumTransferMax(double energy)
{
	constexpr double k = 0.0806554; //constant for momentum transfer
	return energy * k;
}


void AttenuationLut::generateFFdata(const std::vector<Material>& materials)
{

	//se http://rcwww.kek.jp/research/egs/egs5_manual/slac730-150228.pdf

	auto momtMax = momentumTransferMax(m_maxEnergy);
	m_momtMaxSqr = momtMax * momtMax;
	m_momtStepSqr = m_momtMaxSqr / (m_energyResolution - 1);


	m_coherData.resize(m_energyResolution + m_energyResolution * m_nMaterials);
	std::fill(m_coherData.begin(), m_coherData.end(), 0.0);
	for (std::size_t i = 0; i < m_energyResolution; ++i)
		m_coherData[i] = m_momtStepSqr * i;


	std::size_t integratorResolution = 1;
	while (integratorResolution * m_energyResolution < 512)
		++integratorResolution;

	std::vector<double> integratorX(integratorResolution*m_energyResolution); // momentum transfer
	std::vector<double> integratorXSQR(integratorX.size()); // momentum transfer squared
	const double integratorStep = momtMax / (integratorX.size() - 1);
	const double integratorStepSqr = m_momtMaxSqr / (integratorX.size() - 1);
	for (std::size_t i = 0; i < integratorX.size(); ++i)
	{
		integratorX[i] = i * integratorStep;
		integratorXSQR[i] = i * integratorStepSqr;
	}
	const double integratorStepSQR = integratorXSQR[1] - integratorXSQR[0];

	for (std::size_t m = 0; m < m_nMaterials; ++m)
	{
		const auto integratorY = materials[m].getFormFactorSquared(integratorX); // squared form factor

		const auto offset = m_energyResolution + m_energyResolution * m ;

		double integral = 0;
		for (std::size_t i = 1; i < integratorResolution*m_energyResolution; ++i)
		{
			integral = integral + integratorStepSQR * (integratorY[i - 1] + integratorY[i]) / 2.0;
			if (i % integratorResolution == 0)
			{
				m_coherData[i / integratorResolution + offset] = integral;
			}
		}
	}

}
