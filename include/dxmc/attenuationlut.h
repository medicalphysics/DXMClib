/*This file is part of OpenDXMC.

OpenDXMC is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

OpenDXMC is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with OpenDXMC. If not, see < https://www.gnu.org/licenses/>.

Copyright 2019 Erlend Andersen
*/


#pragma once

#include "dxmc/material.h"
#include <vector>
#include <array>


/**
 * @brief Attenuation look up tables. 
 * 
 * Class for generation and fast access of total and partial attenuation coefficients for materials specified.
 */
class AttenuationLut
{
public:
	/**
	 * @brief Construct a new Attenuation Lut object
	 * 
	 */
	AttenuationLut() {};
	
	/**
	 * @brief Energy resolution of the lookup table
	 * 
	 * Sets the energy resolution of the lookup table, minimum value is 0.1 keV
	 * @param keV Energy resolution in keV, default value is 1.0 keV
	 */
	void setEnergyResolution(double keV) { m_energyStep = keV > 0.1 ? keV : 0.1; }

	/**
	 * @brief Energy resolution in keV
	 * 
	 * @return double 
	 */
	double energyResolution()const { return m_energyStep; }
	
	/**
	 * @brief Generate lookup tables for attenuation data
	 * The order of materials matters, i.e the first material will get index 0, second has index 1 etc...
	 * @param materials Vector of materials to generate attenuation table
	 * @param minEnergy Minimum energy to consider in keV
	 * @param maxEnergy Maximum energy to consider in keV
	 */
	void generate(const std::vector<Material>& materials, double minEnergy=0.0, double maxEnergy = 150.0);
	
	/**
	 * @brief Generate lookup tables for attenuation data
	 * The order of materials matters, i.e the first material will get index 0, second has index 1 etc...
	 * @param materials Vector of materials to generate attenuation table
	 * @param energies Vector of energies to generate in lookup table (in keV)
	 */
	void generate(const std::vector<Material>& materials, const std::vector<double>& energies);
	
	/**
	 * @brief Return the maximum mass attenuation value for all materials at spesific photon energy
	 * 
	 * @param energy Photon energy for maximum mass attenuation in table 
	 * @return double 
	 */
	double maxMassTotalAttenuation(double energy) const;
	
	/**
	 * @brief Calculate max total mass attenuation for a material and density array
	 * This templated function will iterate over material indices and density indices and calculate maximum total mass attenuation value for each photon energy in the lookup table
	 * 
	 * @param materialIndexBegin Pointer or iterator pointing to the first index in the material array
	 * @param materialIndexEnd Pointer or iterator pointing to the last + 1 index in the material array
	 * @param densityBegin Pointer or iterator pointing to the first index in the density array
	 */
	template<typename It1, typename It2>
	void generateMaxMassTotalAttenuation(It1 materialIndexBegin, It1 materialIndexEnd, It2 densityBegin);

	/**
	 * @brief Total mass attenuation for a material at specified photon energy
	 * 
	 * @param material material index
	 * @param energy photon energy in keV
	 * @return double 
	 */
	double totalAttenuation(std::size_t material, double energy) const;

	/**
	 * @brief Rayleigh mass attenuation for a material at specified photon energy
	 * 
	 * @param material material index
	 * @param energy photon energy in keV
	 * @return double 
	 */
	double rayleightAttenuation(std::size_t material, double energy) const;

	/**
	 * @brief Photoelectric mass attenuation for a material at specified photon energy
	 * 
	 * @param material material index
	 * @param energy photon energy in keV
	 * @return double 
	 */
	double photoelectricAttenuation(std::size_t material, double energy) const;

	/**
	 * @brief Compton mass attenuation for a material at specified photon energy
	 * 
	 * @param material material index
	 * @param energy photon energy in keV
	 * @return double 
	 */
	double comptonAttenuation(std::size_t material, double energy) const;

	/**
	 * @brief Shortcut to returning photoelectric, compton and rayleigh, respectivly, mass attenuation for a material at specified photon energy
	 * 
	 * @param material material index
	 * @param energy photon energy in keV
	 * @return std::array<double, 3>
	 */
	std::array<double, 3> photoComptRayAttenuation(std::size_t material, double energy) const;

	/**
	 * @brief Iterator to the first element of photon energy data
	 * 
	 * @return std::vector<double>::iterator 
	 */
	std::vector<double>::iterator energyBegin() { return m_attData.begin(); }

	/**
	 * @brief Iterator to the end of photon energy data
	 * 
	 * @return std::vector<double>::iterator 
	 */
	std::vector<double>::iterator energyEnd() { return m_attData.begin() + m_energyResolution; }

	/**
	 * @brief Iterator to the first element of total attenuation data for specified material
	 * 
	 * @param material Material index 
	 * @return std::vector<double>::iterator 
	 */
	std::vector<double>::iterator attenuationTotalBegin(std::size_t material) { return m_attData.begin() + m_energyResolution + m_energyResolution * 4 * material; } // todo error checking
	
	/**
	 * @brief Iterator to the end of total attenuation data for specified material
	 * 
	 * @param material Material index 
	 * @return std::vector<double>::iterator 
	 */
	std::vector<double>::iterator attenuationTotalEnd(std::size_t material) { return m_attData.begin() + m_energyResolution + m_energyResolution * 4 * material + m_energyResolution; } // todo error checking

	/**
	 * @brief Calculate squared momentumtransfer for a material at a cummulative squared atomic form factor
	 * 
	 * @param material Material index
	 * @param cumFormFactorSquared cummulative atomic form factor
	 * @return double 
	 */
	double momentumTransferSquared(std::size_t material, double cumFormFactorSquared) const;

	/**
	 * @brief Calculate cummulative atomic form factor squared from the squared momentum transfer for a specified material
	 * 
	 * @param material Material index
	 * @param momentumTransferSquared squared of momentum transfer 
	 * @return double 
	 */
	double cumFormFactorSquared(std::size_t material, double momentumTransferSquared) const;
	
	/**
	 * @brief Momentum transfer for Rayleigh scattering
	 * 
	 * @param energy Photon energy in keV
	 * @param angle Scattering angle in radians
	 * @return double 
	 */
	static double momentumTransfer(double energy, double angle);

	/**
	 * @brief Max of momentum transfer for a Rayleigh scattering 
	 * 
	 * @param energy Photon energy in keV
	 * @return double 
	 */
	static double momentumTransferMax(double energy);

protected:
	/**
	 * @brief Generate atomic form factors for materials specified.
	 * 
	 * @param materials Vector of Materials
	 */
	void generateFFdata(const std::vector<Material>& materials);
private:
	double m_minEnergy = 0;
	double m_maxEnergy = 150.0;
	double m_energyStep = 1.0;
	double m_momtMaxSqr = 0;
	double m_momtStepSqr = 0;
	std::size_t m_energyResolution = 150;
	std::size_t m_materials = 0;
	std::vector<double> m_attData; // energy, array-> total, photo, compton, rauleight
	std::vector<double> m_coherData; //qsquared, array-> A(qsquared)
	std::vector<double> m_maxMassAtt;
};



template<typename It1, typename It2>
void AttenuationLut::generateMaxMassTotalAttenuation(It1 materialIndexBegin, It1 materialIndexEnd, It2 densityBegin)
{
	std::vector<double> maxDens(m_materials, 0.0);
	
	while (materialIndexBegin != materialIndexEnd)
	{
		maxDens[*materialIndexBegin] = std::max(maxDens[*materialIndexBegin], *densityBegin);
		++materialIndexBegin;
		++densityBegin;
	}

	m_maxMassAtt.resize(m_energyResolution);
	std::fill(m_maxMassAtt.begin(), m_maxMassAtt.end(), 0.0);
	for (std::size_t material = 0; material < m_materials; ++material)
	{
		for (std::size_t i = 0; i < m_energyResolution; ++i)
		{
			//m_maxMassAtt[i] = std::max(m_maxMassAtt[i], totalAttenuation(material, m_coherData[i])*maxDens[material]);
			m_maxMassAtt[i] = std::max(m_maxMassAtt[i], totalAttenuation(material, m_attData[i])*maxDens[material]);
		}
	}
}