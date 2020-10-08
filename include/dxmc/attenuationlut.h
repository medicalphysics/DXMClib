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

#pragma once

#include "dxmc/material.h"
#include <array>
#include <vector>

/**
 * @brief Attenuation look up tables.
 *
 * Class for generation and fast access of total and partial attenuation
 * coefficients for materials specified. Valid energies are from 1 keV up to
 * about 1000 keV (before pair production in electron fields may occur). This
 * class is typical used internally in the World class.
 */
namespace dxmc {

class AttenuationLut {
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
    double energyResolution() const { return m_energyStep; }

    /**
   * @brief Generate lookup tables for attenuation data
   * The order of materials matters, i.e the first material will get index 0,
   * second has index 1 etc...
   * @param materials Vector of materials to generate attenuation table
   * @param minEnergy Minimum energy to consider in keV
   * @param maxEnergy Maximum energy to consider in keV
   */
    void generate(const std::vector<Material>& materials, double minEnergy = 0.0,
        double maxEnergy = 150.0);

    /**
   * @brief Return the maximum mass attenuation value for all materials at
   * spesific photon energy
   *
   * @param energy Photon energy for maximum mass attenuation in table
   * @return double
   */
    double maxMassTotalAttenuation(double energy) const;

    /**
   * @brief Calculate max total mass attenuation for a material and density
   * array This templated function will iterate over material indices and
   * density indices and calculate maximum total mass attenuation value for each
   * photon energy in the lookup table
   *
   * @param materialIndexBegin Pointer or iterator pointing to the first index
   * in the material array
   * @param materialIndexEnd Pointer or iterator pointing to the last + 1 index
   * in the material array
   * @param densityBegin Pointer or iterator pointing to the first index in the
   * density array
   */
    template <typename It1, typename It2>
    void generateMaxMassTotalAttenuation(It1 materialIndexBegin,
        It1 materialIndexEnd, It2 densityBegin);

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
   * @brief Photoelectric mass attenuation for a material at specified photon
   * energy
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
   * @brief Shortcut to returning photoelectric, compton and rayleigh,
   * respectivly, mass attenuation for a material at specified photon energy
   *
   * @param material material index
   * @param energy photon energy in keV
   * @return std::array<double, 3>
   */
    std::array<double, 3> photoComptRayAttenuation(std::size_t material,
        double energy) const;

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
    std::vector<double>::iterator energyEnd()
    {
        return m_attData.begin() + m_energyResolution;
    }

    /**
   * @brief Iterator to the first element of total attenuation data for
   * specified material
   *
   * @param material Material index
   * @return std::vector<double>::iterator
   */
    std::vector<double>::iterator attenuationTotalBegin(std::size_t material)
    {
        return m_attData.begin() + m_energyResolution + m_energyResolution * 4 * material;
    } // todo error checking

    /**
   * @brief Iterator to the end of total attenuation data for specified material
   *
   * @param material Material index
   * @return std::vector<double>::iterator
   */
    std::vector<double>::iterator attenuationTotalEnd(std::size_t material)
    {
        return m_attData.begin() + m_energyResolution + m_energyResolution * 4 * material + m_energyResolution;
    } // todo error checking

    /**
   * @brief Calculate squared momentumtransfer for a material at a cummulative
   * squared atomic form factor
   *
   * @param material Material index
   * @param cumFormFactorSquared cummulative atomic form factor
   * @return double
   */
    double momentumTransfer(std::size_t material,
        double cumFormFactorSquared) const;

    /**
   * @brief Calculate cummulative atomic form factor squared from the squared
   * momentum transfer for a specified material
   *
   * @param material Material index
   * @param momentumTransferSquared momentum transfer
   * @return double
   */
    double cumFormFactorSquared(std::size_t material,
        double momentumTransfer) const;

    /**
   * @brief Calculate compton scatter function (low energy approximation) for
   * current material and photon energy (momentum transfer)
   *
   * @param material Material index
   * @param momentumTransfer momentum transfer
   * @return double
   */
    double comptonScatterFactor(std::size_t material,
        double momentumTransfer) const;

    /**
   * @brief Momentum transfer for photon
   *
   * @param energy Photon energy in keV
   * @param angle Scattering angle in radians
   * @return double
   */
    static double momentumTransfer(double energy, double angle);

    /**
   * @brief Momentum transfer for photon
   *
   * @param energy Photon energy in keV
   * @param angle Cosine of scattering angle in radians
   * @return double
   */
    static double momentumTransferFromCos(double energy, double cosAngle);

    /**
   * @brief Max possible momentum transfer for photon
   *
   * @param energy Photon energy in keV
   * @return double
   */
    static double momentumTransferMax(double energy);

    /**
   * @brief Cosine of scatter angle for a photon
   *
   * @param energy  Photon energy in keV
   * @param momentumTransfer momentumtransfer in $Å^{-1}$
   * @return double
   */
    static double cosAngle(double energy, double momentumTransfer);

protected:
    /**
   * @brief Generate atomic form factors for materials specified.
   *
   * @param materials Vector of Materials
   */
    void generateFFdata(const std::vector<Material>& materials);

    /**
   * @brief Generate atomic scatter factors for materials specified.
   *
   * @param materials Vector of Materials
   */
    void generateSFdata(const std::vector<Material>& materials);

    /**
   * @brief Template function for attenuation data look-up.
   *
   * @tparam N Type of interaction N= 0:total attenuation, 1:photoelectric
   * effect, 2:compton interaction, 3:rayleight scatter
   * @param materialIdx Material index
   * @param energy Photon energy in keV
   * @return double Mass-attenuation coefficient in cm2/g
   */
    template <std::size_t N>
    double attenuationDataFromTypeNumber(std::size_t materialIdx,
        double energy) const;

private:
    double m_minEnergy = 0;
    double m_maxEnergy = 150.0;
    double m_energyStep = 1.0;
    std::size_t m_energyResolution = 150;
    std::size_t m_momtranfResolution = 128;
    double m_momtranfStep = 0;

    std::size_t m_nMaterials = 0;
    std::vector<double> m_attData; // energy, array-> total, photo, compton, rauleight
    std::vector<double> m_rayleighFormFactorSqr; // qsquared, array-> A(qsquared)
    std::vector<double> m_comptonScatterFactor;
    std::vector<double> m_maxMassAtt;
};

/**
 * @brief Compute maximum mass-attenuation for this world
 * @tparam It1 Iterator type
 * @tparam It2 Iterator type
 * @param materialIndexBegin Material index start iterator 
 * @param materialIndexEnd Material index end iterator
 * @param densityBegin Density array start iterator
*/
template <typename It1, typename It2>
void AttenuationLut::generateMaxMassTotalAttenuation(It1 materialIndexBegin,
    It1 materialIndexEnd,
    It2 densityBegin)
{
    std::vector<double> maxDens(m_nMaterials, 0.0);

    while (materialIndexBegin != materialIndexEnd) {
        maxDens[*materialIndexBegin] = std::max(maxDens[*materialIndexBegin], *densityBegin);
        ++materialIndexBegin;
        ++densityBegin;
    }

    m_maxMassAtt.resize(m_energyResolution);
    std::fill(m_maxMassAtt.begin(), m_maxMassAtt.end(), 0.0);
    for (std::size_t material = 0; material < m_nMaterials; ++material) {
        for (std::size_t i = 0; i < m_energyResolution; ++i) {
            m_maxMassAtt[i] = std::max(m_maxMassAtt[i], totalAttenuation(material, m_attData[i]) * maxDens[material]);
        }
    }
}
}