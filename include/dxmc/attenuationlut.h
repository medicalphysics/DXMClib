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

#include "dxmc/constants.h"
#include "dxmc/dxmcrandom.h"
#include "dxmc/floating.h"
#include "dxmc/interpolation.h"
#include "dxmc/material.h"
#include "dxmc/world.h"

#include <array>
#include <vector>
namespace dxmc {

/**
 * @brief Attenuation look up tables.
 *
 * Class for generation and fast access of total and partial attenuation
 * coefficients for materials specified. Valid energies are from 1 keV up to
 * about 1000 keV (before pair production in electron fields may occur). This
 * class is typical used internally in the World class.
 */

template <Floating T = double>
class AttenuationLut {
public:
    static constexpr T MAX_PHOTON_ENERGY() { return 2 * ELECTRON_REST_MASS<T>(); }
    /**
   * @brief Construct a new Attenuation Lut object
   *
   */
    AttenuationLut() {};
    AttenuationLut(const World<T>& world, T minEnergy = 1.0, T maxEnergy = 150, T energyStep = 1.0)
    {
        generate(world, minEnergy, maxEnergy, energyStep);
    }

    /**
   * @brief Energy resolution of the lookup table
   *
   * Sets the energy resolution of the lookup table, minimum value is 0.1 keV
   * @param keV Energy resolution in keV, default value is 1.0 keV
   */
    void setEnergyResolution(T keV) { m_energyStep = keV > T { 0.1 } ? keV : T { 0.1 }; }

    /**
   * @brief Energy resolution in keV
   *
   * @return double
   */
    T energyResolution() const { return m_energyStep; }

    /**
   * @brief Generate lookup tables for attenuation data
   * Generate lookup tables from a valid world
   * @param world valid World object to generate attenuation table from
   * @param minEnergy Minimum energy to consider in keV, minimum 1.0 keV
   * @param maxEnergy Maximum energy to consider in keV, minimum minEnergy+energyStep
   * @param energyStep Energy stepsize in lut table, minimum 0.1 keV
   */
    void generate(const World<T>& world, T minEnergy = 1.0, T maxEnergy = 150, T energyStep = 1.0)
    {
        const auto& materials = world.materialMap();
        generate(materials, minEnergy, maxEnergy, energyStep);

        auto mat_beg = world.materialIndexArray()->cbegin();
        auto mat_end = world.materialIndexArray()->cend();
        auto dens_beg = world.densityArray()->cbegin();
        generateMaxMassTotalAttenuation(mat_beg, mat_end, dens_beg);
    }

    /**
   * @brief Generate lookup tables for attenuation data
   * The order of materials matters, i.e the first material will get index 0,
   * second has index 1 etc...
   * @param materials Vector of materials to generate attenuation table
   * @param minEnergy Minimum energy to consider in keV, minimum 1.0 keV
   * @param maxEnergy Maximum energy to consider in keV, minimum minEnergy+energyStep, maximum 2mc^2 keV
   * @param energyStep Energy stepsize in lut table, minimum 0.1 keV
   */
    void generate(const std::vector<Material>& materials, T minEnergy = 1, T maxEnergy = 150, T energyStep = 1.0)
    {
        m_energyStep = std::max(T { 0.1 }, energyStep);
        m_minEnergy = std::max(T { 0 }, std::min(maxEnergy, minEnergy));
        m_maxEnergy = std::min(MAX_PHOTON_ENERGY(), std::max(maxEnergy, minEnergy));
        m_maxEnergy = std::max(m_minEnergy + m_energyStep, m_maxEnergy);

        m_energyResolution = 1 + static_cast<std::size_t>(std::ceil((m_maxEnergy - m_minEnergy) / m_energyStep));
        m_maxEnergy = m_minEnergy + (m_energyResolution - 1) * m_energyStep;

        m_nMaterials = materials.size();

        m_attData.resize(m_energyResolution * m_nMaterials * 4 + m_energyResolution);
        m_meanBindingEnergy.resize(m_nMaterials);
        // generating energy table
        for (std::size_t i = 0; i < m_energyResolution; ++i)
            m_attData[i] = m_minEnergy + i * m_energyStep;

        for (std::size_t m = 0; m < m_nMaterials; ++m) {
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

            //getting mean binding energy
            m_meanBindingEnergy[m] = static_cast<T>(materials[m].getMeanBindingEnergy());
        }

        generateFFdata(materials);
        generateSFdata(materials);

        //getting electron configurations
        m_electronShellConfiguration.reserve(materials.size());
        for (const auto m : materials) {
            m_electronShellConfiguration.push_back(m.getElectronConfiguration<T>());
        }

        // maxenergies
        std::vector<T> densArray;
        std::vector<std::uint8_t> matArray;
        for (std::size_t i = 0; i < m_nMaterials; ++i) {
            const auto dens = static_cast<T>(materials[i].standardDensity());
            densArray.push_back(dens);
            matArray.push_back(static_cast<std::uint8_t>(i));
        }
        generateMaxMassTotalAttenuation(matArray.begin(), matArray.end(),
            densArray.begin());
    }

    /**
   * @brief Return the maximum mass attenuation value for all materials at
   * spesific photon energy
   *
   * @param energy Photon energy for maximum mass attenuation in table
   * @return double
   */
    T maxMassTotalAttenuation(T energy) const
    {
        if (energy <= m_minEnergy)
            return m_maxMassAtt[0];

        const std::size_t idx = static_cast<std::size_t>((energy - m_minEnergy) / m_energyStep);

        if (idx >= m_energyResolution - 1)
            return m_maxMassAtt[m_energyResolution - 1];
        return interp(m_attData[idx], m_attData[idx + 1], m_maxMassAtt[idx], m_maxMassAtt[idx + 1], energy);
    }

    T meanBindingEnergy(std::size_t materialIdx) const
    {
        return m_meanBindingEnergy[materialIdx];
    }

    /**
   * @brief Total mass attenuation for a material at specified photon energy
   *
   * @param material material index
   * @param energy photon energy in keV
   * @return double
   */
    T totalAttenuation(std::size_t material, T energy) const
    {
        return attenuationDataFromTypeNumber<0>(material, energy);
    }

    /**
   * @brief Rayleigh mass attenuation for a material at specified photon energy
   *
   * @param material material index
   * @param energy photon energy in keV
   * @return double
   */
    T rayleightAttenuation(std::size_t material, T energy) const
    {
        return attenuationDataFromTypeNumber<3>(material, energy);
    }

    /**
   * @brief Photoelectric mass attenuation for a material at specified photon
   * energy
   *
   * @param material material index
   * @param energy photon energy in keV
   * @return double
   */
    T photoelectricAttenuation(std::size_t material, T energy) const
    {
        return attenuationDataFromTypeNumber<1>(material, energy);
    }

    /**
   * @brief Compton mass attenuation for a material at specified photon energy
   *
   * @param material material index
   * @param energy photon energy in keV
   * @return double
   */
    T comptonAttenuation(std::size_t material, T energy) const
    {
        return attenuationDataFromTypeNumber<2>(material, energy);
    }

    /**
   * @brief Shortcut to returning photoelectric, compton and rayleigh,
   * respectivly, mass attenuation for a material at specified photon energy
   *
   * @param material material index
   * @param energy photon energy in keV
   * @return std::array<double, 3>
   */
    std::array<T, 3> photoComptRayAttenuation(std::size_t material,
        T energy) const
    {
        const std::size_t offset = m_energyResolution + m_energyResolution * material * 4 + m_energyResolution;
        const std::size_t idx = static_cast<std::size_t>((energy - m_minEnergy) / m_energyStep);

        std::array<T, 3> att;
        if (energy <= m_minEnergy) {
            for (std::size_t i = 0; i < 3; ++i) {
                att[i] = m_attData[offset + m_energyResolution * i];
            }
        } else if (idx >= m_energyResolution - 1) {
            for (std::size_t i = 0; i < 3; ++i) {
                att[i] = m_attData[offset + m_energyResolution * (i + 1) - 1];
            }
        } else {
            for (std::size_t i = 0; i < 3; ++i) {
                att[i] = interp(&m_attData[idx],
                    &m_attData[offset + m_energyResolution * i + idx], energy);
            }
        }
        return att;
    }

    /**
   * @brief Iterator to the first element of photon energy data
   *
   * @return std::vector<double>::iterator
   */
    typename std::vector<T>::iterator energyBegin() { return m_attData.begin(); }

    /**
   * @brief Iterator to the end of photon energy data
   *
   * @return std::vector<double>::iterator
   */
    typename std::vector<T>::iterator energyEnd()
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
    typename std::vector<T>::iterator attenuationTotalBegin(std::size_t material)
    {
        return m_attData.begin() + m_energyResolution + m_energyResolution * 4 * material;
    } // todo error checking

    /**
   * @brief Iterator to the end of total attenuation data for specified material
   *
   * @param material Material index
   * @return std::vector<double>::iterator
   */
    typename std::vector<T>::iterator attenuationTotalEnd(std::size_t material)
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
    T momentumTransfer(std::size_t material,
        T cumFormFactorSquared) const
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

        return interp(&m_rayleighFormFactorSqr[offset + distance - 1],
            &m_rayleighFormFactorSqr[distance - 1], cumFormFactorSquared);
    }

    /**
   * @brief Calculate cummulative atomic form factor squared from the squared
   * momentum transfer for a specified material
   *
   * @param material Material index
   * @param momentumTransferSquared momentum transfer
   * @return double
   */
    T cumFormFactorSquared(std::size_t material,
        T momentumTransfer) const
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
        return interp(&m_rayleighFormFactorSqr[distance - 1],
            &m_rayleighFormFactorSqr[offset + distance - 1],
            momentumTransfer);
    }

    /**
   * @brief Calculate compton scatter function (low energy approximation) for
   * current material and photon energy (momentum transfer)
   *
   * @param material Material index
   * @param momentumTransfer momentum transfer
   * @return double
   */
    T comptonScatterFactor(std::size_t material,
        T momentumTransfer) const
    {
        const std::size_t offset = m_momtranfResolution + m_momtranfResolution * material;
        if (momentumTransfer <= 0.0)
            return m_comptonScatterFactor[offset];

        const std::size_t idx = static_cast<std::size_t>(momentumTransfer / m_momtranfStep);
        if (idx >= m_momtranfResolution - 1)
            return m_comptonScatterFactor[offset + m_momtranfResolution - 1];
        return interp(&m_comptonScatterFactor[idx],
            &m_comptonScatterFactor[idx + offset], momentumTransfer);
    }

    /**
   * @brief Electron shell configuration for material
   *
   * @param materialIdx material index   
   * @return vector of ElectronShellConfiguration structs for each shell in material.
   */
    const std::array<ElectronShellConfiguration<T>, 12>& electronShellConfiguration(std::size_t materialIdx) const
    {
        return m_electronShellConfiguration[materialIdx];
    }

    /**
   * @brief Momentum transfer for photon
   *
   * @param energy Photon energy in keV
   * @param angle Scattering angle in radians
   * @return double
   */
    static T momentumTransfer(T energy, T angle)
    {
        constexpr T k = T { 1.0 } / KEV_TO_ANGSTROM<T>(); // constant for momentum transfer
        return energy * std::sin(angle / T { 2.0 }) * k;
    }

    /**
   * @brief Momentum transfer for photon
   *
   * @param energy Photon energy in keV
   * @param angle Cosine of scattering angle in radians
   * @return double
   */
    static T momentumTransferFromCos(T energy, T cosAngle)
    {
        constexpr T k = T { 1.0 } / KEV_TO_ANGSTROM<T>(); // constant for momentum transfer
        return energy * k * std::sqrt(T { 0.5 } - cosAngle * T { 0.5 });
    }

    /**
   * @brief Max possible momentum transfer for photon
   *
   * @param energy Photon energy in keV
   * @return double
   */
    static T momentumTransferMax(T energy)
    {
        constexpr T k = T { 1.0 } / KEV_TO_ANGSTROM<T>(); // constant for momentum transfer
        return energy * k;
    }

    /**
   * @brief Cosine of scatter angle for a photon
   *
   * @param energy  Photon energy in keV
   * @param momentumTransfer momentumtransfer in $Å^{-1}$
   * @return double
   */
    static T cosAngle(T energy, T momentumTransfer)
    {
        const auto invE = KEV_TO_ANGSTROM<T>() / energy;

        return T { 1.0 } - T { 2.0 } * momentumTransfer * momentumTransfer * invE * invE;
    }

protected:
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
        It1 materialIndexEnd, It2 densityBegin)
    {
        std::vector<T> maxDens(m_nMaterials, T { 0 });

        for (std::size_t i = 0; i < m_nMaterials; ++i) {
            maxDens[i] = std::transform_reduce(
                std::execution::par_unseq, materialIndexBegin, materialIndexEnd, densityBegin, T { 0 }, [](auto d1, auto d2) { return std::max(d1, d2); }, [=](auto m, auto d) -> T { return m == i ? d : 0; });
        }

        m_maxMassAtt.resize(m_energyResolution);
        std::fill(m_maxMassAtt.begin(), m_maxMassAtt.end(), T { 0.0 });
        for (std::size_t material = 0; material < m_nMaterials; ++material) {
            for (std::size_t i = 0; i < m_energyResolution; ++i) {
                m_maxMassAtt[i] = std::max(m_maxMassAtt[i], totalAttenuation(material, m_attData[i]) * maxDens[material]);
            }
        }
    }
    /**
   * @brief Generate atomic form factors for materials specified.
   *
   * @param materials Vector of Materials
   */
    void generateFFdata(const std::vector<Material>& materials)
    {
        // se http://rcwww.kek.jp/research/egs/egs5_manual/slac730-150228.pdf

        const auto globalMomMax = momentumTransferMax(m_maxEnergy);

        std::vector<T> qvalues(m_momtranfResolution);
        for (std::size_t i = 0; i < m_momtranfResolution; ++i)
            qvalues[i] = (i * globalMomMax / (m_momtranfResolution - 1)) * (i * globalMomMax / (m_momtranfResolution - 1)) / globalMomMax; // nonlinear integration steps

        m_rayleighFormFactorSqr.resize(m_momtranfResolution + materials.size() * m_momtranfResolution);
        std::copy(qvalues.cbegin(), qvalues.cend(), m_rayleighFormFactorSqr.begin());

        for (std::size_t i = 0; i < materials.size(); ++i) {
            auto ffmaxsqr = materials[i].getRayleightFormFactorSquared(qvalues);
            auto ffmaxsqrint = trapz(ffmaxsqr, qvalues);
            auto start = m_rayleighFormFactorSqr.begin() + m_momtranfResolution + i * m_momtranfResolution;
            std::copy(ffmaxsqrint.cbegin(), ffmaxsqrint.cend(), start);
        }
    }

    /**
   * @brief Generate atomic scatter factors for materials specified.
   *
   * @param materials Vector of Materials
   */
    void generateSFdata(const std::vector<Material>& materials)
    {
        // Compton scatter corrections

        // se http://rcwww.kek.jp/research/egs/egs5_manual/slac730-150228.pdf

        const auto globalMomMax = momentumTransferMax(m_maxEnergy);
        std::vector<T> qvalues(m_momtranfResolution);
        m_momtranfStep = globalMomMax / (m_momtranfResolution - 1);
        for (std::size_t i = 0; i < m_momtranfResolution; ++i)
            qvalues[i] = i * m_momtranfStep;

        m_comptonScatterFactor.resize(m_momtranfResolution + materials.size() * m_momtranfResolution);
        std::copy(qvalues.cbegin(), qvalues.cend(), m_comptonScatterFactor.begin());
        for (std::size_t i = 0; i < materials.size(); ++i) {
            auto sfmax = materials[i].getComptonNormalizedScatterFactor(qvalues);
            auto start = m_comptonScatterFactor.begin() + m_momtranfResolution + i * m_momtranfResolution;
            std::copy(sfmax.cbegin(), sfmax.cend(), start);
        }
    }

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
    T attenuationDataFromTypeNumber(std::size_t materialIdx, T energy) const
    {
        const std::size_t offset = m_energyResolution + m_energyResolution * materialIdx * 4 + m_energyResolution * N;
        if (energy <= m_minEnergy)
            return m_attData[offset];

        const std::size_t idx = static_cast<std::size_t>((energy - m_minEnergy) / m_energyStep);
        if (idx >= m_energyResolution - 1)
            return m_attData[m_energyResolution - 1 + offset];

        return interp(&m_attData[idx], &m_attData[idx + offset], energy);
    }

private:
    T m_minEnergy = 0;
    T m_maxEnergy = 150.0;
    T m_energyStep = 1.0;
    std::size_t m_energyResolution = 150;
    std::size_t m_momtranfResolution = 128;
    T m_momtranfStep = 0;

    std::size_t m_nMaterials = 0;
    std::vector<T> m_attData; // energy, array-> total, photo, compton, rauleight
    std::vector<T> m_rayleighFormFactorSqr; // qsquared, array-> A(qsquared)
    std::vector<T> m_comptonScatterFactor;
    std::vector<T> m_maxMassAtt;
    std::vector<T> m_meanBindingEnergy;
    std::vector<std::array<ElectronShellConfiguration<T>, 12>> m_electronShellConfiguration;
};
}