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

#include "dxmc/attenuationinterpolator.h"
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
    static constexpr T MIN_PHOTON_ENERGY() { return T { 0.5 }; }
    /**
   * @brief Construct a new Attenuation Lut object
   *
   */
    AttenuationLut() { }
    AttenuationLut(const World<T>& world, T maxEnergy = 150, T minEnergy = 1)
    {
        generate(world, maxEnergy, minEnergy);
    }

    /**
   * @brief Generate lookup tables for attenuation data
   * Generate lookup tables from a valid world
   * @param world valid World object to generate attenuation table from
   * @param minEnergy Minimum energy to consider in keV, minimum 1.0 keV
   * @param maxEnergy Maximum energy to consider in keV, minimum minEnergy+energyStep
   */
    void generate(const World<T>& world, T maxEnergy = 150, T minEnergy = 1)
    {
        const auto& materials = world.materialMap();
        generate(materials, maxEnergy, minEnergy, false);
        m_attenuationData = AttenuationLutInterpolator<T>(world, m_maxEnergy, m_minEnergy);
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
    void generate(const std::vector<Material>& materials, T maxEnergy = 150, T minEnergy = 1, bool generatePhotonData = true)
    {
        m_minEnergy = std::max(MIN_PHOTON_ENERGY(), std::min(maxEnergy, minEnergy));
        m_maxEnergy = std::min(MAX_PHOTON_ENERGY(), std::max(maxEnergy, minEnergy));

        generateFFdata(materials);
        generateSFdata(materials);

        //getting electron configurations
        m_electronShellConfiguration.reserve(materials.size());
        for (const auto& m : materials) {
            m_electronShellConfiguration.push_back(m.getElectronConfiguration<T>());
        }
        if (generatePhotonData) {
            m_attenuationData = AttenuationLutInterpolator<T>(materials, m_maxEnergy, m_minEnergy);
        }
    }

    /**
   * @brief Return the maximum mass attenuation value for all materials at
   * spesific photon energy
   *
   * @param energy Photon energy for maximum mass attenuation in table
   * @return double
   */
    T maxTotalAttenuationInverse(T energy) const
    {
        return m_attenuationData.maxAttenuationInverse(energy);
    }

    /**
   * @brief Returning photoelectric, compton and rayleigh,
   * respectivly, mass attenuation for a material at specified photon energy
   *
   * @param material material index
   * @param energy photon energy in keV
   * @return std::array<double, 3>
   */
    std::array<T, 3> photoComptRayAttenuation(std::size_t material, T energy) const
    {

        return m_attenuationData(material, energy);
    }

    T momentumTransferFromFormFactor(std::size_t material, const T momentumTransferMax, RandomState& state) const
    {

        const T res = m_formFactor[material](state, momentumTransferMax);
        return res;
    }

    /**
   * @brief Calculate compton scatter function (low energy approximation) for
   * current material and photon energy (momentum transfer)
   *
   * @param material Material index
   * @param momentumTransfer momentum transfer
   * @return double
   */
    inline T comptonScatterFactor(std::size_t material, T momentumTransfer) const
    {
        return m_comptonScatterFactor[material](momentumTransfer);
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
        constexpr T k = 1 / KEV_TO_ANGSTROM<T>(); // constant for momentum transfer
        return energy * std::sin(angle * T { 0.5 }) * k;
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
        constexpr T k = 1 / KEV_TO_ANGSTROM<T>(); // constant for momentum transfer
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
        constexpr T k = 1 / KEV_TO_ANGSTROM<T>(); // constant for momentum transfer
        return energy * k;
    }

    /**
   * @brief Cosine of scatter angle for a photon
   *
   * @param energy  Photon energy in keV
   * @param momentumTransferSquared momentumtransfer squared in $Å^{-2}$
   * @return double
   */
    static inline T cosAngle(const T energy, const T momentumTransferSquared)
    {
        const auto invE = KEV_TO_ANGSTROM<T>() / energy;
        return 1 - 2 * momentumTransferSquared * invE * invE;
    }

protected:
    /**
   * @brief Generate atomic form factors for materials specified.
   *
   * @param materials Vector of Materials
   */
    void generateFFdata(const std::vector<Material>& materials)
    {
        m_formFactor.clear();
        m_formFactor.reserve(materials.size());
        const auto qmax = momentumTransferMax(m_maxEnergy);
        const auto qmax_squared = qmax * qmax;
        for (const auto& m : materials) {

            T qmax_applicable = 1;
            T ff_applicable = m.getRayleightFormFactorSquared(qmax_applicable);
            while (qmax_applicable < qmax_squared && ff_applicable > T { 0.001 }) {
                if (ff_applicable > T { 0.5 })
                    qmax_applicable += T { 0.5 };
                else
                    qmax_applicable += T { 0.1 };
                ff_applicable = m.getRayleightFormFactorSquared(qmax_applicable);
            }

            auto f = [&](auto q) { return m.getRayleightFormFactorSquared(std::sqrt(q)); };
            m_formFactor.emplace_back(T { 0 }, qmax_applicable, f);
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
        // find max applica
        m_comptonScatterFactor.clear();
        m_comptonScatterFactor.reserve(materials.size());
        for (const auto& m : materials) {
            // find best applicable max momentumtransfer
            const T qmax_energy = momentumTransferMax(m_maxEnergy);
            T qmax_applicable = 0.5;
            T sf_applicable = m.getComptonNormalizedScatterFactor(qmax_applicable);
            while (sf_applicable < T { 0.999 } && qmax_applicable < qmax_energy) {
                sf_applicable = m.getComptonNormalizedScatterFactor(qmax_applicable);
                if (sf_applicable < T { 0.5 })
                    qmax_applicable += T { 0.5 };
                else
                    qmax_applicable += T { 0.1 };
            }

            auto f = [&](const T q) -> T { return m.getComptonNormalizedScatterFactor(q); };

            m_comptonScatterFactor.emplace_back(T { 0 }, qmax_applicable, f);
        }
    }

private:
    T m_minEnergy = 0;
    T m_maxEnergy = 150.0;
    std::vector<CubicSplineInterpolator<T, 16>> m_comptonScatterFactor;
    std::vector<RITA<T, 56>> m_formFactor;
    AttenuationLutInterpolator<T> m_attenuationData;

    std::vector<std::array<ElectronShellConfiguration<T>, 12>> m_electronShellConfiguration;
};
}