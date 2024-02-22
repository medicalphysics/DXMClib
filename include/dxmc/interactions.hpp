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

Copyright 2023 Erlend Andersen
*/

#pragma once

#include "dxmc/constants.hpp"
#include "dxmc/dxmcrandom.hpp"
#include "dxmc/material/material.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"

namespace dxmc {
namespace interactions {

    constexpr static double russianRuletteProbability()
    {
        return 0.9;
    }

    constexpr static double russianRuletteWeightThreshold()
    {
        return 0.1;
    }

    template <std::size_t Nshells, int Lowenergycorrection = 2>
    void rayleightScatter(Particle& particle, const Material<Nshells>& material, RandomState& state) noexcept
    {
        if constexpr (Lowenergycorrection == 0) {
            bool reject;
            double theta;
            do {
                constexpr auto extreme = (4 * std::numbers::sqrt2_v<double>) / (3 * std::numbers::sqrt3_v<double>);
                const auto r1 = state.randomUniform<double>(0, extreme);
                theta = state.randomUniform<double>(0, PI_VAL<double>());
                const auto sinang = std::sin(theta);
                reject = r1 > ((2 - sinang * sinang) * sinang);
            } while (reject);
            // calc angle and add randomly 90 degrees since dist i symetrical
            const auto phi = state.randomUniform(PI_VAL() + PI_VAL());
            particle.dir = vectormath::peturb(particle.dir, theta, phi);

        } else {
            // theta is scattering angle
            // see http://rcwww.kek.jp/research/egs/egs5_manual/slac730-150228.pdf

            // finding qmax
            const auto qmax = material.momentumTransferMax(particle.energy);
            const auto qmax_squared = qmax * qmax;

            double cosAngle;
            do {
                const auto q_squared = material.sampleSquaredMomentumTransferFromFormFactorSquared(qmax_squared, state);
                cosAngle = 1 - 2 * q_squared / qmax_squared;
            } while ((1 + cosAngle * cosAngle) * 0.5 < state.randomUniform());

            const auto phi = state.randomUniform(PI_VAL() + PI_VAL());
            const auto theta = std::acos(cosAngle);
            particle.dir = vectormath::peturb(particle.dir, theta, phi);
        }
    }

    template <std::size_t Nshells>
    auto comptonScatterIA(Particle& particle, const Material<Nshells>& material, RandomState& state) noexcept
    {
        // Penelope model for compton scattering. Note spelling error in manual for sampling of pz.
        // In addition we use hartree scatter factors instead of integrating all compton profiles when sampling cosTheta
        const auto k = particle.energy / ELECTRON_REST_MASS();

        // sample cosangle and e
        double e, cosTheta;
        {
            const auto emin = 1 / (1 + 2 * k);
            bool rejected;
            const auto gmaxInv = emin / (1 + emin * emin);
            do {
                const auto r1 = state.randomUniform();
                e = r1 + (1 - r1) * emin;
                const auto t = std::min((1 - e) / (k * e), 2.0); // to prevent rounding errors with t > 2 (better way?)
                cosTheta = 1 - t;
                const auto sinThetaSqr = 1 - cosTheta * cosTheta;
                const auto g = (1 / e + e - sinThetaSqr) * gmaxInv;
                const auto q = material.momentumTransferCosAngle(particle.energy, cosTheta);
                const auto scatterFactor = material.scatterFactor(q);
                rejected = state.randomUniform(material.effectiveZ()) > (g * scatterFactor);
            } while (rejected);
        }
        const auto kc = k * e;
        const auto qc = std::sqrt(k * k + kc * kc - 2 * k * kc * cosTheta);
        // sample shell and pz
        const auto nia = [](const auto pz, const auto J0) -> double {
            constexpr auto d1 = 1 / std::numbers::sqrt2_v<double>;
            constexpr auto d2 = std::numbers::sqrt2_v<double>;
            const auto p1 = d1 + d2 * J0 * std::abs(pz);
            const auto p2 = 0.5 * std::exp(std::max(0.5 - p1 * p1, 0.0));
            return pz > 0 ? 1 - p2 : p2;
        };
        std::uint_fast16_t shellIdx;
        double pz;
        {
            std::array<double, Nshells + 1> shell_nia;
            std::transform(std::execution::unseq, material.shells().cbegin(), material.shells().cend(), shell_nia.begin(), [&nia, cosTheta, k](const auto& shell) -> double {
                const auto U = shell.bindingEnergy / ELECTRON_REST_MASS();
                if (U > k)
                    return 0.0;
                const auto pi_max = (k * (k - U) * (1 - cosTheta) - U) / std::sqrt(2 * k * (k - U) * (1 - cosTheta) + U * U);
                const auto ni = nia(pi_max, shell.HartreeFockOrbital_0);
                return ni;
            });
            std::array<double, Nshells + 1> shell_probs;
            std::transform(shell_nia.cbegin(), shell_nia.cend(), material.shells().cbegin(), shell_probs.begin(), [](const auto n, const auto& shell) -> double {
                return n * shell.numberOfElectrons;
            });
            const auto shell_probs_sum = std::reduce(std::execution::unseq, shell_probs.cbegin(), shell_probs.cend(), 0.0);
            bool rejected;
            do {
                do {
                    const auto shell_probs_r = state.randomUniform(shell_probs_sum);
                    shellIdx = 0;
                    auto accum = shell_probs[shellIdx];
                    while (shell_probs_r > accum) {
                        ++shellIdx;
                        accum += shell_probs[shellIdx];
                    }

                    const auto nia_shell = shell_nia[shellIdx];
                    const auto A = state.randomUniform(nia_shell);
                    constexpr auto d1 = 1 / std::numbers::sqrt2_v<double>;
                    constexpr auto d2 = std::numbers::sqrt2_v<double>;
                    const auto J0 = material.shell(shellIdx).HartreeFockOrbital_0;
                    if (A < 0.5) {
                        pz = (d1 - std::sqrt(0.5 - std::log(2 * A))) / (d2 * J0);
                    } else {
                        pz = (std::sqrt(0.5 - std::log(2 * (1 - A))) - d1) / (d2 * J0);
                    }
                } while (pz < -1);

                const auto F_eval = [k, kc, qc, cosTheta](const auto p) {
                    const auto pn = std::clamp(p, -0.2, 0.2);
                    const auto f = 1 + qc * pn * (1 + kc * (kc - k * cosTheta) / (qc * qc)) / k;
                    return std::max(f, 0.0);
                };
                const auto Fmax = std::max(F_eval(0.2), F_eval(-0.2));
                const auto F = F_eval(pz);
                rejected = state.randomUniform(Fmax) > F;
            } while (rejected);
        }
        const auto t = pz * pz;
        const auto eb1 = 1 - t * e * e;
        const auto eb2 = 1 - t * e * cosTheta;
        const auto sign = pz > 0 ? 1 : -1;

        const auto sqrt_factor = std::max(eb2 * eb2 - eb1 * (1 - t), 0.0); // may have small negative facors, fix

        const auto eb = e / eb1 * (eb2 + sign * std::sqrt(sqrt_factor));

        const auto E = particle.energy;
        particle.energy *= eb;
        const auto phi = state.randomUniform(PI_VAL() + PI_VAL());
        const auto theta = std::acos(cosTheta);
        particle.dir = vectormath::peturb(particle.dir, theta, phi);
        auto Ei = (E - particle.energy) * particle.weight;
        return Ei;
    }

    template <std::size_t Nshells, int Lowenergycorrection = 2>
    auto comptonScatter(Particle& particle, const Material<Nshells>& material, RandomState& state) noexcept
    // see http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/PhysicsReferenceManual/fo/PhysicsReferenceManual.pdf
    // and
    // https://nrc-cnrc.github.io/EGSnrc/doc/pirs701-egsnrc.pdf
    {
        if constexpr (Lowenergycorrection == 2) {
            return comptonScatterIA(particle, material, state);
        } else {
            const auto k = particle.energy / ELECTRON_REST_MASS();
            const auto emin = 1 / (1 + 2 * k);
            const auto gmaxInv = emin / (1 + emin * emin);

            double e, cosTheta;
            bool rejected;
            do {
                const auto r1 = state.randomUniform();
                e = r1 + (1 - r1) * emin;
                const auto t = std::min((1 - e) / (k * e), 2.0); // to prevent rounding errors with t > 2 (better way?)
                cosTheta = 1 - t;
                const auto sinThetaSqr = 1 - cosTheta * cosTheta;
                const auto g = (1 / e + e - sinThetaSqr) * gmaxInv;
                if constexpr (Lowenergycorrection == 0) {
                    rejected = state.randomUniform() > g;
                } else { // Livermore
                    const auto q = material.momentumTransferCosAngle(particle.energy, cosTheta);
                    const auto scatterFactor = material.scatterFactor(q);
                    rejected = state.randomUniform(material.effectiveZ()) > (g * scatterFactor);
                }
            } while (rejected);

            const auto phi = state.randomUniform(PI_VAL() + PI_VAL());
            const auto theta = std::acos(cosTheta);
            particle.dir = vectormath::peturb(particle.dir, theta, phi);

            const auto E = particle.energy;
            particle.energy *= e;
            return (E - particle.energy) * particle.weight;
        }
    }

    template <std::size_t Nshells>
    auto photoelectricEffectIA(const double totalPhotoCrossSection, Particle& particle, const Material<Nshells>& material, RandomState& state) noexcept
    {
        // finding shell based on photoelectric cross section
        const std::uint_fast8_t max_shell = material.numberOfShells();
        std::uint_fast8_t shell = 0;
        auto prob = state.randomUniform();
        bool next;
        do {
            const auto& sh = material.shell(shell);
            if (sh.bindingEnergy < MIN_ENERGY()) {
                shell = max_shell;
            }
            next = shell != max_shell;
            if (next && sh.bindingEnergy < particle.energy) {
                const auto shellCS = material.attenuationPhotoelectricShell(shell, particle.energy);
                const auto shellProb = shellCS / totalPhotoCrossSection;
                prob -= shellProb;
                next = prob > 0;
            }
            if (next)
                ++shell;
        } while (next);

        auto E = particle.energy * particle.weight;
        particle.energy = 0;
        if (shell != max_shell) {
            const auto& s = material.shell(shell);
            if (s.energyOfPhotonsPerInitVacancy > MIN_ENERGY()) {
                particle.energy = s.energyOfPhotonsPerInitVacancy;
                E -= particle.energy * particle.weight;
                const auto theta = state.randomUniform(PI_VAL());
                const auto phi = state.randomUniform(PI_VAL() + PI_VAL());
                particle.dir = vectormath::peturb(particle.dir, theta, phi);
            }
        }
        return E;
    }

    template <int Nshells, int Lowenergycorrection = 2>
    auto photoelectricEffect(const double totalPhotoCrossSection, Particle& particle, const Material<Nshells>& material, RandomState& state) noexcept
    {
        if constexpr (Lowenergycorrection == 2) {
            return photoelectricEffectIA(totalPhotoCrossSection, particle, material, state);
        } else {
            const auto E = particle.energy * particle.weight;
            particle.energy = 0;
            return E;
        }
    }

    struct InteractionResult {
        double energyImparted = 0;
        bool particleAlive = true;
        bool particleEnergyChanged = false;
        bool particleDirectionChanged = false;
    };

    template <std::size_t Nshells, int Lowenergycorrection = 2>
    InteractionResult interact(const AttenuationValues& attenuation, Particle& particle, const Material<Nshells>& material, RandomState& state)
    {
        InteractionResult res;
        const auto r2 = state.randomUniform(attenuation.sum());
        if (r2 < attenuation.photoelectric) {
            const auto Ei = interactions::photoelectricEffect<Nshells, Lowenergycorrection>(attenuation.photoelectric, particle, material, state);
            res.energyImparted = Ei;
            res.particleEnergyChanged = true;
        } else if (r2 < (attenuation.photoelectric + attenuation.incoherent)) {
            const auto Ei = interactions::comptonScatter<Nshells, Lowenergycorrection>(particle, material, state);
            res.energyImparted = Ei;
            res.particleEnergyChanged = true;
            res.particleDirectionChanged = true;
        } else {
            interactions::rayleightScatter<Nshells, Lowenergycorrection>(particle, material, state);
            res.particleDirectionChanged = true;
        }
        if (particle.energy < MIN_ENERGY()) {
            res.particleAlive = false;
            res.energyImparted += particle.energy * particle.weight;
        } else {
            if (particle.weight < interactions::russianRuletteWeightThreshold() && res.particleAlive) {
                if (state.randomUniform() < interactions::russianRuletteProbability()) {
                    res.particleAlive = false;
                    particle.energy = 0;
                } else {
                    constexpr auto factor = 1 / (1 - interactions::russianRuletteProbability());
                    particle.weight *= factor;
                    res.particleAlive = true;
                }
            }
        }
        return res;
    }

    template <std::size_t NMaterialShells, int LOWENERGYCORRECTION = 2>
    InteractionResult interactForced(double maxStepLen, double materialDensity, const AttenuationValues& attenuation, Particle& particle, const Material<NMaterialShells>& material, RandomState& state)
    {
        InteractionResult intRes;
        const auto relativePeProbability = attenuation.photoelectric / attenuation.sum();
        const auto attSum = attenuation.sum() * materialDensity;
        const auto probNotInteraction = std::exp(-attSum * maxStepLen);
        // Forced photoelectric effect
        intRes.energyImparted = particle.energy * particle.weight * (1 - probNotInteraction) * relativePeProbability;

        // Remainder probability
        const auto p1 = state.randomUniform();
        if (p1 > probNotInteraction) {
            if (state.randomUniform() > relativePeProbability) { // scatter interaction happens
                // Translate particle to interaction point
                const auto stepLen = -std::log(p1) / attSum;
                particle.translate(stepLen);

                // Decide what scatter interaction
                const auto r2 = state.randomUniform(attenuation.incoherent + attenuation.coherent);
                if (r2 < attenuation.incoherent) {
                    const auto Ei = interactions::comptonScatter<NMaterialShells, LOWENERGYCORRECTION>(particle, material, state);
                    intRes.energyImparted += Ei;
                    intRes.particleEnergyChanged = true;
                    intRes.particleDirectionChanged = true;
                } else {
                    interactions::rayleightScatter<NMaterialShells, LOWENERGYCORRECTION>(particle, material, state);
                    intRes.particleDirectionChanged = true;
                }

                // Handle cutoff values
                if (particle.energy < MIN_ENERGY()) {
                    intRes.particleAlive = false;
                    intRes.energyImparted += particle.energy;
                } else {
                    if (particle.weight < interactions::russianRuletteWeightThreshold() && intRes.particleAlive) {
                        if (state.randomUniform() < interactions::russianRuletteProbability()) {
                            intRes.particleAlive = false;
                        } else {
                            constexpr auto factor = 1 / (1 - interactions::russianRuletteProbability());
                            particle.weight *= factor;
                            intRes.particleAlive = true;
                        }
                    }
                }
            } else {
                // real photoelectric effect event
                // we don't score energy since it's already done. But terminates particle to prevent bias.
                intRes.particleAlive = false;
                particle.energy = 0;
            }
        } else {
            // No interaction event, transport particle to border
            particle.border_translate(maxStepLen);
        }
        return intRes;
    }

    template <std::size_t Nshells, int Lowenergycorrection = 2>
    InteractionResult interactForced(double maxStepLenght, double density, Particle& particle, const Material<Nshells>& material, RandomState& state)
    {
        const auto attenuation = material.attenuationValues(particle.energy);
        return interactForced<Nshells, Lowenergycorrection>(maxStepLenght, density, attenuation, particle, material, state);
    }

}
}
