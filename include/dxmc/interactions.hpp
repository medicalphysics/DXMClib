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
#include "dxmc/floating.hpp"
#include "dxmc/material/material.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"

namespace dxmc {
namespace interactions {

    template <Floating T>
    constexpr static T russianRuletteProbability()
    {
        return T { 0.9 };
    }
    template <Floating T>
    constexpr static T russianRuletteWeightThreshold()
    {
        return T { 0.1 };
    }

    template <Floating T, int Nshells, int Lowenergycorrection = 3>
    void rayleightScatter(Particle<T>& particle, const Material2<T, Nshells>& material, RandomState& state) noexcept
    {
        if constexpr (Lowenergycorrection == 0) {
            bool reject;
            T theta;
            do {
                constexpr T extreme = (4 * std::numbers::sqrt2_v<T>) / (3 * std::numbers::sqrt3_v<T>);
                const auto r1 = state.randomUniform<T>(T { 0 }, extreme);
                theta = state.randomUniform(T { 0 }, PI_VAL<T>());
                const auto sinang = std::sin(theta);
                reject = r1 > ((2 - sinang * sinang) * sinang);
            } while (reject);
            // calc angle and add randomly 90 degrees since dist i symetrical
            const auto phi = state.randomUniform<T>(PI_VAL<T>() + PI_VAL<T>());
            vectormath::peturb<T>(particle.dir, theta, phi);
        } else {
            // theta is scattering angle
            // see http://rcwww.kek.jp/research/egs/egs5_manual/slac730-150228.pdf

            // finding qmax
            const auto qmax = material.momentumTransferMax(particle.energy);
            const auto qmax_squared = qmax * qmax;

            T cosAngle;
            do {
                const auto q_squared = material.sampleSquaredMomentumTransferFromFormFactorSquared(qmax_squared, state);
                cosAngle = T { 1 } - 2 * q_squared / qmax_squared;
            } while ((1 + cosAngle * cosAngle) * T { 0.5 } < state.randomUniform<T>());

            const auto theta = std::acos(cosAngle);
            const auto phi = state.randomUniform<T>(PI_VAL<T>() + PI_VAL<T>());
            vectormath::peturb<T>(particle.dir, theta, phi);
        }
    }

    template <Floating T, int Nshells>
    T comptonScatterIA(Particle<T>& particle, const Material2<T, Nshells>& material, RandomState& state) noexcept
    {
        // Penelope model for compton scattering. Note spelling error in manual for sampling of pz.
        // In addition we use hartree scatter factors instead of integrating all compton profiles when sampling cosTheta
        const auto k = particle.energy / ELECTRON_REST_MASS<T>();

        // sample cosangle and e
        T e, cosTheta;
        {
            const auto emin = 1 / (1 + 2 * k);
            bool rejected;
            const auto gmaxInv = emin / (1 + emin * emin);
            do {
                const auto r1 = state.randomUniform<T>();
                e = r1 + (1 - r1) * emin;
                const auto t = std::min((1 - e) / (k * e), T { 2 }); // to prevent rounding errors with t > 2 (better way?)
                cosTheta = 1 - t;
                const auto sinThetaSqr = 1 - cosTheta * cosTheta;
                const auto g = (1 / e + e - sinThetaSqr) * gmaxInv;
                const auto q = material.momentumTransferCosAngle(particle.energy, cosTheta);
                const auto scatterFactor = material.scatterFactor(q);
                rejected = state.randomUniform<T>(material.effectiveZ()) > (g * scatterFactor);
            } while (rejected);
        }
        const auto kc = k * e;
        const auto qc = std::sqrt(k * k + kc * kc - 2 * k * kc * cosTheta);
        const auto alpha = qc / k + (kc * kc - kc * k * cosTheta) / (k * qc * qc);
        // sample shell and pz
        const auto nia = [](const T pz, const T J0) -> T {
            constexpr auto d1 = T { 1 } / std::numbers::sqrt2_v<T>;
            constexpr auto d2 = std::numbers::sqrt2_v<T>;
            const T p1 = d1 + d2 * J0 * std::abs(pz);
            const auto p2 = T { 0.5 } * std::exp(std::max(T { 0.5 } - p1 * p1, T { 0 }));
            return pz > 0 ? 1 - p2 : p2;
        };
        std::uint_fast16_t shellIdx;
        T pz;
        {
            std::array<T, Nshells + 1> shell_nia;
            std::transform(std::execution::unseq, material.shells().cbegin(), material.shells().cend(), shell_nia.begin(), [&nia, cosTheta, k](const auto& shell) {
                const auto U = shell.bindingEnergy / ELECTRON_REST_MASS<T>();
                if (U > k)
                    return T { 0 };
                const auto pi_max = (k * (k - U) * (1 - cosTheta) - U) / std::sqrt(2 * k * (k - U) * (1 - cosTheta) + U * U);
                const auto ni = nia(pi_max, shell.HartreeFockOrbital_0);
                return ni;
            });
            std::array<T, Nshells + 1> shell_probs;
            std::transform(shell_nia.cbegin(), shell_nia.cend(), material.shells().cbegin(), shell_probs.begin(), [](const auto n, const auto& shell) {
                return n * shell.numberOfElectrons;
            });
            const auto shell_probs_sum = std::reduce(std::execution::unseq, shell_probs.cbegin(), shell_probs.cend(), T { 0 });
            bool rejected;
            do {
                do {
                    const auto shell_probs_r = state.randomUniform(shell_probs_sum);
                    shellIdx = 0;
                    T accum = shell_probs[shellIdx];
                    while (shell_probs_r > accum) {
                        ++shellIdx;
                        accum += shell_probs[shellIdx];
                    }

                    const auto nia_shell = shell_nia[shellIdx];
                    const auto A = state.randomUniform(nia_shell);
                    constexpr auto d1 = T { 1 } / std::numbers::sqrt2_v<T>;
                    constexpr auto d2 = std::numbers::sqrt2_v<T>;
                    const auto J0 = material.shell(shellIdx).HartreeFockOrbital_0;
                    if (A < T { 0.5 }) {
                        pz = (d1 - std::sqrt(T { 0.5 } - std::log(2 * A))) / (d2 * J0);
                    } else {
                        pz = (std::sqrt(T { 0.5 } - std::log(2 * (1 - A))) - d1) / (d2 * J0);
                    }
                } while (pz < -1);

                const auto F_eval = [k, kc, qc, cosTheta](const auto p) {
                    T pn;
                    if (p < T { -0.2 })
                        pn = T { -0.2 };
                    else if (p < T { 0.2 })
                        pn = p;
                    else
                        pn = T { 0.2 };
                    const auto f = 1 + qc * pn * (1 + kc * (kc - k * cosTheta) / (qc * qc)) / k;
                    return std::max(f, T { 0 });
                };
                const auto Fmax = std::max(F_eval(T { 0.2 }), F_eval(T { -0.2 }));
                const auto F = F_eval(pz);
                rejected = state.randomUniform(Fmax) > F;
            } while (rejected);
        }
        const auto t = pz * pz;
        const auto eb1 = 1 - t * e * e;
        const auto eb2 = 1 - t * e * cosTheta;
        const auto sign = pz > 0 ? 1 : -1;
        const auto eb = e / eb1 * (eb2 + sign * std::sqrt(eb2 * eb2 - eb1 * (1 - t)));
        const auto E = particle.energy;
        particle.energy *= eb;
        const auto theta = std::acos(cosTheta);
        const auto phi = state.randomUniform(PI_VAL<T>() * 2);
        vectormath::peturb(particle.dir, theta, phi);
        return (E - particle.energy) * particle.weight;
    }
    template <Floating T, int Nshells, int Lowenergycorrection = 3>
    T comptonScatter(Particle<T>& particle, const Material2<T, Nshells>& material, RandomState& state) noexcept
    // see http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/PhysicsReferenceManual/fo/PhysicsReferenceManual.pdf
    // and
    // https://nrc-cnrc.github.io/EGSnrc/doc/pirs701-egsnrc.pdf
    {
        if constexpr (Lowenergycorrection == 2) {
            return comptonScatterIA(particle, material, state);
        } else {
            const auto k = particle.energy / ELECTRON_REST_MASS<T>();
            const auto emin = 1 / (1 + 2 * k);
            const auto gmaxInv = emin / (1 + emin * emin);

            T e, cosTheta;
            bool rejected;
            do {
                const auto r1 = state.randomUniform<T>();
                e = r1 + (1 - r1) * emin;
                const auto t = std::min((1 - e) / (k * e), T { 2 }); // to prevent rounding errors with t > 2 (better way?)
                cosTheta = 1 - t;
                const auto sinThetaSqr = 1 - cosTheta * cosTheta;
                const auto g = (1 / e + e - sinThetaSqr) * gmaxInv;
                if constexpr (Lowenergycorrection == 0) {
                    rejected = state.randomUniform<T>() > g;
                } else { // Livermore
                    const auto q = material.momentumTransferCosAngle(particle.energy, cosTheta);
                    const auto scatterFactor = material.scatterFactor(q);
                    rejected = state.randomUniform<T>(material.effectiveZ()) > (g * scatterFactor);
                }
            } while (rejected);

            const auto theta = std::acos(cosTheta);
            const auto phi = state.randomUniform<T>(PI_VAL<T>() + PI_VAL<T>());
            vectormath::peturb<T>(particle.dir, theta, phi);

            const auto E = particle.energy;
            particle.energy *= e;
            return (E - particle.energy) * particle.weight;
        }
    }

    template <Floating T, int Nshells>
    T photoelectricEffectIA(const T totalPhotoCrossSection, Particle<T>& particle, const Material2<T, Nshells>& material, RandomState& state) noexcept
    {
        // finding shell based on photoelectric cross section
        const std::uint_fast8_t max_shell = material.numberOfShells();
        std::uint_fast8_t shell = 0;
        T prob = state.randomUniform<T>();
        bool next;
        do {
            const auto& sh = material.shell(shell);
            if (sh.bindingEnergy < MIN_ENERGY<T>()) {
                shell = max_shell;
            }
            next = shell != max_shell;
            if (next && sh.bindingEnergy < particle.energy) {
                const auto shellCS = material.attenuationPhotoelectricShell(shell, particle.energy);
                const auto shellProb = shellCS / totalPhotoCrossSection;
                prob -= shellProb;
                next = prob > T { 0 };
            }
            if (next)
                ++shell;
        } while (next);

        T E = particle.energy * particle.weight;
        particle.energy = 0;
        if (shell != max_shell) {
            const auto& s = material.shell(shell);
            if (s.energyOfPhotonsPerInitVacancy > MIN_ENERGY<T>()) {
                particle.energy = s.energyOfPhotonsPerInitVacancy;
                E -= particle.energy * particle.weight;
                particle.weight *= s.numberOfPhotonsPerInitVacancy;
                const auto theta = state.randomUniform(PI_VAL<T>());
                const auto phi = state.randomUniform(PI_VAL<T>() + PI_VAL<T>());
                vectormath::peturb(particle.dir, theta, phi);
            }
        }
        return E;
    }
    template <Floating T, int Nshells, int Lowenergycorrection = 3>
    T photoelectricEffect(const T totalPhotoCrossSection, Particle<T>& particle, const Material2<T, Nshells>& material, RandomState& state) noexcept
    {
        if constexpr (Lowenergycorrection == 3) {
            return photoelectricEffectIA(totalPhotoCrossSection, particle, material, state);
        } else {
            const T E = particle.energy * particle.weight;
            particle.energy = 0;
            return E;
        }
    }

    template<Floating T>
    struct InteractionResult {
        T energyImparted;
        bool particleAlive;
        bool particleEnergyChanged;
    };


    template<Floating T, int Nshells, int Lowenergycorrection=3>
    InteractionResult<T> interact(const AttenuationValues<T>& attenuation, Particle<T>& particle, const Material2<T, Nshells>& material, RandomState& state)
    {        
        const auto r2 = state.randomUniform<T>(attenuation.sum());
        InteractionResult<T> res;
        if (r2 < attenuation.photoelectric) {
            const auto Ei = interactions::photoelectricEffect(attenuation.photoelectric, particle, material, state);
            res.energyImparted=Ei; 
            res.particleEnergyChanged = true;
        } else if (r2 < (attenuation.photoelectric + attenuation.coherent)) {
            const auto Ei = interactions::comptonScatter(particle, material, state);
            res.energyImparted = Ei;
            res.particleEnergyChanged = true;
        } else {
            interactions::rayleightScatter(particle, material, state);
        }
        if (p.energy < MIN_ENERGY<T>()) {
            res.particleAlive = false;
        } else {
            if (p.weight < interactions::russianRuletteWeightThreshold<T>()) {
                if (state.randomUniform<T>() < interactions::russianRuletteProbability<T>()) {
                    res.particleAlive = false;
                } else {
                    p.weight /= (T { 1 } - interactions::russianRuletteProbability<T>());
                }
            }
        }
        return res;
    }

    


}
}
