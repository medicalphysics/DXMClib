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

    template <Floating T, int Lowenergycorrection, int Nshells>
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

    template <Floating T, int Lowenergycorrection, int Nshells>
    T comptonScatter(Particle<T>& particle, const Material2<T, Nshells>& material, RandomState& state) noexcept
    // see http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/PhysicsReferenceManual/fo/PhysicsReferenceManual.pdf
    // and
    // https://nrc-cnrc.github.io/EGSnrc/doc/pirs701-egsnrc.pdf
    {
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
            } else { // Livermore and IA process
                const auto q = material.momentumTransferCosAngle(particle.energy, cosTheta);
                const auto scatterFactor = material.scatterFactor(q);
                // normalize scatterfactor
                rejected = state.randomUniform<T>(material.effectiveZ()) > (g * scatterFactor);
            }
        } while (rejected);

        if constexpr (Lowenergycorrection == 2) { // only IA process

            const auto kc = k * e;
            const auto cqc = std::sqrt(k * k + kc * kc - 2 * k * kc * cosTheta);
            std::array<T, Nshells + 1> nia_pz_max;
            std::transform(std::execution::unseq, material.shells().cbegin(), material.shells().cend(), nia_pz_max.begin(), [cosTheta, k](const auto& shell) -> T {
                const auto ku = shell.bindingEnergy / ELECTRON_REST_MASS<T>();
                const auto pz_part1 = k * (k - ku) * (1 - cosTheta);
                const auto pz_imax = (pz_part1 - ku) / std::sqrt(2 * pz_part1 + ku * ku);

                constexpr auto d2 = std::numbers::sqrt2_v<T>;
                constexpr auto d1 = T { 1 } / std::numbers::sqrt2_v<T>;

                const auto exp_part1 = d1 + d2 * shell.HartreeFockOrbital_0 * std::abs(pz_imax);
                const auto exp = std::exp(T { 0.5 } - exp_part1 * exp_part1);
                const auto ni_pz = pz_imax > T { 0 } ? 1 - T { 0.5 } * exp : T { 0.5 } * exp;
                return ni_pz;
            });

            // selecting shell for interaction
            std::array<T, Nshells + 1> shellProbs;
            std::transform(std::execution::unseq, material.shells().cbegin(), material.shells().cend(), nia_pz_max.cbegin(), shellProbs.begin(), [&particle](const auto& shell, const auto nia) -> T {
                return shell.bindingEnergy > particle.energy ? shell.numberOfElectronsFraction * nia : T { 0 };
            });
            // cummulative sum
            std::partial_sum(shellProbs.cbegin(), shellProbs.cend(), shellProbs.begin());
            T pz, F;
            do {
                do {
                    // select shell
                    const auto r3 = state.randomUniform<T>(shellProbs.back());
                    auto idx = std::lower_bound(shellProbs.cbegin(), shellProbs.cend(), r3);
                    const auto i = std::distance(shellProbs.cbegin(), idx);
                    const auto A = state.randomUniform<T>() * nia_pz_max[i];

                    if (A < T { 0.5 }) {
                        constexpr auto d2 = std::numbers::sqrt2_v<T>;
                        constexpr auto d1 = T { 1 } / std::numbers::sqrt2_v<T>;
                        constexpr auto ln2 = std::numbers::ln2_v<T>;
                        pz = (d1 - std::sqrt(T { 0.5 } - ln2 * A)) / (d2 * material.shell(i).HartreeFockOrbital_0);
                    } else {
                        constexpr auto d2 = std::numbers::sqrt2_v<T>;
                        constexpr auto d1 = T { 1 } / std::numbers::sqrt2_v<T>;
                        constexpr auto ln2 = std::numbers::ln2_v<T>;
                        pz = (std::sqrt(T { 0.5 } - ln2 * (1 - A)) - d1) / (d2 * material.shell(i).HartreeFockOrbital_0);
                    }
                } while (pz < -1);

                F = 1 + cqc * (1 + (kc * (kc - k * cosTheta)) / (cqc * cqc)) * pz / k;

            } while (state.randomUniform<T>(T { 0.2 }) < F);

            // correcting value of e
            const auto t = pz * pz;
            const auto epart1 = 1 - t * e * cosTheta;
            const auto epart2 = 1 - t * e * e;
            const auto sign = (pz > 0) - (pz < 0);
            e = e * (epart1 + sign * std::sqrt(epart1 * epart1 - epart2 * (1 - t))) / epart2;
        }

        const auto theta = std::acos(cosTheta);
        const auto phi = state.randomUniform<T>(PI_VAL<T>() + PI_VAL<T>());
        vectormath::peturb<T>(particle.dir, theta, phi);

        const auto E = particle.energy;
        particle.energy *= e;
        return E - particle.energy;

        /*
        if constexpr (Lowenergycorrection == 2) {
            return comptonScatterNRC(particle, materialIdx, state);
        } else {

            const auto E = particle.energy;
            const auto k = E / ELECTRON_REST_MASS<T>();
            const auto emin = 1 / (1 + 2 * k);
            const auto gmax_inv = 1 / (1 / emin + emin);
            T e, cosAngle, Ee;
            bool rejected;
            do {
                const auto r1 = state.randomUniform<T>();
                e = r1 + (1 - r1) * emin;
                const auto t = (1 - e) / (k * e);
                const auto sinthetasqr = t * (2 - t);
                cosAngle = 1 - t;
                const auto g = (1 / e + e - sinthetasqr) * gmax_inv;
                const auto r2 = state.randomUniform<T>();
                if constexpr (Lowenergycorrection == 1) {
                    // simple correction with scatterfactor to supress forward scattering due to binding energy
                    const auto q = m_attenuationLut.momentumTransferFromCos(E, cosAngle);
                    const auto scatterFactor = m_attenuationLut.comptonScatterFactor(materialIdx, q);
                    rejected = r2 > g * scatterFactor;
                } else {
                    rejected = r2 > g;
                }

            } while (rejected);
            const auto theta = std::acos(cosAngle);
            const auto phi = state.randomUniform<T>(PI_VAL<T>() + PI_VAL<T>());
            vectormath::peturb<T>(particle.dir, theta, phi);
            particle.energy *= e;
            return E - particle.energy;
        }*/
    }
    /*
    T comptonScatterNRC(Particle<T>& particle, std::uint8_t materialIdx, RandomState& state) const noexcept
    {

        const auto& eConfig = m_attenuationLut.electronShellConfiguration(materialIdx);
        std::array<T, 12> shellProbs;

        // calculating shell probs
        std::transform(std::execution::unseq, eConfig.cbegin(), eConfig.cend(), shellProbs.begin(), [&](const auto& e) -> T {
            return e.bindingEnergy < particle.energy && e.hartreeFockOrbital_0 > 0 ? e.numberElectrons : 0;
        });
        // cummulative sum of shell probs
        std::partial_sum(shellProbs.cbegin(), shellProbs.cend(), shellProbs.begin());
        // selectring shell
        int shellIdx = 0;
        const auto shellIdxSample = shellProbs.back() * state.randomUniform<T>();
        while (shellProbs[shellIdx] < shellIdxSample && shellIdx < 11) {
            ++shellIdx;
        }

        const auto U = eConfig[shellIdx].bindingEnergy / ELECTRON_REST_MASS<T>();
        const auto p = std::sqrt(2 * U + U * U);
        const auto J0 = eConfig[shellIdx].hartreeFockOrbital_0;

        const auto k = particle.energy / ELECTRON_REST_MASS<T>();

        // Rejecting interaction if binding energy is larger than particle energy
        // Rejecting the interaction will correct the compton cross section
        if (U > k) {
            return 0;
        }

        const auto emin = 1 / (1 + 2 * k);
        const auto gmax_inv = 1 / (1 / emin + emin);

        T e, cosAngle;
        bool rejected;
        do {
            const auto r1 = state.randomUniform<T>();
            e = r1 + (1 - r1) * emin;
            const auto t = (1 - e) / (k * e);
            const auto sinAngleSqr = t * (2 - t);
            cosAngle = 1 - t;
            const auto g = (1 / e + e - sinAngleSqr) * gmax_inv;
            const auto r2 = state.randomUniform<T>();
            rejected = r2 > g;

            if (!rejected) {
                const auto pi = (k * (k - U) * (1 - cosAngle) - U) / std::sqrt(2 * k * (k - U) * (1 - cosAngle) + U * U);
                const auto kc = k * e;
                const auto qc = std::sqrt(k * k + kc * kc - 2 * k * kc * cosAngle);
                const auto alpha = qc * (1 + kc * (kc - k * cosAngle) / (qc * qc)) / k;
                const auto b_part = 1 + 2 * J0 * std::abs(pi);
                const auto b = (b_part * b_part + 1) / 2;
                const auto expb = std::exp(-b);

                // calculating S
                T S;
                if (pi <= -p) {
                    S = (1 - alpha * p) * expb / 2;
                } else if (pi <= 0) {
                    const auto pabs = std::abs(p);
                    const auto piabs = std::abs(pi);
                    constexpr auto sp2 = 1 / (std::numbers::inv_sqrtpi_v<T> / std::numbers::sqrt2_v<T>);
                    const auto part1 = alpha * sp2 / (4 * J0);
                    constexpr auto a1 = T { 0.34802 };
                    constexpr auto a2 = T { -0.0958798 };
                    constexpr auto a3 = T { 0.7478556 };
                    const auto sqrte = std::sqrt(std::numbers::e_v<T>);
                    const auto tp = 1 / (1 + T { 0.332673 } * (1 + 2 * J0 * pabs));
                    const auto tpi = 1 / (1 + T { 0.332673 } * (1 + 2 * J0 * piabs));
                    const auto part2p = sqrte - expb * tp * (a1 + a2 * tp + a3 * tp * tp);
                    const auto part2pi = sqrte - expb * tpi * (a1 + a2 * tpi + a3 * tpi * tpi);
                    S = (1 - alpha * pi) * expb / 2 - part1 * (part2p - part2pi);
                } else if (pi < p) {
                    const auto pabs = std::abs(p);
                    const auto piabs = std::abs(pi);
                    constexpr auto sp2 = 1 / (std::numbers::inv_sqrtpi_v<T> / std::numbers::sqrt2_v<T>);
                    const auto part1 = alpha * sp2 / (4 * J0);
                    constexpr auto a1 = T { 0.34802 };
                    constexpr auto a2 = T { -0.0958798 };
                    constexpr auto a3 = T { 0.7478556 };
                    const auto sqrte = std::sqrt(std::numbers::e_v<T>);
                    const auto tp = 1 / (1 + T { 0.332673 } * (1 + 2 * J0 * pabs));
                    const auto tpi = 1 / (1 + T { 0.332673 } * (1 + 2 * J0 * piabs));
                    const auto part2p = sqrte - expb * tp * (a1 + a2 * tp + a3 * tp * tp);
                    const auto part2pi = sqrte - expb * tpi * (a1 + a2 * tpi + a3 * tpi * tpi);
                    S = 1 - (1 - alpha * pi) * expb / 2 - part1 * (part2p - part2pi);
                } else {
                    S = 1 - (1 - alpha * p) * expb / 2;
                }

                const auto r3 = state.randomUniform<T>();
                rejected = r3 > S;
                if (!rejected) {
                    // sampling pz
                    T Fmax;
                    if (pi <= -p) {
                        Fmax = 1 - alpha * p;
                    } else if (pi >= p) {
                        Fmax = 1 + alpha * p;
                    } else {
                        Fmax = 1 + alpha * pi;
                    }
                    const auto r4 = state.randomUniform<T>();
                    const auto r_bar2 = 2 * r4 * expb;
                    T pz;
                    if (r_bar2 < 1) {
                        const auto part = sqrt(1 - 2 * std::log(r_bar2));
                        pz = (1 - part) / (2 * J0) / ELECTRON_REST_MASS<T>();
                    } else {
                        const auto part = sqrt(1 - 2 * std::log(2 - r_bar2));
                        pz = (part - 1) / (2 * J0) / ELECTRON_REST_MASS<T>();
                    }

                    T Fpz;
                    if (pz <= -p) {
                        Fpz = 1 - alpha * p;
                    } else if (pz >= p) {
                        Fpz = 1 + alpha * p;
                    } else {
                        Fpz = 1 + alpha * pz;
                    }

                    const auto r5 = state.randomUniform<T>();
                    rejected = r5 > Fpz / Fmax;

                    if (!rejected) {
                        // calculate new energy
                        const auto part = std::sqrt(1 - 2 * e * cosAngle + e * e * (1 - pz * pz * sinAngleSqr));
                        const auto k_bar = kc / (1 - pz * pz * e * e) * (1 - pz * pz * e * cosAngle + pz * part);
                        e = k_bar / k;
                    }
                }
            }
        } while (rejected);
        const auto theta = std::acos(cosAngle);
        const auto phi = state.randomUniform<T>(PI_VAL<T>() + PI_VAL<T>());
        vectormath::peturb<T>(particle.dir, theta, phi);
        const auto E = particle.energy;
        particle.energy *= e;
        return E - particle.energy;
    }
    */
}
}
