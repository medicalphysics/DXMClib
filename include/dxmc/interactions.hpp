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

    template <Floating T, int Nshells>
    T comptonScatterPenelope(Particle<T>& particle, const Material2<T, Nshells>& material, RandomState& state) noexcept
    {
        const auto k = particle.energy / ELECTRON_REST_MASS<T>();

        const auto Smax = std::transform_reduce(std::execution::unseq, material.shells().cbegin(), material.shells().cend(), T { 0 }, std::plus<>(), [k](const auto& shell) -> T {
            const auto ku = shell.bindingEnergy / ELECTRON_REST_MASS<T>();
            if (ku > k)
                return T { 0 };
            const auto pi_max = (k * (k - ku) * 2 - ku) / std::sqrt(4 * k * (k - ku) + ku * ku);
            const auto J0 = shell.HartreeFockOrbital_0;
            const auto expb = std::exp(T { 0.5 } - (T { 0.5 } + 2 * J0 * std::abs(pi_max) + 2 * J0 * J0 * pi_max * pi_max));
            const auto nia = pi_max < 0 ? T { 0.5 } * expb : 1 - T { 0.5 } * expb;
            return shell.numberOfElectrons * nia;
        });

        T e, cosTheta, S;
        std::array<T, Nshells + 1> nia, S_arr;
        const auto emin = 1 / (1 + 2 * k);
        const auto gmaxInv = emin / (1 + emin * emin);
        bool rejected;
        do {
            const auto r1 = state.randomUniform<T>();
            e = r1 + (1 - r1) * emin;
            const auto t = std::min((1 - e) / (k * e), T { 2 }); // to prevent rounding errors with t > 2 (better way?)
            cosTheta = 1 - t;
            const auto sinThetaSqr = 1 - cosTheta * cosTheta;
            const auto g = (1 / e + e - sinThetaSqr) * gmaxInv;

            std::transform(std::execution::unseq, material.shells().cbegin(), material.shells().cend(), nia.begin(), [k, cosTheta](const auto& shell) {
                const auto ku = shell.bindingEnergy / ELECTRON_REST_MASS<T>();
                if (ku > k)
                    return T { 0 };
                const auto pi_max = (k * (k - ku) * (1 - cosTheta) - ku) / std::sqrt(2 * k * (k - ku) * (1 - cosTheta) + ku * ku);
                const auto J0 = shell.HartreeFockOrbital_0;
                const auto expb = std::exp(T { 0.5 } - (T { 0.5 } + 2 * J0 * std::abs(pi_max) + 2 * J0 * J0 * pi_max * pi_max));
                const auto nia = pi_max < 0 ? T { 0.5 } * expb : 1 - T { 0.5 } * expb;
                return nia;
            });
            std::transform(std::execution::unseq, material.shells().cbegin(), material.shells().cend(), nia.cbegin(), S_arr.begin(), [](const auto& shell, const auto n) {
                return shell.numberOfElectrons * n;
            });
            S = std::reduce(std::execution::unseq, S_arr.cbegin(), S_arr.cend(), T { 0 });

            rejected = state.randomUniform<T>() > g * S / Smax;

        } while (rejected);

        // picking shell
        T pz, F;
        const auto kc = k * e;
        const auto qc = std::sqrt(k * k + kc * kc - 2 * k * kc * cosTheta);
        constexpr auto p = T { 0.2 };
        const auto Fmax_part = 1 + (kc * (kc - k * cosTheta) / (qc * qc));

        const auto Fmax = Fmax_part < T { 0 } ? 1 - qc * Fmax_part * p / k : 1 + qc * Fmax_part * p / k;

        do {
            do {
                const auto Sr = state.randomUniform(S);
                std::uint_fast16_t shellIdx = 0;
                T shell_cum_prob = S_arr[shellIdx];
                while (shell_cum_prob < Sr) {
                    ++shellIdx;
                    shell_cum_prob += S_arr[shellIdx];
                }

                const auto& shell = material.shell(shellIdx);
                const auto A = state.randomUniform<T>() * nia[shellIdx];

                if (A < T { 0.5 }) {
                    pz = (1 / std::numbers::sqrt2_v<T> - std::sqrt(T { 0.5 } - std::numbers::ln2_v<T> * A)) / (std::numbers::sqrt2_v<T> * shell.HartreeFockOrbital_0);
                } else {
                    pz = (std::sqrt(T { 0.5 } - std::numbers::ln2_v<T> * (1 - A)) - 1 / std::numbers::sqrt2_v<T>) / (std::numbers::sqrt2_v<T> * shell.HartreeFockOrbital_0);
                }
            } while (pz < T { -1 });

            if (std::abs(pz) < p) {
                F = 1 + qc / k * (1 + (kc * (kc - k * cosTheta)) / (qc * qc)) * pz;
            } else {
                F = 1 + qc / k * (1 + (kc * (kc - k * cosTheta)) / (qc * qc)) * p;
            }
            F = std::max(F, T { 0 }); // may be negative due to approximation

        } while (state.randomUniform(Fmax) > F);
        const auto t = pz * pz;
        const auto e_part1 = 1 - t * e * e;
        const auto e_part2 = 1 - e * t * cosTheta;
        const auto pz_sign = pz > 0 ? 1 : -1;

        const auto e_b = e / e_part1 * (e_part2 + pz_sign * std::sqrt(e_part2 * e_part2 - e_part1 * (1 - t)));

        const auto E = particle.energy;
        particle.energy *= e_b;
        const auto theta = std::acos(cosTheta);
        const auto phi = state.randomUniform<T>(PI_VAL<T>() + PI_VAL<T>());
        vectormath::peturb<T>(particle.dir, theta, phi);

        return E - particle.energy;
    }

    template <Floating T, int Nshells>
    T comptonScatterIA(Particle<T>& particle, const Material2<T, Nshells>& material, RandomState& state) noexcept
    {
        const auto E = particle.energy;
        std::uint_fast32_t shellIdx = 0;

        const auto r1 = state.randomUniform<T>();
        T shell_cum_prob = material.shell(shellIdx).numberOfElectronsFraction;
        while (shell_cum_prob < r1) {
            ++shellIdx;
            shell_cum_prob += material.shell(shellIdx).numberOfElectronsFraction;
        }

        const auto& shell = material.shell(shellIdx);

        if (shell.bindingEnergy > particle.energy) {
            // we reject interaction
            return T { 0 };
        }

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
            rejected = state.randomUniform<T>() > g;
        } while (rejected);

        const auto U = shell.bindingEnergy / ELECTRON_REST_MASS<T>();
        const auto p_i = (k * (k - U) * (1 - cosTheta) - U) / std::sqrt(2 * k * (k - U) * (1 - cosTheta) + U * U);
        constexpr auto p = T { 0.15 };

        const auto kc = k * e;
        const auto qc = std::sqrt(k * k + kc * kc - 2 * k * kc * cosTheta);
        const auto alpha = qc * (1 + kc * (kc - k * cosTheta) / (qc * qc)) / k;

        const auto J0 = shell.HartreeFockOrbital_0;
        const auto exp_b_part = 1 + 2 * J0 * std::abs(p_i);

        const auto exp_b = std::exp(-(T { 0.5 } * (exp_b_part * exp_b_part - 1)));

        auto errf = [exp_b, J0](const T pz) -> T {
            const auto t = T { 1 } / (1 + T { 0.332673 } * (1 + 2 * J0 * std::abs(pz)));
            return T { 1.64872127 } - exp_b * t * (T { 0.34802 } - T { 0.0958798 } * t + T { 0.7478556 } * t * t);
        };
        T Si;
        if (p_i < -p) {
            Si = (1 - alpha * p) * exp_b * T { 0.5 };
        } else if (p_i <= 0) {
            Si = (1 + alpha * p_i) * exp_b * T { 0.5 } - alpha / (4 * J0 * std::numbers::inv_sqrtpi_v<T> * std::numbers::sqrt2_v<T>)*(errf(p) - errf(p_i));
        } else if (p_i < p) {
            Si = 1 - (1 + alpha * p_i) * exp_b * T { 0.5 } - alpha / (4 * J0 * std::numbers::inv_sqrtpi_v<T> * std::numbers::sqrt2_v<T>)*(errf(p) - errf(p_i));
        } else {
            Si = 1 - (1 - alpha * p) * exp_b * T { 0.5 };
        }

        if (state.randomUniform<T>() > Si) {
            return T { 0 };
        }
        T Fmax;
        if (p_i < -p) {
            Fmax = 1 - alpha * p;
        } else if (p_i > p) {
            Fmax = 1 + alpha * p;
        } else {
            Fmax = 1 + alpha * p_i;
        }

        T pz, Fpz;
        do {

            const auto r_b = state.randomUniform(exp_b);
            if (r_b < T { 0.5 }) {
                pz = (1 - std::sqrt(1 - 2 * std::numbers::ln2_v<T> * r_b)) / (2 * J0);
                // pz = (1 - std::sqrt(1 - 2 * std::log(2 * r_b))) / (2 * J0);
            } else {
                pz = (std::sqrt(1 - 2 * std::numbers::ln2_v<T> * (1 - r_b)) - 1) / (2 * J0);
                // pz = (std::sqrt(1 - 2 * std::log(2 * (1 - r_b))) - 1) / (2 * J0);
            }

            if (pz < -p) {
                Fpz = 1 - alpha * p;
            } else if (pz > p) {
                Fpz = 1 + alpha * p;
            } else {
                Fpz = 1 + alpha * pz;
            }
        } while (state.randomUniform<T>(Fmax) > Fpz);

        const auto eb = e / (1 - pz * pz * e * e) * (1 - pz * pz * e * cosTheta + pz * std::sqrt(1 - 2 * e * cosTheta + e * e * (1 - pz * pz * (1 - cosTheta * cosTheta))));

        particle.energy *= eb;
        const auto theta = std::acos(cosTheta);
        const auto phi = state.randomUniform<T>(PI_VAL<T>() + PI_VAL<T>());
        vectormath::peturb<T>(particle.dir, theta, phi);

        return E - particle.energy;
    }

    template <Floating T, int Lowenergycorrection, int Nshells>
    T comptonScatter(Particle<T>& particle, const Material2<T, Nshells>& material, RandomState& state) noexcept
    // see http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/PhysicsReferenceManual/fo/PhysicsReferenceManual.pdf
    // and
    // https://nrc-cnrc.github.io/EGSnrc/doc/pirs701-egsnrc.pdf
    {
        if constexpr (Lowenergycorrection == 2) {
            return comptonScatterIA(particle, material, state);
        } else if constexpr (Lowenergycorrection == 3) {
            return comptonScatterPenelope(particle, material, state);
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
                } else { // Livermore and IA process
                    const auto q = material.momentumTransferCosAngle(particle.energy, cosTheta);
                    const auto scatterFactor = material.scatterFactor(q);
                    // normalize scatterfactor
                    rejected = state.randomUniform<T>(material.effectiveZ()) > (g * scatterFactor);
                }
            } while (rejected);

            const auto theta = std::acos(cosTheta);
            const auto phi = state.randomUniform<T>(PI_VAL<T>() + PI_VAL<T>());
            vectormath::peturb<T>(particle.dir, theta, phi);

            const auto E = particle.energy;
            particle.energy *= e;
            return E - particle.energy;
        }
    }

}
}
