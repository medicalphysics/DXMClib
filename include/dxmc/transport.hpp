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

#include "dxmc/attenuationlut.hpp"
#include "dxmc/dxmcrandom.hpp"
#include "dxmc/exposure.hpp"
#include "dxmc/lowenergycorrectionmodel.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/progressbar.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <functional>
#include <list>
#include <memory>
#include <optional>
#include <thread>

namespace dxmc {

template <Floating T = double>
struct Result {
    std::vector<T> dose;
    std::vector<std::uint_fast32_t> nEvents;
    std::vector<T> variance;
    std::vector<std::uint_fast32_t> nEntered;
    std::uint64_t numberOfHistories { 0 };
    std::chrono::duration<float> simulationTime { 0 };
    std::string_view dose_units = "";

    Result() { }
    Result(std::size_t size)
    {
        dose.resize(size, 0);
        nEvents.resize(size, 0);
        variance.resize(size, 0);
        nEntered.resize(size, 0);
    }
    [[nodiscard]] std::vector<T> relativeError() const
    {
        std::vector<T> err(dose.size());
        std::transform(std::execution::par_unseq, dose.cbegin(), dose.cend(), variance.cbegin(), err.begin(),
            [=](const auto d, const auto v) { return d > 0 ? std::sqrt(v) / d : 0; });
        return err;
    }
};

class EvenTaskPool {
    std::list<std::pair<std::uint64_t, std::uint64_t>> m_intervals;
    std::list<std::pair<std::uint64_t, std::uint64_t>>::iterator m_idx;
    std::mutex m_mutex;

public:
    EvenTaskPool(std::uint64_t start, std::uint64_t stop, std::uint64_t nworkers)
    {
        nworkers = std::max(std::uint64_t { 1 }, nworkers);
        const auto step = std::max((std::max(start, stop) - std::min(start, stop)) / nworkers,
            std::uint64_t { 1 });

        const auto max = std::max(start, stop);
        for (auto min = std::min(start, stop); min < max; min += step) {
            m_intervals.push_back(std::make_pair(min, std::min(min + step, max)));
        }
        m_idx = m_intervals.begin();
    }

    std::optional<std::uint64_t> getJob() noexcept
    {
        std::lock_guard lock(m_mutex);

        if (m_idx == m_intervals.end())
            return {};

        if (m_idx->first >= m_idx->second) {
            m_idx = m_intervals.erase(m_idx);
        } else {
            m_idx++;
        }
        if (m_idx == m_intervals.end()) {
            m_idx = m_intervals.begin();
            if (m_idx == m_intervals.end())
                return {};
        }
        const auto job = (m_idx->first)++;
        return job;
    }
};

// forward declaration
template <Floating T>
class Source;
template <Floating T>
class CTSource;

template <Floating T = double>
class Transport {
public:
    enum class OUTPUTMODE {
        EV_PER_HISTORY,
        DOSE
    };

    Transport()
    {
        std::uint64_t nthreads = std::thread::hardware_concurrency();
        this->m_nThreads = std::max(nthreads, std::uint64_t { 1 });
    }

    void setNumberOfWorkers(std::uint64_t n) { this->m_nThreads = std::max(n, std::uint64_t { 1 }); }
    std::size_t numberOfWorkers() const { return m_nThreads; }
    void setLowEnergyCorrectionModel(LOWENERGYCORRECTION model) { m_lowenergyCorrection = model; }
    LOWENERGYCORRECTION lowEnergyCorrectionModel() const { return m_lowenergyCorrection; }
    void setOutputMode(Transport<T>::OUTPUTMODE mode) { m_outputmode = mode; }
    Transport<T>::OUTPUTMODE outputMode() const { return m_outputmode; }

    template <typename U>
    requires std::is_base_of<World<T>, U>::value
        Result<T>
    operator()(const U& world, Source<T>* source, ProgressBar<T>* progressbar = nullptr, bool useSourceDoseCalibration = true)
    {
        Result<T> result(world.size());

        if (!world.isValid())
            return result;

        if (!source)
            return result;

        source->updateFromWorld(world);
        source->validate();
        if (!source->isValid())
            return result;

        result.numberOfHistories = source->historiesPerExposure() * source->totalExposures();

        const auto maxEnergy = source->maxPhotonEnergyProduced();

        m_attenuationLut.generate(world, maxEnergy);

        const std::uint64_t totalExposures = source->totalExposures();

        const std::uint64_t nJobs = std::max(m_nThreads, std::uint64_t { 1 });
        if (progressbar) {
            progressbar->setTotalExposures(totalExposures);
            progressbar->setDoseData(result.dose.data(), world.dimensions(), world.spacing());
        }
        const auto start = std::chrono::system_clock::now();
        if (m_lowenergyCorrection == LOWENERGYCORRECTION::NONE) {
            parallellRun<0>(world, source, result, 0, totalExposures, nJobs, progressbar);
        } else if (m_lowenergyCorrection == LOWENERGYCORRECTION::LIVERMORE) {
            parallellRun<1>(world, source, result, 0, totalExposures, nJobs, progressbar);
        } else {
            parallellRun<2>(world, source, result, 0, totalExposures, nJobs, progressbar);
        }
        result.simulationTime = std::chrono::system_clock::now() - start;

        if (progressbar) {
            progressbar->clearDoseData();
            if (progressbar->cancel()) {
                std::fill(result.dose.begin(), result.dose.end(), T { 0.0 });
                std::fill(result.nEvents.begin(), result.nEvents.end(), 0);
                std::fill(result.variance.begin(), result.variance.end(), T { 0.0 });
                result.numberOfHistories = 0;
                return result;
            }
        }

        normalizeScoring(result);
        if (m_outputmode == OUTPUTMODE::DOSE) {
            if (useSourceDoseCalibration) {
                const T calibrationValue = source->getCalibrationValue(m_lowenergyCorrection, progressbar);
                // energy imparted to dose
                energyImpartedToDose(world, result, calibrationValue);
                result.dose_units = "mGy";
            } else {
                energyImpartedToDose(world, result);
                result.dose_units = "eV/kg/history";
            }
        } else {
            result.dose_units = "eV/history";
        }
        return result;
    }

    const AttenuationLut<T>& attenuationLut() const { return m_attenuationLut; }
    AttenuationLut<T>& attenuationLut() { return m_attenuationLut; }

protected:
    template <typename U>
    requires std::is_arithmetic_v<U>
    inline void safeValueAdd(U& value, const U addValue) const noexcept
    {
        std::atomic_ref value_l(value);
        value_l.fetch_add(addValue);
    }

    template <int Lowenergycorrection>
    T photoAbsorption(Particle<T>& particle, std::uint8_t materialIdx, RandomState& state) const noexcept
    {
        const auto E = particle.energy;
        particle.energy = 0;
        if constexpr (Lowenergycorrection < 2) {
            return E;
        } else {
            // select shells where binding energy is greater than photon energy
            const auto& elConf = m_attenuationLut.electronShellConfiguration(materialIdx);
            std::array<T, 12> shell_probs;
            std::transform(std::execution::unseq, elConf.cbegin(), elConf.cend(), shell_probs.begin(), [=](const auto& c) { return E > c.bindingEnergy ? c.photoIonizationProbability : 0; });
            // cummulative sum of shell probs
            std::partial_sum(shell_probs.cbegin(), shell_probs.cend(), shell_probs.begin());
            // selectring shell
            std::size_t idx = 0;
            const auto shellIdxSample = shell_probs.back() * state.randomUniform<T>();
            while (shell_probs[idx] < shellIdxSample && idx < 11) {
                ++idx;
            }

            if (elConf[idx].bindingEnergy > E || elConf[idx].bindingEnergy < ENERGY_CUTOFF_THRESHOLD()) {
                // no suitable shell, all energy is deposited locally
                return E;
            }

            // Determine randomly if a fluorscence photon is emitted from corrected yield
            const auto r1 = state.randomUniform<T>();
            if (r1 <= elConf[idx].fluorescenceYield) {
                // select line transition from the three most probable lines
                std::size_t lineIdx = 0;
                auto r3 = state.randomUniform<T>() - elConf[idx].fluorLineProbabilities[lineIdx];
                while (r3 > 0 && lineIdx < 2) {
                    ++lineIdx;
                    r3 -= elConf[idx].fluorLineProbabilities[lineIdx];
                }
                // give the "new" particle the line energy and an random isotropic direction
                particle.energy = elConf[idx].fluorLineEnergies[lineIdx];
                const auto theta = state.randomUniform<T>(PI_VAL<T>());
                const auto phi = state.randomUniform<T>(PI_VAL<T>() + PI_VAL<T>());
                vectormath::peturb(particle.dir, theta, phi);
                // return energy imparted locally
                return E - particle.energy;
            } else {
                return E;
            }
        }
    }
    template <int Lowenergycorrection>
    void rayleightScatter(Particle<T>& particle, std::uint8_t materialIdx, RandomState& state) const noexcept
    {
        if constexpr (Lowenergycorrection == 0) {
            bool reject = true;
            T theta;
            while (reject) {
                constexpr T extreme = (4 * std::numbers::sqrt2_v<T>) / (3 * std::numbers::sqrt3_v<T>);
                const auto r1 = state.randomUniform<T>(T { 0 }, extreme);
                theta = state.randomUniform(T { 0 }, PI_VAL<T>());
                const auto sinang = std::sin(theta);
                reject = r1 > ((2 - sinang * sinang) * sinang);
            }
            // calc angle and add randomly 90 degrees since dist i symetrical
            const auto phi = state.randomUniform<T>(PI_VAL<T>() + PI_VAL<T>());
            vectormath::peturb<T>(particle.dir, theta, phi);
        } else {
            // theta is scattering angle
            // see http://rcwww.kek.jp/research/egs/egs5_manual/slac730-150228.pdf

            // finding qmax
            const auto qmax = m_attenuationLut.momentumTransferMax(particle.energy);
            const auto qmax_squared = qmax * qmax;

            T cosAngle;
            do {
                const auto q_squared = m_attenuationLut.momentumTransferFromFormFactor(materialIdx, qmax_squared, state);
                cosAngle = m_attenuationLut.cosAngle(particle.energy, q_squared);
            } while ((1 + cosAngle * cosAngle) * T { 0.5 } < state.randomUniform<T>());

            const auto theta = std::acos(cosAngle);
            const auto phi = state.randomUniform<T>(PI_VAL<T>() + PI_VAL<T>());
            vectormath::peturb<T>(particle.dir, theta, phi);
        }
    }

    template <int Lowenergycorrection>
    T comptonScatter(Particle<T>& particle, std::uint8_t materialIdx, RandomState& state) const noexcept
    // see http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/PhysicsReferenceManual/fo/PhysicsReferenceManual.pdf
    // and
    // https://nrc-cnrc.github.io/EGSnrc/doc/pirs701-egsnrc.pdf
    {
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
        }
    }

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

    inline bool particleInsideWorld(const World<T>& world, const std::array<T, 3>& pos) const noexcept
    {
        const auto& extent = world.matrixExtentSafe();
        return (pos[0] > extent[0] && pos[0] < extent[1]) && (pos[1] > extent[2] && pos[1] < extent[3]) && (pos[2] > extent[4] && pos[2] < extent[5]);
    }

    inline bool particleInsideWorld(const World<T>& world, const Particle<T>& particle) const noexcept
    {
        return particleInsideWorld(world, particle.pos);
    }

    inline std::array<std::size_t, 3> gridIndexFromPosition(const std::array<T, 3>& pos, const World<T>& world) const noexcept
    {
        // assumes particle is inside world
        const auto& wpos = world.matrixExtentSafe();
        const auto& wspac = world.spacing();
        std::array<std::size_t, 3> arraypos = {
            static_cast<std::size_t>((pos[0] - wpos[0]) / wspac[0]),
            static_cast<std::size_t>((pos[1] - wpos[2]) / wspac[1]),
            static_cast<std::size_t>((pos[2] - wpos[4]) / wspac[2]),
        };
        return arraypos;
    }

    inline std::size_t indexFromPosition(const std::array<T, 3>& pos, const World<T>& world) const noexcept
    {
        // assumes particle is inside world
        const auto& wdim = world.dimensions();
        const auto arraypos = gridIndexFromPosition(pos, world);
        const auto idx = arraypos[2] * wdim[0] * wdim[1] + arraypos[1] * wdim[0] + arraypos[0];
        return idx;
    }

    inline std::size_t indexFromPosition(const Particle<T>& particle, const World<T>& world) const noexcept
    {
        return indexFromPosition(particle.pos, world);
    }

    template <int Lowenergycorrection>
    bool computeInteractionsForced(const T eventProbability, const std::array<T, 3>& attenuation, Particle<T>& p, const std::uint8_t matIdx, Result<T>& result, const std::size_t resultBufferIdx, RandomState& state, bool& updateMaxAttenuation) const noexcept
    {
        const auto attPhoto = attenuation[0];
        const auto attCompt = attenuation[1];
        const auto attRayl = attenuation[2];
        const auto attenuationTotal = std::reduce(std::execution::unseq, attenuation.cbegin(), attenuation.cend(), T { 0 });

        const auto photoEventProbability = attPhoto / attenuationTotal;
        const auto weightCorrection = eventProbability * photoEventProbability;
        {
            // virtual event, we ignore characteristic radiation in case of forced events (this may slightly overestimate dose in high Z materials)
            auto p_forced = p; // make a copyin case of characteristic x-ray
            const auto e_forced = photoAbsorption<Lowenergycorrection>(p_forced, matIdx, state);
            if (p_forced.energy < ENERGY_CUTOFF_THRESHOLD())
                [[unlikely]] {
                const auto energyImparted = (e_forced + p_forced.energy) * p_forced.weight * weightCorrection;
                safeValueAdd(result.dose[resultBufferIdx], energyImparted);
                safeValueAdd(result.nEvents[resultBufferIdx], std::uint_fast32_t { 1 });
                safeValueAdd(result.variance[resultBufferIdx], energyImparted * energyImparted);
            } else
                [[likely]] {
                const auto energyImparted = e_forced * p_forced.weight * weightCorrection;
                safeValueAdd(result.dose[resultBufferIdx], energyImparted);
                safeValueAdd(result.nEvents[resultBufferIdx], std::uint_fast32_t { 1 });
                safeValueAdd(result.variance[resultBufferIdx], energyImparted * energyImparted);
            }
        }
        const auto r1 = state.randomUniform<T>();
        if (r1 < eventProbability * (1 - photoEventProbability)) {
            // A real event happend
            const auto r2 = state.randomUniform(attCompt + attRayl);
            if (r2 < attCompt) // Compton event
            {
                const auto e = comptonScatter<Lowenergycorrection>(p, matIdx, state);
                if (p.energy < ENERGY_CUTOFF_THRESHOLD())
                    [[unlikely]] {
                    const auto energyImparted = (e + p.energy) * p.weight;
                    safeValueAdd(result.dose[resultBufferIdx], energyImparted);
                    safeValueAdd(result.nEvents[resultBufferIdx], std::uint_fast32_t { 1 });
                    safeValueAdd(result.variance[resultBufferIdx], energyImparted * energyImparted);
                    p.energy = 0;
                    return false;
                } else
                    [[likely]] {
                    const auto energyImparted = e * p.weight;
                    safeValueAdd(result.dose[resultBufferIdx], energyImparted);
                    safeValueAdd(result.nEvents[resultBufferIdx], std::uint_fast32_t { 1 });
                    safeValueAdd(result.variance[resultBufferIdx], energyImparted * energyImparted);
                    updateMaxAttenuation = true;
                }
            } else // Rayleigh scatter event
            {
                rayleightScatter<Lowenergycorrection>(p, matIdx, state);
            }
        }
        p.weight *= (1 - weightCorrection); // to prevent bias (weightcorrection is the real probability of a photelectric event)
        return true;
    }

    template <int Lowenergycorrection>
    bool computeInteractions(const std::array<T, 3>& attenuation, Particle<T>& p, const std::uint8_t matIdx, Result<T>& result, const std::size_t resultBufferIdx, RandomState& state, bool& updateMaxAttenuation) const noexcept
    {
        const auto attPhoto = attenuation[0];
        const auto attCompt = attenuation[1];
        const auto attRayl = attenuation[2];

        const auto attenuationTotal = std::reduce(std::execution::unseq, attenuation.cbegin(), attenuation.cend(), T { 0 });

        const auto r3 = state.randomUniform(attenuationTotal);
        if (r3 < attPhoto) // Photoelectric event
        {
            const auto e = photoAbsorption<Lowenergycorrection>(p, matIdx, state);
            if (p.energy < ENERGY_CUTOFF_THRESHOLD())
                [[unlikely]] {
                const auto energyImparted = (e + p.energy) * p.weight;
                safeValueAdd(result.dose[resultBufferIdx], energyImparted);
                safeValueAdd(result.nEvents[resultBufferIdx], std::uint_fast32_t { 1 });
                safeValueAdd(result.variance[resultBufferIdx], energyImparted * energyImparted);
                p.energy = 0;
                return false;
            } else
                [[likely]] {
                const auto energyImparted = e * p.weight;
                safeValueAdd(result.dose[resultBufferIdx], energyImparted);
                safeValueAdd(result.nEvents[resultBufferIdx], std::uint_fast32_t { 1 });
                safeValueAdd(result.variance[resultBufferIdx], energyImparted * energyImparted);
                updateMaxAttenuation = true;
            }

        } else if (r3 < (attPhoto + attCompt)) // Compton event
        {
            const auto e = comptonScatter<Lowenergycorrection>(p, matIdx, state);
            if (p.energy < ENERGY_CUTOFF_THRESHOLD())
                [[unlikely]] {
                const auto energyImparted = (e + p.energy) * p.weight;
                safeValueAdd(result.dose[resultBufferIdx], energyImparted);
                safeValueAdd(result.nEvents[resultBufferIdx], std::uint_fast32_t { 1 });
                safeValueAdd(result.variance[resultBufferIdx], energyImparted * energyImparted);
                p.energy = 0;
                return false;
            } else
                [[likely]] {
                const auto energyImparted = e * p.weight;
                safeValueAdd(result.dose[resultBufferIdx], energyImparted);
                safeValueAdd(result.nEvents[resultBufferIdx], std::uint_fast32_t { 1 });
                safeValueAdd(result.variance[resultBufferIdx], energyImparted * energyImparted);
                updateMaxAttenuation = true;
            }
        } else // Rayleigh scatter event
        {
            rayleightScatter<Lowenergycorrection>(p, matIdx, state);
        }

        return true;
    }

    template <int Lowenergycorrection>
    void woodcockParticleTracking(const World<T>& world, Particle<T>& p, RandomState& state, Result<T>& result) const noexcept
    {
        const auto densityBuffer = world.densityArray()->data();
        const auto materialBuffer = world.materialIndexArray()->data();
        const auto measurementBuffer = world.measurementMapArray()->data();

        T maxAttenuationInv;
        bool updateMaxAttenuation = true;
        bool continueSampling = true;
        while (continueSampling) {
            if (updateMaxAttenuation) {
                maxAttenuationInv = m_attenuationLut.maxTotalAttenuationInverse(p.energy);
                updateMaxAttenuation = false;
            }
            const auto r1 = state.randomUniform<T>();
            const auto stepLenght = -std::log(r1) * maxAttenuationInv * T { 10 }; // cm -> mm
            for (std::size_t i = 0; i < 3; i++)
                p.pos[i] += p.dir[i] * stepLenght;

            if (particleInsideWorld(world, p)) {
                const std::size_t bufferIdx = indexFromPosition(p, world);
                const auto matIdx = materialBuffer[bufferIdx];
                const auto density = densityBuffer[bufferIdx];
                const auto measurement = measurementBuffer[bufferIdx];

                const auto attenuation = m_attenuationLut.photoComptRayAttenuation(matIdx, p.energy);
                const auto attenuationTotal = std::reduce(std::execution::unseq, attenuation.cbegin(), attenuation.cend(), T { 0 }) * density;
                const auto eventProbability = attenuationTotal * maxAttenuationInv;

                safeValueAdd(result.nEntered[bufferIdx], std::uint_fast32_t { 1 });

                if (measurement == 0)
                    [[likely]] // naive sampling
                {
                    const auto r2 = state.randomUniform<T>();
                    if (r2 < eventProbability) // an event will happen
                    {
                        continueSampling = computeInteractions<Lowenergycorrection>(attenuation, p, matIdx, result, bufferIdx, state, updateMaxAttenuation);
                    }
                } else
                    [[unlikely]] // forced photoelectric effect
                {
                    continueSampling = computeInteractionsForced<Lowenergycorrection>(eventProbability, attenuation, p, matIdx, result, bufferIdx, state, updateMaxAttenuation);
                }

                if (continueSampling) {
                    if ((p.energy * p.weight < RUSSIAN_RULETTE_ENERGY_THRESHOLD())) {
                        const auto r4 = state.randomUniform<T>();
                        if (r4 < RUSSIAN_RULETTE_PROBABILITY()) {
                            continueSampling = false;
                        } else {
                            constexpr T factor = T { 1 } / (T { 1 } - RUSSIAN_RULETTE_PROBABILITY());
                            p.weight *= factor;
                        }
                    }
                }
            } else // particle not inside world
            {
                continueSampling = false;
            }
        }
    }

    bool transportParticleToWorld(const World<T>& world, Particle<T>& particle) const noexcept
    // returns false if particle do not intersect world
    {
        const bool isInside = particleInsideWorld(world, particle);
        if (isInside)
            return true;

        auto amin = std::numeric_limits<T>::min();
        auto amax = std::numeric_limits<T>::max();

        const auto& extent = world.matrixExtentSafe();
        for (std::size_t i = 0; i < 3; i++) {
            if (std::abs(particle.dir[i]) > N_ERROR()) {
                const auto a0 = (extent[i * 2] - particle.pos[i]) / particle.dir[i];
                const auto an = (extent[i * 2 + 1] - particle.pos[i]) / particle.dir[i];
                amin = std::max(amin, std::min(a0, an));
                amax = std::min(amax, std::max(a0, an));
            }
        }
        if (amin < amax && amin > 0) {
            for (std::size_t i = 0; i < 3; i++) {
                particle.pos[i] += amin * particle.dir[i]; // making sure particle is inside world
            }
            return true;
        }
        return false;
    }
    template <int Lowenergycorrection>
    void transport(const World<T>& world, const Exposure<T>& exposure, RandomState& state, Result<T>& result) const noexcept
    {
        const auto nHistories = exposure.numberOfHistories();
        for (std::size_t i = 0; i < nHistories; ++i) {
            // Draw a particle
            auto particle = exposure.sampleParticle(state);

            // Is particle intersecting with world
            const bool isInWorld = transportParticleToWorld(world, particle);
            if (isInWorld) {
                woodcockParticleTracking<Lowenergycorrection>(world, particle, state, result);
            }
        }
    }

    template <int Lowenergycorrection>
    void workerRun(const World<T>& w, const Source<T>* source, Result<T>& result,
        EvenTaskPool& taskpool, ProgressBar<T>* progressbar) const noexcept
    {
        RandomState state;
        const auto& worldBasis = w.directionCosines();
        auto job = taskpool.getJob();
        while (job) {
            if (progressbar)
                if (progressbar->cancel())
                    return;
            auto exposure = source->getExposure(job.value());
            exposure.alignToDirectionCosines(worldBasis);
            transport<Lowenergycorrection>(w, exposure, state, result);
            if (progressbar)
                progressbar->exposureCompleted();
            job = taskpool.getJob();
        }
    }

    template <int Lowenergycorrection>
    void parallellRun(const World<T>& w, const Source<T>* source, Result<T>& result,
        const std::uint64_t expBeg, const std::uint64_t expEnd, std::uint64_t nJobs, ProgressBar<T>* progressbar) const
    {
        EvenTaskPool taskpool(expBeg, expEnd, nJobs);
        std::vector<std::thread> jobs;
        jobs.reserve(nJobs - 1);
        for (std::size_t i = 0; i < nJobs - 1; ++i) {
            jobs.emplace_back(&Transport<T>::workerRun<Lowenergycorrection>, this, std::cref(w), source, std::ref(result), std::ref(taskpool), progressbar);
        }
        workerRun<Lowenergycorrection>(w, source, result, taskpool, progressbar);
        for (auto& job : jobs)
            job.join();
    }

    void normalizeScoring(Result<T>& res) const noexcept
    {
        // Here we normalize dose and precision to eV per history
        // Computing precision by: Var[X] = E[X**2] - E[X]**2
        // and Var[A]+ Var[A] = 2Var[A]

        const auto h_inv = T { 1 } / (res.numberOfHistories - 1);
        const auto h_d_inv = T { 1E3 } / res.numberOfHistories;
        const auto h_v_inv = T { 1E6 } / res.numberOfHistories;
        std::transform(std::execution::par_unseq, res.dose.cbegin(), res.dose.cend(), res.nEvents.cbegin(), res.dose.begin(),
            [=](const auto d, const auto n) { return d * h_d_inv; });

        std::transform(std::execution::par_unseq, res.dose.cbegin(), res.dose.cend(), res.variance.cbegin(), res.variance.begin(),
            [=](const auto d, const auto dd) -> T { return (dd * h_v_inv - d * d) * h_inv; });
    }

    void energyImpartedToDose(const World<T>& world, Result<T>& res, const T calibrationValue = 1) noexcept
    {
        const auto& spacing = world.spacing();
        const auto voxelVolume = spacing[0] * spacing[1] * spacing[2] / T { 1000 }; // cm3
        auto density = world.densityArray();

        std::transform(
            std::execution::par_unseq, res.dose.cbegin(), res.dose.cend(), density->cbegin(), res.dose.begin(),
            [=](auto ei, auto de) -> auto{
                const auto voxelMass = de * voxelVolume * T { 0.001 }; // g/cm3 * cm3 / 1000 = kg
                const auto factor = calibrationValue / voxelMass;
                return de > T { 0.0 } ? ei * factor : T { 0.0 };
            });
        std::transform(
            std::execution::par_unseq, res.variance.cbegin(), res.variance.cend(), density->cbegin(), res.variance.begin(),
            [=](auto var, auto de) -> auto{
                const auto voxelMass = de * voxelVolume * T { 0.001 }; // g/cm3 * cm3 / 1000 = kg
                const auto factor = calibrationValue / voxelMass;
                return de > T { 0.0 } ? var * factor * factor : T { 0.0 }; // for variance we must multiply by square
            });
    }

private:
    static constexpr T RUSSIAN_RULETTE_PROBABILITY()
    {
        return T { 0.8 };
    }
    static constexpr T RUSSIAN_RULETTE_ENERGY_THRESHOLD() // keV
    {
        return T { 5 };
    }
    static constexpr T ENERGY_CUTOFF_THRESHOLD() // keV
    {
        return T { 1 };
    }
    static constexpr T N_ERROR()
    {
        return T { 1.0e-9 };
    }

    AttenuationLut<T> m_attenuationLut;
    std::uint64_t m_nThreads;
    OUTPUTMODE m_outputmode = OUTPUTMODE::DOSE;
    LOWENERGYCORRECTION m_lowenergyCorrection = LOWENERGYCORRECTION::LIVERMORE;
};
}