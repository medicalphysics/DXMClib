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

#include "dxmc/attenuationlut.h"
#include "dxmc/dxmcrandom.h"
#include "dxmc/exposure.h"
#include "dxmc/lowenergycorrectionmodel.h"
#include "dxmc/particle.h"
#include "dxmc/progressbar.h"
#include "dxmc/vectormath.h"
#include "dxmc/world.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <functional>
#include <memory>
#include <optional>

namespace dxmc {

template <Floating T = double>
struct Result {
    std::vector<T> dose;
    std::vector<std::uint32_t> nEvents;
    std::vector<T> variance;
    std::uint64_t numberOfHistories { 0 };
    std::chrono::duration<float> simulationTime { 0 };
    std::string_view dose_units = "";

    Result() {}
    Result(std::size_t size)
    {
        dose.resize(size, 0);
        nEvents.resize(size, 0);
        variance.resize(size, 0);
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
                //energy imparted to dose
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

        if constexpr (Lowenergycorrection < 2) {
            particle.energy = 0;
            return E;
        } else {
            // select shells where binding energy is greater than photon energy
            const auto& elConf = m_attenuationLut.electronShellConfiguration(materialIdx);
            std::size_t idx = 0;
            while (elConf[idx].bindingEnergy > E && idx < 11) {
                idx++;
            }
            if (elConf[idx].bindingEnergy > E || elConf[idx].bindingEnergy > ENERGY_CUTOFF_THRESHOLD()) {
                // no suitable shell, all energy is deposited locally
                particle.energy = 0;
                return E;
            }

            //select shell based on partial photoionization crossection
            auto r1 = state.randomUniform<T>();
            while (r1 >= elConf[idx].photoIonizationProbability && idx < 11) {
                r1 = (1 - r1) / (1 - elConf[idx].photoIonizationProbability);
                ++idx;
            }
            if (r1 >= elConf[idx].photoIonizationProbability) {
                particle.energy = 0;
                return E;
            }

            //Determine randomly if a fluorscence photon is emitted from corrected yield
            const auto r2 = state.randomUniform<T>();
            if (r2 <= elConf[idx].fluorescenceYield) {
                // select line transition from the three most probable lines
                std::size_t lineIdx = 0;
                auto r3 = state.randomUniform<T>() - elConf[idx].fluorLineProbabilities[lineIdx];
                while (r3 > 0) {
                    ++lineIdx;
                    r3 -= elConf[idx].fluorLineProbabilities[lineIdx];
                }
                // give the "new" particle the line energy and an random isotropic direction
                particle.energy = elConf[idx].fluorLineEnergies[lineIdx];
                const auto theta = state.randomUniform<T>(PI_VAL<T>());
                const auto phi = state.randomUniform<T>(PI_VAL<T>() + PI_VAL<T>());
                vectormath::peturb(particle.dir.data(), theta, phi);
                //return energy imparted locally
                return E - particle.energy;
            } else {
                particle.energy = 0;
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
            vectormath::peturb<T>(particle.dir.data(), theta, phi);
        } else {
            // theta is scattering angle
            // see http://rcwww.kek.jp/research/egs/egs5_manual/slac730-150228.pdf

            //finding qmax
            const auto qmax = m_attenuationLut.momentumTransferMax(particle.energy);
            const auto qmax_squared = qmax * qmax;

            T cosAngle;
            do {
                const auto q_squared = m_attenuationLut.momentumTransferFromFormFactor(materialIdx, qmax_squared, state);
                cosAngle = m_attenuationLut.cosAngle(particle.energy, q_squared);
            } while ((1 + cosAngle * cosAngle) * T { 0.5 } < state.randomUniform<T>());

            const auto theta = std::acos(cosAngle);
            const auto phi = state.randomUniform<T>(PI_VAL<T>() + PI_VAL<T>());
            vectormath::peturb<T>(particle.dir.data(), theta, phi);
        }
    }

    template <int Lowenergycorrection>
    T comptonScatter(Particle<T>& particle, std::uint8_t materialIdx, RandomState& state) const noexcept
    // see http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/PhysicsReferenceManual/fo/PhysicsReferenceManual.pdf
    // and
    // https://nrc-cnrc.github.io/EGSnrc/doc/pirs701-egsnrc.pdf
    {
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

            if constexpr (Lowenergycorrection == 1) {
                //simple correction with scatterfactor to supress forward scattering due to binding energy
                const auto q = m_attenuationLut.momentumTransferFromCos(E, cosAngle);
                const auto scatterFactor = m_attenuationLut.comptonScatterFactor(materialIdx, q);
                const auto r2 = state.randomUniform<T>();
                const auto g = (1 / e + e - sinthetasqr) * gmax_inv;
                rejected = r2 > g * scatterFactor;
            } else {
                const auto g = (1 / e + e - sinthetasqr) * gmax_inv;
                const auto r2 = state.randomUniform<T>();
                rejected = r2 > g;
            }
            if constexpr (Lowenergycorrection == 2) {
                rejected = rejected || (t == T { 0 });
                // See penelope-2018__a_code_system_for_monte_carlo_simulation_of_electron_and_photon_transport
                // for details
                if (!rejected) {
                    //selecting a electron shell by random, should be improved?
                    bool acceptshell;
                    T pz;
                    const auto EC = E * e;
                    const auto cqc = std::sqrt(E * E + EC * EC - 2 * E * EC * cosAngle);
                    const auto F_pz_part = cqc * (1 + EC * (EC - E * cosAngle) / (cqc * cqc)) / E;
                    const auto F_zp_max = std::max(1 + F_pz_part * T { 0.2 }, 1 - F_pz_part * T { 0.2 });

                    const auto& electronConfigurations = m_attenuationLut.electronShellConfiguration(materialIdx);
                    std::array<T, 12> pz_max, n_pz_max, shell_probs;
                    //calculating pz_max and n_pz_max for all shells
                    for (std::size_t i = 0; i < pz_max.size(); ++i) {
                        const auto U = electronConfigurations[i].bindingEnergy;
                        const auto cq = std::sqrt(E * E + (E - U) * (E - U) - 2 * E * (E - U) * cosAngle);
                        pz_max[i] = E * (E - U - EC) / (EC * cq) * ELECTRON_REST_MASS<T>();
                        //equation 2.49
                        constexpr T d2 = std::numbers::sqrt2_v<T>;
                        constexpr T d1 = T { 1 } / d2;
                        const auto J = electronConfigurations[i].hartreeFockOrbital_0;
                        if (pz_max[i] < 0) {
                            const auto nom = (d1 - d2 * J * pz_max[i]);
                            n_pz_max[i] = T { 0.5 } * std::exp(d1 * d1 - nom * nom);
                        } else {
                            const auto nom = (d1 + d2 * J * pz_max[i]);
                            n_pz_max[i] = 1 - T { 0.5 } * std::exp(d1 * d1 - nom * nom);
                        }
                    }
                    // calculating shell probability
                    std::transform(std::execution::unseq,
                        n_pz_max.cbegin(), n_pz_max.cend(), electronConfigurations.cbegin(), shell_probs.begin(),
                        [=](const auto n_pz, const auto& c) -> T { return E > c.bindingEnergy && c.hartreeFockOrbital_0 > T { 0 } ? n_pz * c.numberElectrons : 0; });
                    const auto shell_probs_sum_inv = 1 / std::reduce(std::execution::unseq, shell_probs.cbegin(), shell_probs.cend());
                    std::transform(std::execution::unseq, shell_probs.cbegin(), shell_probs.cend(), shell_probs.begin(), [=](const auto p) { return p * shell_probs_sum_inv; });
                    std::partial_sum(shell_probs.cbegin(), shell_probs.cend(), shell_probs.begin());
                    
                    do {
                        //sampling electron shell
                        const auto r3 = state.randomUniform<T>();
                        const auto shell_pos = std::upper_bound(shell_probs.cbegin(), shell_probs.cend(), r3);
                        //const auto electronIdx = std::min(std::distance(shell_probs.cbegin(), shell_pos), 11);
                        const auto electronIdx = std::distance(shell_probs.cbegin(), shell_pos);

                        //sampling pz analytically
                        const auto r4 = state.randomUniform<T>();
                        const auto A = r4 * n_pz_max[electronIdx];
                        // eq 2.59
                        if (A < T { 0.5 }) {
                            constexpr T ln2 = std::numbers::ln2_v<T>;
                            constexpr T d2 = std::numbers::sqrt2_v<T>;
                            constexpr T d1 = T { 1 } / d2;
                            pz = (d1 - std::sqrt(d1 * d1 - ln2 * A)) / (d2 * electronConfigurations[electronIdx].hartreeFockOrbital_0);
                        } else {
                            constexpr T ln2 = std::numbers::ln2_v<T>;
                            constexpr T d2 = std::numbers::sqrt2_v<T>;
                            constexpr T d1 = T { 1 } / d2;
                            pz = (std::sqrt(d1 * d1 - ln2 * (1 - A)) - d1) / (d2 * electronConfigurations[electronIdx].hartreeFockOrbital_0);
                        }
                        //normalizing pz to atomic units
                        pz *= (1 / ELECTRON_REST_MASS<T>());

                        //pz must be larger than  minus electron rest mass
                        acceptshell = pz > -1;

                        //Sampling
                        const auto F_zp = 1 + F_pz_part * pz;
                        const auto r5 = state.randomUniform<T>();
                        acceptshell = acceptshell && r5 * F_zp_max > F_zp;

                    } while (!acceptshell);
                    // calculating doppler shifted e
                    const auto t_s = pz * pz;
                    const auto nom = (1 - t_s * e * cosAngle);
                    if (pz > 0) {
                        e = e / (1 - t_s * e * e) * (nom + std::sqrt(nom * nom - (1 - t_s * e * e) * (1 - t_s)));
                    } else {
                        e = e / (1 - t_s * e * e) * (nom - std::sqrt(nom * nom - (1 - t_s * e * e) * (1 - t_s)));
                    }
                    Ee = E * (1 - e);
                }
            }

        } while (rejected);
        const auto theta = std::acos(cosAngle);
        const auto phi = state.randomUniform<T>(PI_VAL<T>() + PI_VAL<T>());
        vectormath::peturb<T>(particle.dir.data(), theta, phi);

        particle.energy *= e;

        if constexpr (Lowenergycorrection < 2) {
            Ee = E * (1 - e);
        }
        return Ee;
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
        //assumes particle is inside world
        const auto& wpos = world.matrixExtent();
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
        //assumes particle is inside world
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
        const auto attenuationTotal = attPhoto + attCompt + attRayl;

        const auto weightCorrection = eventProbability * attPhoto / attenuationTotal;
        //copy of particle for forced event

        const auto r2 = state.randomUniform<T>();
        if (r2 < eventProbability) {
            //A real event happend
            const auto r3 = state.randomUniform(attenuationTotal);
            if (r3 <= attPhoto) // Photoelectric event
            {
                const auto e = photoAbsorption<Lowenergycorrection>(p, matIdx, state);
                if (p.energy < ENERGY_CUTOFF_THRESHOLD())
                    [[likely]]
                    {
                        const auto energyImparted = (e + p.energy) * p.weight;
                        safeValueAdd(result.dose[resultBufferIdx], energyImparted);
                        safeValueAdd(result.nEvents[resultBufferIdx], std::uint32_t { 1 });
                        safeValueAdd(result.variance[resultBufferIdx], energyImparted * energyImparted);
                        p.energy = 0;
                        return false;
                    }
                else
                    [[unlikely]]
                    {
                        const auto energyImparted = e * p.weight;
                        safeValueAdd(result.dose[resultBufferIdx], energyImparted);
                        safeValueAdd(result.nEvents[resultBufferIdx], std::uint32_t { 1 });
                        safeValueAdd(result.variance[resultBufferIdx], energyImparted * energyImparted);
                        updateMaxAttenuation = true;
                    }

            } else if (r3 <= attPhoto + attCompt) // Compton event
            {
                const auto e = comptonScatter<Lowenergycorrection>(p, matIdx, state);
                if (p.energy < ENERGY_CUTOFF_THRESHOLD())
                    [[unlikely]]
                    {
                        const auto energyImparted = (e + p.energy) * p.weight;
                        safeValueAdd(result.dose[resultBufferIdx], energyImparted);
                        safeValueAdd(result.nEvents[resultBufferIdx], std::uint32_t { 1 });
                        safeValueAdd(result.variance[resultBufferIdx], energyImparted * energyImparted);
                        p.energy = 0;
                        return false;
                    }
                else
                    [[likely]]
                    {
                        const auto energyImparted = e * p.weight;
                        safeValueAdd(result.dose[resultBufferIdx], energyImparted);
                        safeValueAdd(result.nEvents[resultBufferIdx], std::uint32_t { 1 });
                        safeValueAdd(result.variance[resultBufferIdx], energyImparted * energyImparted);
                        updateMaxAttenuation = true;
                    }
            } else // Rayleigh scatter event
            {
                rayleightScatter<Lowenergycorrection>(p, matIdx, state);
            }

        } else {
            // virtual event, we ignore characteristic radiation in case of forced events (this may slightly overestimate dose in high Z materials)
            auto p_forced = p; // make a copyin case of characteristic x-ray
            p.weight *= (1 - weightCorrection); // to prevent bias (weightcorrection is the real probability of a photelectric event)
            const auto e = photoAbsorption<Lowenergycorrection>(p_forced, matIdx, state);
            if (p_forced.energy < ENERGY_CUTOFF_THRESHOLD())
                [[likely]]
                {
                    const auto energyImparted = (e + p_forced.energy) * p_forced.weight * weightCorrection;
                    safeValueAdd(result.dose[resultBufferIdx], energyImparted);
                    safeValueAdd(result.nEvents[resultBufferIdx], std::uint32_t { 1 });
                    safeValueAdd(result.variance[resultBufferIdx], energyImparted * energyImparted);
                }
            else
                [[unlikely]]
                {
                    const auto energyImparted = e * p_forced.weight * weightCorrection;
                    safeValueAdd(result.dose[resultBufferIdx], energyImparted);
                    safeValueAdd(result.nEvents[resultBufferIdx], std::uint32_t { 1 });
                    safeValueAdd(result.variance[resultBufferIdx], energyImparted * energyImparted);
                }
        }
        return true;
    }
    template <int Lowenergycorrection>
    bool computeInteractions(const std::array<T, 3>& attenuation, Particle<T>& p, const std::uint8_t matIdx, Result<T>& result, const std::size_t resultBufferIdx, RandomState& state, bool& updateMaxAttenuation) const noexcept
    {
        const auto attPhoto = attenuation[0];
        const auto attCompt = attenuation[1];
        const auto attRayl = attenuation[2];

        const auto r3 = state.randomUniform(attPhoto + attCompt + attRayl);
        if (r3 <= attPhoto) // Photoelectric event
        {
            const auto e = photoAbsorption<Lowenergycorrection>(p, matIdx, state);
            if (p.energy < ENERGY_CUTOFF_THRESHOLD())
                [[likely]]
                {
                    const auto energyImparted = (e + p.energy) * p.weight;
                    safeValueAdd(result.dose[resultBufferIdx], energyImparted);
                    safeValueAdd(result.nEvents[resultBufferIdx], std::uint32_t { 1 });
                    safeValueAdd(result.variance[resultBufferIdx], energyImparted * energyImparted);
                    p.energy = 0;
                    return false;
                }
            else
                [[unlikely]]
                {
                    const auto energyImparted = e * p.weight;
                    safeValueAdd(result.dose[resultBufferIdx], energyImparted);
                    safeValueAdd(result.nEvents[resultBufferIdx], std::uint32_t { 1 });
                    safeValueAdd(result.variance[resultBufferIdx], energyImparted * energyImparted);
                    updateMaxAttenuation = true;
                }

        } else if (r3 <= attPhoto + attCompt) // Compton event
        {
            const auto e = comptonScatter<Lowenergycorrection>(p, matIdx, state);
            if (p.energy < ENERGY_CUTOFF_THRESHOLD())
                [[unlikely]]
                {
                    const auto energyImparted = (e + p.energy) * p.weight;
                    safeValueAdd(result.dose[resultBufferIdx], energyImparted);
                    safeValueAdd(result.nEvents[resultBufferIdx], std::uint32_t { 1 });
                    safeValueAdd(result.variance[resultBufferIdx], energyImparted * energyImparted);
                    p.energy = 0;
                    return false;
                }
            else
                [[likely]]
                {
                    const auto energyImparted = e * p.weight;
                    safeValueAdd(result.dose[resultBufferIdx], energyImparted);
                    safeValueAdd(result.nEvents[resultBufferIdx], std::uint32_t { 1 });
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

                const auto attenuation = m_attenuationLut.photoComptRayAttenuation(matIdx, p.energy);
                const auto attenuationTotal = std::reduce(std::execution::unseq, attenuation.cbegin(), attenuation.cend(), T { 0 }) * density;
                const auto eventProbability = attenuationTotal * maxAttenuationInv;

                if (measurementBuffer[bufferIdx] == 0)
                    [[likely]] // naive sampling
                    {
                        const auto r2 = state.randomUniform<T>();
                        if (r2 < eventProbability) // an event will happen
                        {
                            continueSampling = computeInteractions<Lowenergycorrection>(attenuation, p, matIdx, result, bufferIdx, state, updateMaxAttenuation);
                        }
                    }
                else
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
    //returns false if particle do not intersect world
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
            //Draw a particle
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
        //Here we normalize dose and precision to eV per history
        //Computing precision by: Var[X] = E[X**2] - E[X]**2
        //and Var[A]+ Var[A] = 2Var[A]

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
            [=](auto ei, auto de) -> auto {
                const auto voxelMass = de * voxelVolume * T { 0.001 }; //g/cm3 * cm3 / 1000 = kg
                const auto factor = calibrationValue / voxelMass;
                return de > T { 0.0 } ? ei * factor : T { 0.0 };
            });
        std::transform(
            std::execution::par_unseq, res.variance.cbegin(), res.variance.cend(), density->cbegin(), res.variance.begin(),
            [=](auto var, auto de) -> auto {
                const auto voxelMass = de * voxelVolume * T { 0.001 }; //g/cm3 * cm3 / 1000 = kg
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