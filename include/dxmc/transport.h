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

    Result() {};
    Result(std::size_t size)
    {
        dose.resize(size, 0);
        nEvents.resize(size, 0);
        variance.resize(size, 0);
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
    Transport()
    {
        std::uint64_t nthreads = std::thread::hardware_concurrency();
        this->m_nThreads = std::max(nthreads, std::uint64_t { 1 });
    }

    void setNumberOfWorkers(std::uint64_t n) { this->m_nThreads = std::max(n, std::uint64_t { 1 }); }
    std::size_t numberOfWorkers() const { return m_nThreads; }
    void setLivermoreComptonModel(bool use) { m_useComptonLivermore = use; }
    bool livermoreComptonModel() const { return m_useComptonLivermore; }
    void setBindingEnergyCorrection(bool use) { m_useBindingEnergyCorrection = use; }
    bool bindingEnergyCorrection() const { return m_useBindingEnergyCorrection; }
    void setSiddonTracking(bool use) { m_useSiddonTracking = use; }
    bool siddonTracking() const { return m_useSiddonTracking; }

    template <typename U>
    requires std::is_base_of<World<T>, U>::value
        Result<T>
        operator()(const U& world, Source<T>* source, ProgressBar<T>* progressbar = nullptr, bool calculateDose = true)
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
        constexpr T minEnergy { 1 };
        const auto energyStep = std::min(maxEnergy / 100, T { 1 });
        m_attenuationLut.generate(world, minEnergy, maxEnergy, energyStep);

        const std::uint64_t totalExposures = source->totalExposures();

        const std::uint64_t nJobs = std::max(m_nThreads, std::uint64_t { 1 });
        if (progressbar) {
            progressbar->setTotalExposures(totalExposures);
            progressbar->setDoseData(result.dose.data(), world.dimensions(), world.spacing());
        }
        const auto start = std::chrono::system_clock::now();
        if (m_useComptonLivermore) {
            if (m_useBindingEnergyCorrection) {
                if (m_useSiddonTracking)
                    parallellRun<true, true, true>(world, source, result, 0, totalExposures, nJobs, progressbar);
                else
                    parallellRun<true, true, false>(world, source, result, 0, totalExposures, nJobs, progressbar);
            } else {
                if (m_useSiddonTracking)
                    parallellRun<true, false, true>(world, source, result, 0, totalExposures, nJobs, progressbar);
                else
                    parallellRun<true, false, false>(world, source, result, 0, totalExposures, nJobs, progressbar);
            }
        } else {
            if (m_useBindingEnergyCorrection) {
                if (m_useSiddonTracking)
                    parallellRun<false, true, true>(world, source, result, 0, totalExposures, nJobs, progressbar);
                else
                    parallellRun<false, true, false>(world, source, result, 0, totalExposures, nJobs, progressbar);
            } else {
                if (m_useSiddonTracking)
                    parallellRun<false, false, true>(world, source, result, 0, totalExposures, nJobs, progressbar);
                else
                    parallellRun<false, false, false>(world, source, result, 0, totalExposures, nJobs, progressbar);
            }
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

        //compute variance
        computeResultVariance(result);

        if (calculateDose) {
            const T calibrationValue = source->getCalibrationValue(progressbar);
            //energy imparted to dose
            energyImpartedToDose(world, result, calibrationValue);
            result.dose_units = "mGy";
        } else {
            energyImpartedToDose(world, result);
            result.dose_units = "keV/kg";
        }
        return result;
    }

    const AttenuationLut<T>& attenuationLut() const { return m_attenuationLut; }
    AttenuationLut<T>& attenuationLut() { return m_attenuationLut; }

protected:
    template <typename U>
    inline void safeValueAdd(U& value, const U addValue) const noexcept
    {
        std::atomic_ref value_l(value);
        value_l.fetch_add(addValue);
    }

    void rayleightScatterLivermore(Particle<T>& particle, std::uint8_t materialIdx, RandomState& state) const noexcept
    {
        // theta is scattering angle
        // see http://rcwww.kek.jp/research/egs/egs5_manual/slac730-150228.pdf

        //finding qmax
        const auto qmax = m_attenuationLut.momentumTransferMax(particle.energy);
        const auto amax = m_attenuationLut.cumFormFactorSquared(materialIdx, qmax);
        T cosAngle;
        do {
            const auto r1 = state.randomUniform<T>();
            const auto aatq = amax * r1;
            const auto q = m_attenuationLut.momentumTransfer(static_cast<size_t>(materialIdx), aatq);
            cosAngle = m_attenuationLut.cosAngle(particle.energy, q);
        } while ((T { 0.5 } + cosAngle * cosAngle * T { 0.5 }) < state.randomUniform<T>());

        const auto theta = std::acos(cosAngle);
        const auto phi = state.randomUniform<T>(PI_VAL<T>() * PI_VAL<T>());
        vectormath::peturb<T>(particle.dir.data(), theta, phi);
    }

    template <bool Livermore = true>
    T comptonScatter(Particle<T>& particle, std::uint8_t materialIdx, RandomState& state) const noexcept
    // see http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/PhysicsReferenceManual/fo/PhysicsReferenceManual.pdf
    // and
    // https://nrc-cnrc.github.io/EGSnrc/doc/pirs701-egsnrc.pdf
    {
        const auto k = particle.energy / ELECTRON_REST_MASS<T>();
        const auto emin = 1 / (1 + 2 * k);
        const auto gmax_inv = 1 / (1 / emin + emin);
        T e, cosAngle;
        bool rejected;
        do {
            const auto r1 = state.randomUniform<T>();
            e = r1 + (1 - r1) * emin;
            const T t = (1 - e) / (k * e);
            const auto sinthetasqr = t * (2 - t);

            cosAngle = 1 - t;

            const auto g = (1 / e + e - sinthetasqr) * gmax_inv;
            if constexpr (Livermore) {
                const auto q = m_attenuationLut.momentumTransferFromCos(particle.energy, cosAngle);
                const auto scatterFactor = m_attenuationLut.comptonScatterFactor(static_cast<std::size_t>(materialIdx), q);
                const auto r2 = state.randomUniform<T>();
                rejected = r2 > g * scatterFactor;
            } else {
                const auto r2 = state.randomUniform<T>();
                rejected = r2 > g;
            }

        } while (rejected);
        const auto theta = std::acos(cosAngle);
        const auto phi = state.randomUniform<T>(PI_VAL<T>() + PI_VAL<T>());
        vectormath::peturb<T>(particle.dir.data(), theta, phi);
        const auto E0 = particle.energy;
        particle.energy *= e;

        return E0 * (T { 1 } - e);
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

    template <bool Livermore, bool bindingEnergyCorrection>
    bool computeInteractionsForced(const T eventProbability, Particle<T>& p, const std::uint8_t matIdx, Result<T>& result, std::size_t resultBufferIdx, RandomState& state, bool& updateMaxAttenuation) const noexcept
    {
        const auto atts = m_attenuationLut.photoComptRayAttenuation(matIdx, p.energy);
        const auto attPhoto = atts[0];
        const auto attCompt = atts[1];
        const auto attRayl = atts[2];
        const auto attenuationTotal = attPhoto + attCompt + attRayl;

        const auto weightCorrection = eventProbability * attPhoto / attenuationTotal;
        if constexpr (bindingEnergyCorrection) {
            const auto bindingEnergy = m_attenuationLut.meanBindingEnergy(matIdx);
            const auto energyImparted = (p.energy - bindingEnergy) * p.weight * weightCorrection;

            safeValueAdd(result.dose[resultBufferIdx], energyImparted);
            safeValueAdd(result.nEvents[resultBufferIdx], std::uint32_t { 1 });
            safeValueAdd(result.variance[resultBufferIdx], energyImparted * energyImparted);

        } else {
            const auto energyImparted = p.energy * p.weight * weightCorrection;
            safeValueAdd(result.dose[resultBufferIdx], energyImparted);
            safeValueAdd(result.nEvents[resultBufferIdx], std::uint32_t { 1 });
            safeValueAdd(result.variance[resultBufferIdx], energyImparted * energyImparted);
        }
        p.weight *= (1 - weightCorrection); // to prevent bias

        const auto r2 = state.randomUniform<T>();
        if (r2 < eventProbability) {

            const auto r3 = state.randomUniform(attCompt + attRayl);

            if (r3 <= attCompt) // Compton event
            {
                const auto e = comptonScatter<Livermore>(p, matIdx, state);
                if (p.energy < ENERGY_CUTOFF_THRESHOLD())
                    [[unlikely]]
                    {
                        if constexpr (bindingEnergyCorrection) {
                            const auto bindingEnergy = m_attenuationLut.meanBindingEnergy(matIdx);
                            const auto energyImparted = (e + p.energy - bindingEnergy) * p.weight;

                            safeValueAdd(result.dose[resultBufferIdx], energyImparted);
                            safeValueAdd(result.nEvents[resultBufferIdx], std::uint32_t { 1 });
                            safeValueAdd(result.variance[resultBufferIdx], energyImparted * energyImparted);

                        } else {
                            const auto energyImparted = (e + p.energy) * p.weight;
                            safeValueAdd(result.dose[resultBufferIdx], energyImparted);
                            safeValueAdd(result.nEvents[resultBufferIdx], std::uint32_t { 1 });
                            safeValueAdd(result.variance[resultBufferIdx], energyImparted * energyImparted);
                        }
                        return false;
                    }
                else
                    [[likely]]
                    {
                        if constexpr (bindingEnergyCorrection) {
                            const auto bindingEnergy = m_attenuationLut.meanBindingEnergy(matIdx);
                            const auto energyImparted = (e - bindingEnergy) * p.weight;
                            safeValueAdd(result.dose[resultBufferIdx], energyImparted);
                            safeValueAdd(result.nEvents[resultBufferIdx], std::uint32_t { 1 });
                            safeValueAdd(result.variance[resultBufferIdx], energyImparted * energyImparted);

                        } else {
                            const auto energyImparted = e * p.weight;
                            safeValueAdd(result.dose[resultBufferIdx], energyImparted);
                            safeValueAdd(result.nEvents[resultBufferIdx], std::uint32_t { 1 });
                            safeValueAdd(result.variance[resultBufferIdx], energyImparted * energyImparted);
                        }
                        updateMaxAttenuation = true;
                    }
            } else // Rayleigh scatter event
            {
                rayleightScatterLivermore(p, matIdx, state);
            }
        }
        return true;
    }
    template <bool Livermore, bool bindingEnergyCorrection>
    bool computeInteractions(Particle<T>& p, const std::uint8_t matIdx, Result<T>& result, std::size_t resultBufferIdx, RandomState& state, bool& updateMaxAttenuation) const noexcept
    {
        const auto atts = m_attenuationLut.photoComptRayAttenuation(matIdx, p.energy);
        const auto attPhoto = atts[0];
        const auto attCompt = atts[1];
        const auto attRayl = atts[2];

        const auto r3 = state.randomUniform(attPhoto + attCompt + attRayl);
        if (r3 <= attPhoto) // Photoelectric event
        {
            if constexpr (bindingEnergyCorrection) {
                const auto bindingEnergy = m_attenuationLut.meanBindingEnergy(matIdx);
                const auto energyImparted = (p.energy - bindingEnergy) * p.weight;

                safeValueAdd(result.dose[resultBufferIdx], energyImparted);
                safeValueAdd(result.nEvents[resultBufferIdx], std::uint32_t { 1 });
                safeValueAdd(result.variance[resultBufferIdx], energyImparted * energyImparted);

            } else {
                const auto energyImparted = p.energy * p.weight;
                safeValueAdd(result.dose[resultBufferIdx], energyImparted);
                safeValueAdd(result.nEvents[resultBufferIdx], std::uint32_t { 1 });
                safeValueAdd(result.variance[resultBufferIdx], energyImparted * energyImparted);
            }
            p.energy = 0.0;
            return false;
        } else if (r3 <= attPhoto + attCompt) // Compton event
        {
            const auto e = comptonScatter<Livermore>(p, matIdx, state);
            if (p.energy < ENERGY_CUTOFF_THRESHOLD())
                [[unlikely]]
                {
                    if constexpr (bindingEnergyCorrection) {
                        const auto bindingEnergy = m_attenuationLut.meanBindingEnergy(matIdx);
                        const auto energyImparted = (e + p.energy - bindingEnergy) * p.weight;

                        safeValueAdd(result.dose[resultBufferIdx], energyImparted);
                        safeValueAdd(result.nEvents[resultBufferIdx], std::uint32_t { 1 });
                        safeValueAdd(result.variance[resultBufferIdx], energyImparted * energyImparted);

                    } else {
                        const auto energyImparted = (e + p.energy) * p.weight;
                        safeValueAdd(result.dose[resultBufferIdx], energyImparted);
                        safeValueAdd(result.nEvents[resultBufferIdx], std::uint32_t { 1 });
                        safeValueAdd(result.variance[resultBufferIdx], energyImparted * energyImparted);
                    }
                    return false;
                }
            else
                [[likely]]
                {
                    if constexpr (bindingEnergyCorrection) {
                        const auto bindingEnergy = m_attenuationLut.meanBindingEnergy(matIdx);
                        const auto energyImparted = (e - bindingEnergy) * p.weight;

                        safeValueAdd(result.dose[resultBufferIdx], energyImparted);
                        safeValueAdd(result.nEvents[resultBufferIdx], std::uint32_t { 1 });
                        safeValueAdd(result.variance[resultBufferIdx], energyImparted * energyImparted);

                    } else {
                        const auto energyImparted = e * p.weight;
                        safeValueAdd(result.dose[resultBufferIdx], energyImparted);
                        safeValueAdd(result.nEvents[resultBufferIdx], std::uint32_t { 1 });
                        safeValueAdd(result.variance[resultBufferIdx], energyImparted * energyImparted);
                    }
                    updateMaxAttenuation = true;
                }
        } else // Rayleigh scatter event
        {
            rayleightScatterLivermore(p, matIdx, state);
        }

        return true;
    }

    template <bool Livermore, bool bindingEnergyCorrection>
    void siddonParticleTracking(const World<T>& world, Particle<T>& p, RandomState& state, Result<T>& result) const noexcept
    {
        const auto densBuffer = world.densityArray()->data();
        const auto matBuffer = world.materialIndexArray()->data();
        const auto measBuffer = world.measurementMapArray()->data();

        //parameterized path: pos1 = pos0 + alpha*dir

        const auto& extent = world.matrixExtent();
        const auto& spacing = world.spacing();
        const auto& dim = world.dimensions();
        const std::array<T, 3> offset = { extent[0], extent[2], extent[4] };
        auto index = gridIndexFromPosition(p.pos, world);

        std::array<int, 3> direction = {
            p.dir[0] < 0 ? -1 : 1,
            p.dir[1] < 0 ? -1 : 1,
            p.dir[2] < 0 ? -1 : 1
        };

        std::array<T, 3> alpha = { 0, 0, 0 };
        for (std::size_t i = 0; i < 3; ++i) {
            if (std::abs(p.dir[i]) > N_ERROR()) {
                const auto nextBorder = direction[i] < 0 ? index[i] : index[i] + 1;
                alpha[i] = (offset[i] + nextBorder * spacing[i] - p.pos[i]) / p.dir[i];
            } else {
                alpha[i] = std::numeric_limits<T>::max();
            }
        }

        auto alpha_ind = vectormath::argmin3<std::size_t>(alpha.data());

        auto continueSampling = true;
        auto energyChanged = false;
        do {
            const auto arrIdx = index[0] + index[1] * dim[0] + index[2] * dim[0] * dim[1];
            const auto matIdx = matBuffer[arrIdx];
            const auto attenuation = m_attenuationLut.totalAttenuation(matIdx, p.energy) * densBuffer[arrIdx] * T { 0.1 }; // cm->mm
            const auto interactionProb = std::exp(-alpha[alpha_ind] * attenuation);
            /* Denne virker ikke som den skal, to photo events kan potensielt skje.
            if (measBuffer[arrIdx]) {
                const auto attPhoto = m_attenuationLut.photoelectricAttenuation(matIdx, p.energy);
                const auto weightCorrection = (1 - interactionProb) * attPhoto * densBuffer[arrIdx] * T { 0.1 } / attenuation;
                if constexpr (bindingEnergyCorrection) {
                    const auto bindingEnergy = m_attenuationLut.meanBindingEnergy(matIdx);
                    const auto energyImparted = (p.energy - bindingEnergy) * p.weight * weightCorrection;

                    safeValueAdd(result.dose[arrIdx], energyImparted);
                    safeValueAdd(result.nEvents[arrIdx], std::uint32_t { 1 });
                    safeValueAdd(result.variance[arrIdx], energyImparted * energyImparted);

                } else {
                    const auto energyImparted = p.energy * p.weight * weightCorrection;
                    safeValueAdd(result.dose[arrIdx], energyImparted);
                    safeValueAdd(result.nEvents[arrIdx], std::uint32_t { 1 });
                    safeValueAdd(result.variance[arrIdx], energyImparted * energyImparted);
                }
                p.weight *= (T { 1 } - weightCorrection); // to prevent bias
            }
            */

            const auto r1 = state.randomUniform<T>();
            if (r1 > interactionProb) //event
            {
                const auto step = -std::log(r1) / attenuation;
                for (std::size_t i = 0; i < 3; ++i) {
                    p.pos[i] += p.dir[i] * step;
                }

                continueSampling = computeInteractions<Livermore, bindingEnergyCorrection>(p, matBuffer[arrIdx], result, arrIdx, state, energyChanged);
                if (continueSampling) {

                    //new steplenght
                    direction = {
                        p.dir[0] < 0 ? -1 : 1,
                        p.dir[1] < 0 ? -1 : 1,
                        p.dir[2] < 0 ? -1 : 1
                    };

                    for (std::size_t i = 0; i < 3; ++i) {
                        if (std::abs(p.dir[i]) > N_ERROR()) {
                            const auto nextBorder = direction[i] < 0 ? index[i] : index[i] + 1;
                            alpha[i] = (offset[i] + nextBorder * spacing[i] - p.pos[i]) / p.dir[i];
                        } else {
                            alpha[i] = std::numeric_limits<T>::max();
                        }
                    }
                    alpha_ind = vectormath::argmin3<std::size_t>(alpha.data());
                }
            } else {
                for (std::size_t i = 0; i < 3; ++i) {
                    p.pos[i] += p.dir[i] * alpha[alpha_ind];
                }
                for (std::size_t i = 0; i < 3; ++i) {
                    alpha[i] -= alpha[alpha_ind];
                }

                //update next voxel
                index[alpha_ind] += direction[alpha_ind];
                if (!(index[alpha_ind] < dim[alpha_ind])) {
                    //we have reach the end of the world, note index is unsigned type and will overflow
                    return;
                }

                //update next alpha
                const auto nextBorder = direction[alpha_ind] < 0 ? index[alpha_ind] : index[alpha_ind] + 1;
                alpha[alpha_ind] = (offset[alpha_ind] + nextBorder * spacing[alpha_ind] - p.pos[alpha_ind]) / p.dir[alpha_ind];
                alpha_ind = vectormath::argmin3<std::size_t>(alpha.data());
            }
        } while (continueSampling);
    }

    template <bool Livermore, bool bindingEnergyCorrection>
    void woodcockParticleTracking(const World<T>& world, Particle<T>& p, RandomState& state, Result<T>& result) const noexcept
    {
        const auto densityBuffer = world.densityArray()->data();
        const auto materialBuffer = world.materialIndexArray()->data();
        const auto measurementBuffer = world.measurementMapArray()->data();

        T maxAttenuationInv;
        bool updateMaxAttenuation = true;
        bool continueSampling = true;
        bool ruletteCandidate = true;
        while (continueSampling) {
            if (updateMaxAttenuation) {
                maxAttenuationInv = T { 1 } / m_attenuationLut.maxMassTotalAttenuation(p.energy);
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
                const auto attenuationTotal = m_attenuationLut.totalAttenuation(matIdx, p.energy) * density;
                const auto eventProbability = attenuationTotal * maxAttenuationInv;

                if (measurementBuffer[bufferIdx] == 0)
                    [[likely]] // naive sampling
                    {
                        const auto r2 = state.randomUniform<T>();
                        if (r2 < eventProbability) // an event will happen
                        {
                            continueSampling = computeInteractions<Livermore, bindingEnergyCorrection>(p, matIdx, result, bufferIdx, state, updateMaxAttenuation);
                        }
                    }
                else
                    [[unlikely]] // forced photoelectric effect
                    {
                        continueSampling = computeInteractionsForced<Livermore, bindingEnergyCorrection>(eventProbability, p, matIdx, result, bufferIdx, state, updateMaxAttenuation);
                    }

                if (continueSampling) {
                    if ((p.energy < RUSSIAN_RULETTE_ENERGY_THRESHOLD()) && ruletteCandidate) {
                        ruletteCandidate = false;
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
    template <bool Livermore, bool bindingEnergyCorrection, bool siddonTracking>
    void transport(const World<T>& world, const Exposure<T>& exposure, RandomState& state, Result<T>& result) const noexcept
    {
        const auto nHistories = exposure.numberOfHistories();
        for (std::size_t i = 0; i < nHistories; ++i) {
            //Draw a particle
            auto particle = exposure.sampleParticle(state);

            // Is particle intersecting with world
            const bool isInWorld = transportParticleToWorld(world, particle);
            if (isInWorld) {
                if constexpr (siddonTracking) {
                    siddonParticleTracking<Livermore, bindingEnergyCorrection>(world, particle, state, result);

                } else {
                    woodcockParticleTracking<Livermore, bindingEnergyCorrection>(world, particle, state, result);
                }
            }
        }
    }

    template <bool Livermore, bool bindingEnergyCorrection, bool siddonTracking>
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
            transport<Livermore, bindingEnergyCorrection, siddonTracking>(w, exposure, state, result);
            if (progressbar)
                progressbar->exposureCompleted();
            job = taskpool.getJob();
        }
    }

    template <bool Livermore, bool bindingEnergyCorrection, bool siddonTracking = false>
    void parallellRun(const World<T>& w, const Source<T>* source, Result<T>& result,
        const std::uint64_t expBeg, const std::uint64_t expEnd, std::uint64_t nJobs, ProgressBar<T>* progressbar) const
    {
        EvenTaskPool taskpool(expBeg, expEnd, nJobs);
        std::vector<std::thread> jobs;
        jobs.reserve(nJobs - 1);
        for (std::size_t i = 0; i < nJobs - 1; ++i) {
            jobs.emplace_back(&Transport<T>::workerRun<Livermore, bindingEnergyCorrection, siddonTracking>, this, std::cref(w), source, std::ref(result), std::ref(taskpool), progressbar);
        }
        workerRun<Livermore, bindingEnergyCorrection, siddonTracking>(w, source, result, taskpool, progressbar);
        for (auto& job : jobs)
            job.join();
    }

    void computeResultVariance(Result<T>& res) const noexcept
    {
        //Computing variance by: Var[X] = E[X**2] - E[X]**2
        //and Var[A]+ Var[A] = 2Var[A]

        /*auto eBeg = res.dose.cbegin();
        auto eEnd = res.dose.cend();
        auto tBeg = res.nEvents.cbegin();
        auto vBeg = res.variance.begin();
        while (eBeg != eEnd) {
            if (*tBeg > 0) {
                const auto nEv = *tBeg;
                *vBeg = *vBeg / nEv - (*eBeg) * (*eBeg) / nEv / nEv;
            }
            ++eBeg;
            ++tBeg;
            ++vBeg;
        }*/
        const auto nh = res.numberOfHistories;
        std::transform(std::execution::par_unseq, res.dose.cbegin(), res.dose.cend(), res.variance.cbegin(), res.variance.begin(),
            [=](const auto d, const auto dd) -> T {
                //const auto sd = std::sqrt((dd / nh - (d / nh) * (d / nh)) / nh); // uncertenty per event
                const auto sd = std::sqrt(dd - (d * d / nh)); // uncertenty for total dose in voxel
                return sd / d * 100;
            });
    }

    void energyImpartedToDose(const World<T>& world, Result<T>& res, const T calibrationValue = 1) noexcept
    {
        const auto& spacing = world.spacing();
        const auto voxelVolume = spacing[0] * spacing[1] * spacing[2] / T { 1000.0 }; // cm3
        auto density = world.densityArray();

        std::transform(
            std::execution::par_unseq, res.dose.cbegin(), res.dose.cend(), density->cbegin(), res.dose.begin(),
            [=](auto ei, auto de) -> auto {
                const auto voxelMass = de * voxelVolume * T { 0.001 }; //kg
                return de > T { 0.0 } ? calibrationValue * ei / voxelMass : T { 0.0 };
            });
        std::transform(
            std::execution::par_unseq, res.variance.cbegin(), res.variance.cend(), density->cbegin(), res.variance.begin(),
            [=](auto var, auto de) -> auto {
                const auto voxelMass = de * voxelVolume * T { 0.001 }; //kg
                const auto factor = calibrationValue * calibrationValue / (voxelMass * voxelMass);
                return de > T { 0.0 } ? factor * var : T { 0.0 };
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
    bool m_useComptonLivermore = true;
    bool m_useBindingEnergyCorrection = true;
    bool m_useSiddonTracking = false;
    std::uint64_t m_nThreads;
};
}