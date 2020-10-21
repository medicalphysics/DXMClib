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
#include <memory>

namespace dxmc {
/**
     * @brief A simple holder for atomic locks while we wait for atomic_ref support in all major compilers
    */
struct resultLock {
    std::atomic_flag dose;
    std::atomic_flag nEvents;
    std::atomic_flag variance;
    resultLock()
    {
    }

    resultLock(const resultLock& other)
    {
    }

    resultLock& operator=(const resultLock& other)
    {
    }
};

template <Floating T>
struct Result {
    std::vector<T> dose;
    std::vector<std::uint32_t> nEvents;
    std::vector<T> variance;
    std::vector<resultLock> locks;

    std::chrono::duration<float> simulationTime { 0 };

    Result(std::size_t size)
    {
        dose.resize(size, 0.0);
        nEvents.resize(size, 0);
        variance.resize(size, 0.0);
        locks.resize(size);
    }
};

// forward declaration
template <Floating T>
class Source;
template <Floating T>
class CTSource;

template <Floating T>
class Transport {
public:
    Transport() { }
    Result<T> operator()(const World<T>& world, Source<T>* source, ProgressBar<T>* progressbar = nullptr, bool calculateDose = true);

    Result<T> operator()(const CTDIPhantom<T>& world, CTSource<T>* source, ProgressBar<T>* progressbar = nullptr);
    const AttenuationLut<T>& attenuationLut() const { return m_attenuationLut; }
    AttenuationLut<T>& attenuationLut() { return m_attenuationLut; }

protected:
    template <typename U>
    void safeValueAdd(U& value, const U addValue, std::atomic_flag& lock) const
    {
        while (lock.test_and_set(std::memory_order_acquire)) // acquire lock
            while (lock.test(std::memory_order_relaxed)) // test lock
                ; // spin
        value += addValue;
        lock.clear(std::memory_order_release); // release lock
    }

    void rayleightScatterLivermore(Particle<T>& particle, unsigned char materialIdx, RandomState& state, T& cosAngle) const
    {
        // theta is scattering angle
        // see http://rcwww.kek.jp/research/egs/egs5_manual/slac730-150228.pdf

        //finding qmax
        const auto qmax = m_attenuationLut.momentumTransferMax(particle.energy);
        const auto amax = m_attenuationLut.cumFormFactorSquared(materialIdx, qmax);

        do {
            const auto r1 = state.randomUniform<T>();
            const auto aatq = amax * r1;
            const auto q = m_attenuationLut.momentumTransfer(static_cast<size_t>(materialIdx), aatq);
            cosAngle = m_attenuationLut.cosAngle(particle.energy, q);
        } while ((T { 0.5 } + cosAngle * cosAngle * T { 0.5 }) < state.randomUniform<T>());

        const auto theta = std::acos(cosAngle);
        const auto phi = state.randomUniform<T>(PI_VAL<T>() * PI_VAL<T>());
        vectormath::peturb<T>(particle.dir, theta, phi);
    }
    T comptonScatterLivermore(Particle<T>& particle, unsigned char materialIdx, RandomState& state, T& cosAngle) const
    // see http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/PhysicsReferenceManual/fo/PhysicsReferenceManual.pdf
    // and
    // https://nrc-cnrc.github.io/EGSnrc/doc/pirs701-egsnrc.pdf
    {
        const auto k = particle.energy / ELECTRON_REST_MASS<T>();
        const auto emin = T { 1.0 } / (T { 1.0 } + T { 2.0 } * k);
        const auto gmax = T { 1.0 } / emin + emin;
        T e;
        bool rejected;
        do {
            const auto r1 = state.randomUniform<T>();
            e = r1 + (T { 1.0 } - r1) * emin;
            cosAngle = T { 1.0 } + T { 1.0 } / k - T { 1.0 } / (k * e);
            const auto sinthetasqr = T { 1.0 } - cosAngle * cosAngle;

            const auto g = (T { 1.0 } / e + e - sinthetasqr) / gmax;

            const auto q = m_attenuationLut.momentumTransferFromCos(particle.energy, cosAngle);
            const auto scatterFactor = m_attenuationLut.comptonScatterFactor(static_cast<std::size_t>(materialIdx), q);
            const auto r2 = state.randomUniform<T>();
            rejected = r2 > g * scatterFactor;

        } while (rejected);

        const auto theta = std::acos(cosAngle);
        const auto phi = state.randomUniform<T>(PI_VAL<T>() + PI_VAL<T>());
        vectormath::peturb<T>(particle.dir, theta, phi);

        const auto E0 = particle.energy;
        particle.energy *= e;

        return E0 * (T { 1.0 } - e);
    }

    T comptonScatter(Particle<T>& particle, RandomState& state, T& cosAngle) const
    // see http://geant4-userdoc.web.cern.ch/geant4-userdoc/UsersGuides/PhysicsReferenceManual/fo/PhysicsReferenceManual.pdf
    // and
    // https://nrc-cnrc.github.io/EGSnrc/doc/pirs701-egsnrc.pdf
    {
        const auto k = particle.energy / ELECTRON_REST_MASS<T>();
        const auto emin = T { 1.0 } / (T { 1.0 } + T { 2.0 } * k);
        const auto gmax = T { 1.0 } / emin + emin;

        bool rejected;
        T sinthetasqr, e, t;
        do {
            const auto r1 = state.randomUniform<T>();
            e = r1 + (T { 1.0 } - r1) * emin;

            t = (T { 1.0 } - e) / (k * e);
            sinthetasqr = t * (T { 2.0 } - t);

            const auto g = (T { 1.0 } / e + e - sinthetasqr) / gmax;

            const auto r2 = state.randomUniform<T>();
            rejected = r2 > g;
        } while (rejected);

        cosAngle = T { 1.0 } - t;
        const auto theta = std::acos(cosAngle);
        const auto phi = state.randomUniform<T>(PI_VAL<T>() + PI_VAL<T>());
        vectormath::peturb<T>(particle.dir, theta, phi);

        const auto E0 = particle.energy;
        particle.energy *= e;

        return E0 * (T { 1.0 } - e);
    }

    inline bool particleInsideWorld(const World<T>& world, const Particle<T>& particle) const
    {
        auto extent = world.matrixExtent();
        return (particle.pos[0] > extent[0] && particle.pos[0] < extent[1]) && (particle.pos[1] > extent[2] && particle.pos[1] < extent[3]) && (particle.pos[2] > extent[4] && particle.pos[2] < extent[5]);
    }

    inline std::size_t indexFromPosition(const Particle<T>& particle, const World<T>& world) const
    {
        //assumes particle is inside world
        std::size_t arraypos[3];
        const auto& wpos = world.matrixExtent();
        const auto& wdim = world.dimensions();
        const auto& wspac = world.spacing();

        for (std::size_t i = 0; i < 3; i++)
            arraypos[i] = static_cast<std::size_t>((particle.pos[i] - wpos[i * 2]) / wspac[i]);
        std::size_t idx = arraypos[2] * wdim[0] * wdim[1] + arraypos[1] * wdim[0] + arraypos[0];
        return idx;
    }
    bool computeInteractionsForced(const T eventProbability, Particle<T>& p, const unsigned char matIdx, Result<T>* result, std::size_t resultBufferIdx, RandomState& state, bool& updateMaxAttenuation) const
    {
        const auto atts = m_attenuationLut.photoComptRayAttenuation(matIdx, p.energy);
        const auto attPhoto = atts[0];
        const auto attCompt = atts[1];
        const auto attRayl = atts[2];
        const auto attenuationTotal = attPhoto + attCompt + attRayl;

        const auto weightCorrection = eventProbability * attPhoto / attenuationTotal;
        const auto energyImparted = p.energy * p.weight * weightCorrection;
        safeValueAdd(result->dose[resultBufferIdx], energyImparted, result->locks[resultBufferIdx].dose);
        safeValueAdd(result->nEvents[resultBufferIdx], std::uint32_t { 1 }, result->locks[resultBufferIdx].nEvents);
        safeValueAdd(result->variance[resultBufferIdx], energyImparted * energyImparted, result->locks[resultBufferIdx].variance);
        p.weight = p.weight * (1.0 - weightCorrection); // to prevent bias

        const auto r2 = state.randomUniform<T>();
        if (r2 < eventProbability) {

            const auto r3 = state.randomUniform(attCompt + attRayl);

            if (r3 <= attCompt) // Compton event
            {
                T cosangle;
#ifdef DXMC_USE_LOWENERGY_COMPTON
                const auto e = comptonScatterLivermore(p, matIdx, state, cosangle);
#else
                const auto e = comptonScatter(p, state, cosangle);
#endif
                if (p.energy < ENERGY_CUTOFF_THRESHOLD())
                    [[unlikely]]
                    {
                        const auto energyImparted = (e + p.energy) * p.weight;
                        safeValueAdd(result->dose[resultBufferIdx], energyImparted, result->locks[resultBufferIdx].dose);
                        safeValueAdd(result->nEvents[resultBufferIdx], std::uint32_t { 1 }, result->locks[resultBufferIdx].nEvents);
                        safeValueAdd(result->variance[resultBufferIdx], energyImparted * energyImparted, result->locks[resultBufferIdx].variance);
                        return false;
                    }
                else
                    [[likely]]
                    {
                        const auto energyImparted = e * p.weight;
                        safeValueAdd(result->dose[resultBufferIdx], energyImparted, result->locks[resultBufferIdx].dose);
                        safeValueAdd(result->nEvents[resultBufferIdx], std::uint32_t { 1 }, result->locks[resultBufferIdx].nEvents);
                        safeValueAdd(result->variance[resultBufferIdx], energyImparted * energyImparted, result->locks[resultBufferIdx].variance);
                        updateMaxAttenuation = true;
                    }
            } else // Rayleigh scatter event
            {
                T cosangle;
                rayleightScatterLivermore(p, matIdx, state, cosangle);
            }
        }
        return true;
    }
    bool computeInteractions(Particle<T>& p, const unsigned char matIdx, Result<T>* result, std::size_t resultBufferIdx, RandomState& state, bool& updateMaxAttenuation) const
    {
        const auto atts = m_attenuationLut.photoComptRayAttenuation(matIdx, p.energy);
        const auto attPhoto = atts[0];
        const auto attCompt = atts[1];
        const auto attRayl = atts[2];

        const auto r3 = state.randomUniform(attPhoto + attCompt + attRayl);
        if (r3 <= attPhoto) // Photoelectric event
        {
            const auto energyImparted = p.energy * p.weight;
            safeValueAdd(result->dose[resultBufferIdx], energyImparted, result->locks[resultBufferIdx].dose);
            safeValueAdd(result->nEvents[resultBufferIdx], std::uint32_t { 1 }, result->locks[resultBufferIdx].nEvents);
            safeValueAdd(result->variance[resultBufferIdx], energyImparted * energyImparted, result->locks[resultBufferIdx].variance);
            p.energy = 0.0;
            return false;
        } else if (r3 <= attPhoto + attCompt) // Compton event
        {
            T cosangle;
#ifdef DXMC_USE_LOWENERGY_COMPTON
            const auto e = comptonScatterLivermore(p, matIdx, state, cosangle);
#else
            const auto e = comptonScatter(p, state, cosangle);
#endif
            if (p.energy < ENERGY_CUTOFF_THRESHOLD())
                [[unlikely]]
                {
                    const auto energyImparted = (e + p.energy) * p.weight;
                    safeValueAdd(result->dose[resultBufferIdx], energyImparted, result->locks[resultBufferIdx].dose);
                    safeValueAdd(result->nEvents[resultBufferIdx], std::uint32_t { 1 }, result->locks[resultBufferIdx].nEvents);
                    safeValueAdd(result->variance[resultBufferIdx], energyImparted * energyImparted, result->locks[resultBufferIdx].variance);
                    return false;
                }
            else
                [[likely]]
                {
                    const auto energyImparted = e * p.weight;
                    safeValueAdd(result->dose[resultBufferIdx], energyImparted, result->locks[resultBufferIdx].dose);
                    safeValueAdd(result->nEvents[resultBufferIdx], std::uint32_t { 1 }, result->locks[resultBufferIdx].nEvents);
                    safeValueAdd(result->variance[resultBufferIdx], energyImparted * energyImparted, result->locks[resultBufferIdx].variance);
                    updateMaxAttenuation = true;
                }
        } else // Rayleigh scatter event
        {
            T cosangle;
            rayleightScatterLivermore(p, matIdx, state, cosangle);
        }
        return true;
    }

    void sampleParticleSteps(const World<T>& world, Particle<T>& p, RandomState& state, Result<T>* result) const
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
                maxAttenuationInv = T { 1.0 } / m_attenuationLut.maxMassTotalAttenuation(p.energy);
                updateMaxAttenuation = false;
            }
            const auto r1 = state.randomUniform<T>();
            const auto stepLenght = -std::log(r1) * maxAttenuationInv * T { 10.0 }; // cm -> mm
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
                            continueSampling = computeInteractions(p, matIdx, result, bufferIdx, state, updateMaxAttenuation);
                        }
                    }
                else
                    [[unlikely]] // forced photoelectric effect
                    {
                        continueSampling = computeInteractionsForced(eventProbability, p, matIdx, result, bufferIdx, state, updateMaxAttenuation);
                    }

                if (continueSampling) {
                    if ((p.energy < RUSSIAN_RULETTE_ENERGY_THRESHOLD()) && ruletteCandidate) {
                        ruletteCandidate = false;
                        const auto r4 = state.randomUniform<T>();
                        if (r4 < RUSSIAN_RULETTE_PROBABILITY()) {
                            continueSampling = false;
                        } else {
                            constexpr T factor = T { 1.0 } / (T { 1.0 } - RUSSIAN_RULETTE_PROBABILITY());
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

    bool transportParticleToWorld(const World<T>& world, Particle<T>& particle) const
    //returns false if particle do not intersect world
    {
        const bool isInside = particleInsideWorld(world, particle);
        if (isInside)
            return true;

        auto amin = std::numeric_limits<T>::min();
        auto amax = std::numeric_limits<T>::max();

        const auto& extent = world.matrixExtent();
        for (std::size_t i = 0; i < 3; i++) {
            if (std::abs(particle.dir[i]) > N_ERROR()) {
                T a0, an;
                a0 = (extent[i * 2] - particle.pos[i]) / particle.dir[i];
                an = (extent[i * 2 + 1] - particle.pos[i]) / particle.dir[i];
                amin = std::max(amin, std::min(a0, an));
                amax = std::min(amax, std::max(a0, an));
            }
        }
        if (amin < amax && amin > 0.0) {
            for (std::size_t i = 0; i < 3; i++) {
                particle.pos[i] += (amin + N_ERROR()) * particle.dir[i]; // making sure particle is inside world
            }
            return true;
        }
        return false;
    }

    void transport(const World<T>& world, const Exposure<T>& exposure, RandomState& state, Result<T>* result) const
    {
        Particle<T> particle;
        const auto nHistories = exposure.numberOfHistories();
        for (std::size_t i = 0; i < nHistories; ++i) {
            //Draw a particle
            exposure.sampleParticle(particle, state);

            // Is particle intersecting with world
            const bool isInWorld = transportParticleToWorld(world, particle);
            if (isInWorld)
                sampleParticleSteps(world, particle, state, result);
        }
    }
    void parallellRun(const World<T>& w, const Source<T>* source, Result<T>* result,
        const std::uint64_t expBeg, const std::uint64_t expEnd, std::uint64_t nJobs, ProgressBar<T>* progressbar) const;

    void parallellRunCtdi(const CTDIPhantom<T>& w, const CTSource<T>* source, Result<T>* result,
        const std::uint64_t expBeg, const std::uint64_t expEnd, std::uint64_t nJobs, ProgressBar<T>* progressbar) const;

    void computeResultVariance(Result<T>& res) const
    {
        //Computing variance by: Var[X] = E[X**2] - E[X]**2
        //and Var[A]+ Var[A] = 2Var[A]
        auto eBeg = res.dose.cbegin();
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
        }
    }

    void energyImpartedToDose(const World<T>& world, Result<T>& res, const T calibrationValue = 1)
    {
        auto spacing = world.spacing();
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
                return de > T { 0.0 } ? calibrationValue * var / voxelMass : T { 0.0 };
            });
    }

private:
    static constexpr T RUSSIAN_RULETTE_PROBABILITY()
    {
        return T { 0.8 };
    }
    static constexpr T RUSSIAN_RULETTE_ENERGY_THRESHOLD() // keV
    {
        return T { 10 };
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
};

#include "dxmc/source.h"

template <Floating T>
Result<T> Transport<T>::operator()(const World<T>& world, Source<T>* source, ProgressBar<T>* progressbar, bool calculateDose)
{

    Result<T> result(world.size());

    if (!source)
        return result;

    if (!world.isValid())
        return result;

    source->updateFromWorld(world);
    source->validate();
    if (!source->isValid())
        return result;

    const auto maxEnergy = source->maxPhotonEnergyProduced();
    constexpr T minEnergy { 1 };
    m_attenuationLut.generate(world, minEnergy, maxEnergy);

    const std::uint64_t totalExposures = source->totalExposures();

    const std::uint64_t nThreads = std::thread::hardware_concurrency();
    const std::uint64_t nJobs = std::max(nThreads, static_cast<std::uint64_t>(1));
    if (progressbar) {
        progressbar->setTotalExposures(totalExposures);
        progressbar->setDoseData(result.dose.data(), world.dimensions(), world.spacing());
    }
    const auto start = std::chrono::system_clock::now();
    parallellRun(world, source, &result, 0, totalExposures, nJobs, progressbar);
    result.simulationTime = std::chrono::system_clock::now() - start;

    if (progressbar) {
        progressbar->clearDoseData();
        if (progressbar->cancel()) {
            std::fill(result.dose.begin(), result.dose.end(), T { 0.0 });
            std::fill(result.nEvents.begin(), result.nEvents.end(), 0);
            std::fill(result.variance.begin(), result.variance.end(), T { 0.0 });
            return result;
        }
    }

    //compute variance
    computeResultVariance(result);

    if (calculateDose) {
        const auto calibrationValue = source->getCalibrationValue(progressbar);
        //energy imparted to dose
        energyImpartedToDose(world, result, calibrationValue);
    }
    return result;
}
template <Floating T>
Result<T> Transport<T>::operator()(const CTDIPhantom<T>& world, CTSource<T>* source, ProgressBar<T>* progressbar)
{
    Result<T> result(world.size());

    if (!source)
        return result;

    if (!world.isValid())
        return result;

    if (!source->isValid())
        return result;

    const auto maxEnergy = source->maxPhotonEnergyProduced();
    constexpr T minEnergy { 1 };
    m_attenuationLut.generate(world, minEnergy, maxEnergy);

    const std::uint64_t totalExposures = source->exposuresPerRotatition();

    const std::uint64_t nThreads = std::thread::hardware_concurrency();
    const std::uint64_t nJobs = std::max(nThreads, static_cast<std::uint64_t>(1));
    if (progressbar) {
        progressbar->setTotalExposures(totalExposures, "CTDI Calibration");
        progressbar->setDoseData(result.dose.data(), world.dimensions(), world.spacing(), ProgressBar<T>::Axis::Z);
    }
    const auto start = std::chrono::system_clock::now();
    parallellRunCtdi(world, source, &result, 0, totalExposures, nJobs, progressbar);
    result.simulationTime = std::chrono::system_clock::now() - start;

    //compute variance
    computeResultVariance(result);

    if (progressbar) {
        progressbar->clearDoseData();
    }

    energyImpartedToDose(world, result, T { 1.0 });
    return result;
}
template <Floating T>
void Transport<T>::parallellRun(const World<T>& w, const Source<T>* source, Result<T>* result,
    const std::uint64_t expBeg, const std::uint64_t expEnd, std::uint64_t nJobs, ProgressBar<T>* progressbar) const
{
    const std::uint64_t len = expEnd - expBeg;
    if ((len <= 1) || (nJobs <= 1)) {
        RandomState state;
        Exposure<T> exposure;
        const auto& worldBasis = w.directionCosines();
        for (std::size_t i = expBeg; i < expEnd; i++) {
            if (progressbar)
                if (progressbar->cancel())
                    return;
            source->getExposure(exposure, i);
            exposure.alignToDirectionCosines(worldBasis);
            transport(w, exposure, state, result);
            if (progressbar)
                progressbar->exposureCompleted();
        }
        return;
    }

    const auto threadLen = len / nJobs;
    std::vector<std::thread> jobs;
    jobs.reserve(nJobs - 1);
    std::uint64_t expBegThread = expBeg;
    std::uint64_t expEndThread = expBegThread + threadLen;
    for (std::size_t i = 0; i < nJobs - 1; ++i) {
        jobs.emplace_back(&Transport<T>::parallellRun, this, w, source, result, expBegThread, expEndThread, 1, progressbar);
        expBegThread = expEndThread;
        expEndThread = expBegThread + threadLen;
    }
    expEndThread = expEnd;
    parallellRun(w, source, result, expBegThread, expEndThread, 1, progressbar);
    for (std::size_t i = 0; i < nJobs - 1; ++i)
        jobs[i].join();
}

template <Floating T>
void Transport<T>::parallellRunCtdi(const CTDIPhantom<T>& w, const CTSource<T>* source, Result<T>* result,
    const std::uint64_t expBeg, const std::uint64_t expEnd, std::uint64_t nJobs, ProgressBar<T>* progressbar) const
{
    const std::uint64_t len = expEnd - expBeg;

    if ((len == 1) || (nJobs <= 1)) {
        RandomState state;
        Exposure<T> exposure;
        const auto& worldBasis = w.directionCosines();
        for (std::size_t i = expBeg; i < expEnd; i++) {
            if (progressbar)
                if (progressbar->cancel())
                    return;
            source->getExposure(exposure, i);
            const auto& pos = source->position();
            exposure.subtractPosition(pos); // aligning to center of phantom
            exposure.setPositionZ(0.0);
            exposure.alignToDirectionCosines(worldBasis);
            exposure.setBeamIntensityWeight(1.0);
            transport(w, exposure, state, result);
            if (progressbar)
                progressbar->exposureCompleted();
        }
        return;
    }

    const auto threadLen = len / nJobs;
    std::vector<std::thread> jobs;
    jobs.reserve(nJobs - 1);
    std::uint64_t expBegThread = expBeg;
    std::uint64_t expEndThread = expBegThread + threadLen;
    for (std::size_t i = 0; i < nJobs - 1; ++i) {
        jobs.emplace_back(&Transport<T>::parallellRunCtdi, this, w, source, result, expBegThread, expEndThread, 1, progressbar);
        expBegThread = expEndThread;
        expEndThread = expBegThread + threadLen;
    }
    expEndThread = expEnd;
    parallellRunCtdi(w, source, result, expBegThread, expEndThread, 1, progressbar);
    for (std::size_t i = 0; i < nJobs - 1; ++i)
        jobs[i].join();
}
}