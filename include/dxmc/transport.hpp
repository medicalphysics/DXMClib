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

#include "dxmc/beams/beamtype.hpp"
#include "dxmc/dxmcrandom.hpp"
#include "dxmc/floating.hpp"
#include "dxmc/world/world.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <functional>
#include <string>
#include <thread>

namespace dxmc {
class TransportProgress {
public:
    TransportProgress() { }
    TransportProgress(std::uint64_t n_particles)
    {
        start(n_particles);
    }

    void start(std::uint64_t N)
    {
        m_nParticles = N;
        m_nParticleCount = 0;
        m_start = std::chrono::high_resolution_clock::now();
    }
    void addCompletedNumber(std::uint64_t N)
    {
        auto aref = std::atomic_ref(m_nParticleCount);
        aref.fetch_add(N, std::memory_order_relaxed);

        auto dref = std::atomic_ref(m_elapsed);
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - m_start);
        dref.exchange(duration);

        if (aref.load() >= m_nParticles) {
            m_end = std::chrono::high_resolution_clock::now();
        }
    }

    std::string humanTotalTime() const
    {
        if (m_nParticleCount == m_nParticles) {
            return human_time(std::chrono::duration_cast<std::chrono::milliseconds>(m_end - m_start));
        }
        return "Not done yet";
    }
    std::chrono::milliseconds totalTime() const
    {
        return std::chrono::duration_cast<std::chrono::milliseconds>(m_end - m_start);
    }

    std::string message() const
    {
        const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - m_start);

        const auto p = std::to_string((m_nParticleCount * 100) / m_nParticles);

        std::string message = "Time elapsed: " + human_time(elapsed) + "Estimated remaining time: ";

        if (m_nParticleCount > 0) {
            const auto remaining = (m_elapsed * (m_nParticles - m_nParticleCount)) / m_nParticleCount;
            message += human_time(remaining) + " [" + p + "%]";
        } else {
            message += "NA [" + p + "%]";
        }
        return message;
    }

    static std::string human_time(const std::chrono::milliseconds& time)
    {
        if (time > std::chrono::hours(3))
            return std::to_string(std::chrono::duration_cast<std::chrono::hours>(time).count());
        else if (time > std::chrono::minutes(3))
            return std::to_string(std::chrono::duration_cast<std::chrono::minutes>(time).count());
        else
            return std::to_string(std::chrono::duration_cast<std::chrono::seconds>(time).count());
    }

private:
    std::uint64_t m_nParticles = 0;
    std::uint64_t m_nParticleCount = 0;
    std::chrono::time_point<std::chrono::high_resolution_clock> m_start;
    std::chrono::time_point<std::chrono::high_resolution_clock> m_end;
    std::chrono::milliseconds m_elapsed;
};

class Transport {
public:
    Transport()
    {
        const std::uint64_t hvc = std::thread::hardware_concurrency();
        m_nThreads = std::max(hvc, std::uint64_t { 1 });
    }

    std::uint64_t numberOfThreads() const { return m_nThreads; }
    void setNumberOfThreads(std::uint64_t n) { m_nThreads = std::max(n, std::uint64_t { 1 }); }

    template <Floating T, BeamType<T> B, WorldItemType<T>... Ws>
    auto operator()(World<T, Ws...>& world, const B& beam, TransportProgress* progress = nullptr, bool useBeamCalibration = true) const
    {
        return run(world, beam, progress, useBeamCalibration);
    }

protected:
    template <Floating T, BeamType<T> B, WorldItemType<T>... Ws>
    static void runWorker(World<T, Ws...>& world, const B& beam, RandomState& state, std::atomic<std::uint64_t>& exposureStart, std::uint64_t exposureEnd, TransportProgress* progress = nullptr)
    {
        auto n = exposureStart.fetch_add(1);
        while (n < exposureEnd) {
            const auto exposure = beam.exposure(n);
            const auto nParticles = exposure.numberOfParticles();
            for (std::uint64_t i = 0; i != nParticles; ++i) {
                auto particle = exposure.sampleParticle(state);
                world.transport(particle, state);
            }
            if (progress) {
                progress->addCompletedNumber(nParticles);
            }
            n = exposureStart.fetch_add(1);
        }
    }

    template <Floating T, BeamType<T> B, WorldItemType<T>... Ws>
    void run(World<T, Ws...>& world, const B& beam, TransportProgress* progress = nullptr, bool useBeamCalibration = true) const
    {
        // clearing scored energy before run
        world.clearEnergyScored();

        const auto nExposures = beam.numberOfExposures();
        std::vector<std::jthread> threads;
        threads.reserve(m_nThreads - 1);
        std::vector<RandomState> states(m_nThreads - 1);
        std::atomic<std::uint64_t> start(0);

        if (progress)
            progress->start(beam.numberOfParticles());

        for (std::size_t i = 0; i < m_nThreads - 1; ++i) {
            threads.emplace_back(Transport::template runWorker<T, B, Ws...>, std::ref(world), std::cref(beam), std::ref(states[i]), std::ref(start), nExposures, progress);
        }
        RandomState state;
        runWorker(world, beam, state, start, nExposures, progress);

        for (auto& thread : threads) {
            thread.join();
        }
        if (useBeamCalibration) {
            const auto beamCalibrationFactor = beam.calibrationFactor();
            world.addEnergyScoredToDoseScore(beamCalibrationFactor);
        } else {
            world.addEnergyScoredToDoseScore();
        }
    }

private:
    std::uint64_t m_nThreads = 1;
};

}
