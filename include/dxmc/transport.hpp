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
#include <chrono>
#include <format>
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

        const auto p = (m_nParticleCount * 100) / m_nParticles;
        if (m_nParticleCount > 0) {
            const auto remaining = (m_elapsed * (m_nParticles - m_nParticleCount)) / m_nParticleCount;
            return std::format("Time elapsed: {0} Estimated remaining time: {1} [{2}%]", human_time(elapsed), human_time(remaining), p);
        } else {
            return std::format("Time elapsed: {0} Estimated remaining time: NA [{1}%]", human_time(elapsed), p);
        }
    }

    static std::string human_time(const std::chrono::milliseconds& time)
    {
        if (time > std::chrono::hours(3))
            return std::format("{}", std::chrono::duration_cast<std::chrono::hours>(time));
        else if (time > std::chrono::minutes(3))
            return std::format("{}", std::chrono::duration_cast<std::chrono::minutes>(time));
        else
            return std::format("{}", std::chrono::duration_cast<std::chrono::seconds>(time));
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
    auto operator()(World2<T, Ws...>& world, const B& beam, TransportProgress* progress = nullptr) const
    {
        return run(world, beam, progress);
    }

protected:
    template <Floating T, BeamType<T> B, WorldItemType<T>... Ws>
    static void runWorker(World2<T, Ws...>& world, const B& beam, std::uint64_t exposureBegin, std::uint64_t exposureEnd, TransportProgress* progress = nullptr)
    {
        RandomState state;
        for (std::uint64_t eIdx = exposureBegin; eIdx != exposureEnd; ++eIdx) {
            const auto exposure = beam.exposure(eIdx);
            const auto nParticles = exposure.numberOfParticles();
            for (std::uint64_t i = 0; i != nParticles; ++i) {
                auto particle = exposure.sampleParticle(state);
                world.transport(particle, state);
            }
            if (progress)
                progress->addCompletedNumber(nParticles);
        }
    }

    template <Floating T, BeamType<T> B, WorldItemType<T>... Ws>
    void run(World2<T, Ws...>& world, const B& beam, TransportProgress* progress = nullptr) const
    {
        const auto nExposures = beam.numberOfExposures();
        const auto step = std::max(std::uint64_t { 1 }, nExposures / m_nThreads);
        std::vector<std::jthread> threads;
        std::uint64_t start = 0;
        auto stop = start + step;

        if (progress)
            progress->start(beam.numberOfParticles());

        for (std::size_t i = 0; i < m_nThreads - 1; ++i) {
            if (start < stop) {
                threads.emplace_back(Transport::template runWorker<T, B, Ws...>, std::ref(world), std::cref(beam), start, stop, progress);
                start = stop;
                stop += step;
            }
        }
        runWorker(world, beam, start, nExposures, progress);

        for (auto& thread : threads) {
            thread.join();
        }
    }

private:
    std::uint64_t m_nThreads = 1;
};

}
