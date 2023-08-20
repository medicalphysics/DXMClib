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
#include "dxmc/transportprogress.hpp"
#include "dxmc/world/world.hpp"

#include <algorithm>
#include <atomic>
#include <functional>
#include <thread>

namespace dxmc {

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
            const auto beamCalibrationFactor = beam.calibrationFactor(progress);
            world.addEnergyScoredToDoseScore(beamCalibrationFactor);
        } else {
            world.addEnergyScoredToDoseScore();
        }
    }

private:
    std::uint64_t m_nThreads = 1;
};

}
