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
#include <functional>
#include <thread>
#include <chrono>

namespace dxmc {
template <Floating T>
class Transport {

public:
    Transport()
    {
        const std::uint64_t hvc = std::thread::hardware_concurrency();
        m_nThreads = std::max(hvc, std::uint64_t { 1 });
    }

    std::uint64_t numberOfThreads() const { return m_nThreads; }
    void setNumberOfThreads(std::uint64_t n) { m_nThreads = std::max(n, std::uint64_t { 1 }); }

    template <BeamType<T> B, WorldItemType<T>... Ws>
    auto operator()(World2<T, Ws...>& world, const B& beam) const
    {
        return run(world, beam);
    }

protected:
    template <BeamType<T> B, WorldItemType<T>... Ws>
    static void runWorker(World2<T, Ws...>& world, const B& beam, std::uint64_t exposureBegin, std::uint64_t exposureEnd)
    {
        RandomState state;
        for (std::uint64_t eIdx = exposureBegin; eIdx != exposureEnd; ++eIdx) {
            const auto exposure = beam.exposure(eIdx);
            const auto nParticles = exposure.numberOfParticles();
            for (std::uint64_t i = 0; i != nParticles; ++i) {
                auto particle = exposure.sampleParticle(state);
                world.transport(particle, state);
            }
        }
    }

    template <BeamType<T> B, WorldItemType<T>... Ws>
    auto run(World2<T, Ws...>& world, const B& beam) const
    {
        const auto nExposures = beam.numberOfExposures();
        const auto step = std::max(std::uint64_t { 1 }, nExposures / m_nThreads);
        std::vector<std::jthread> threads;
        std::uint64_t start = 0;
        auto stop = start + step;

        const auto time_start = std::chrono::high_resolution_clock::now();

        for (std::size_t i = 0; i < m_nThreads - 1; ++i) {
            if (start < stop) {
                threads.emplace_back(Transport<T>::template runWorker<B, Ws...>, std::ref(world), std::cref(beam), start, stop);
                start = stop;
                stop += step;
            }
        }
        runWorker(world, beam, start, nExposures);

        for (auto& thread : threads) {
            thread.join();
        }

        const auto time_end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration_cast<std::chrono::milliseconds> (time_end - time_start);
    }

private:
    std::uint64_t m_nThreads = 1;
};

}
