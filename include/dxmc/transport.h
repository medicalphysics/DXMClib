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
#include "dxmc/dxmcrandom.h"
#include "dxmc/exposure.h"
#include "dxmc/particle.h"
#include "dxmc/progressbar.h"
#include "dxmc/source.h"
#include "dxmc/vectormath.h"
#include "dxmc/world.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <memory>

namespace dxmc {

namespace transport {

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

    struct Result {
        std::vector<double> dose;
        std::vector<std::uint32_t> nEvents;
        std::vector<double> variance;
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

    double comptonScatter(Particle<double>& particle, RandomState& seed, double& cosAngle);
    double comptonScatterLivermore(Particle<double>& particle, unsigned char materialIdx, const AttenuationLut& attLut, RandomState& seed, double& cosAngle);
    void rayleightScatterLivermore(Particle<double>& particle, unsigned char materialIdx, const AttenuationLut& attLut, RandomState& seed, double& cosAngle);
    Result run(const World& world, Source* source, ProgressBar* progressBar = nullptr, bool calculateDose = true);
    Result run(const CTDIPhantom& world, CTSource* source, ProgressBar* progressBar = nullptr);
}
}