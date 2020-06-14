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
#include "dxmc/progressbar.h"
#include "dxmc/source.h"
#include "dxmc/vectormath.h"
#include "dxmc/world.h"

#include <algorithm>
#include <chrono>
#include <memory>
#include <mutex>
namespace transport {

struct Result {
    std::vector<double> dose;
    std::vector<std::uint32_t> nEvents;
    std::vector<double> variance;
    std::chrono::duration<float> simulationTime;
};

double comptonScatter(Particle& particle, RandomState &seed, double& cosAngle);
double comptonScatterLivermore(Particle& particle, unsigned char materialIdx, const AttenuationLut& attLut, RandomState &seed, double& cosAngle);
void rayleightScatterLivermore(Particle& particle, unsigned char materialIdx, const AttenuationLut& attLut, RandomState &seed, double& cosAngle);
Result run(const World& world, Source* source, ProgressBar* progressBar = nullptr, bool calculateDose = true);
Result run(const CTDIPhantom& world, CTSource* source, ProgressBar* progressBar = nullptr);
}