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

Copyright 2022 Erlend Andersen
*/

#pragma once
#include "dxmc/dxmcrandom.hpp"
#include "dxmc/material/material.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/world/dosescore.hpp"
#include "dxmc/world/energyscore.hpp"
#include "dxmc/world/visualizationintersectionresult.hpp"
#include "dxmc/world/worldintersectionresult.hpp"

#include <algorithm>
#include <optional>

namespace dxmc {

class WorldItemBase {
public:
    virtual void translate(const std::array<double, 3>& dist) = 0;
    virtual std::array<double, 3> center() const = 0;
    virtual std::array<double, 6> AABB() const = 0;
    virtual WorldIntersectionResult intersect(const Particle& p) const = 0;
    virtual VisualizationIntersectionResult<WorldItemBase> intersectVisualization(const Particle& p) const = 0;
    virtual const EnergyScore& energyScored(std::size_t index = 0) const = 0;
    virtual void clearEnergyScored() = 0;
    virtual void addEnergyScoredToDoseScore(double calibration_factor = 1) = 0;
    virtual const DoseScore& doseScored(std::size_t index = 0) const = 0;
    virtual void clearDoseScored() = 0;
    virtual void transport(Particle& p, RandomState& state) = 0;
};
}