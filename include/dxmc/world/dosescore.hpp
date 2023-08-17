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
#include "dxmc/world/energyscore.hpp"
#include <atomic>

namespace dxmc {

template <Floating T>
class DoseScore {
public:
    DoseScore() { }

    void addScoredEnergy(const EnergyScored<T>& energy, T volume, T density, T calibrationfactor = 1)
    {
        const auto mass = volume * density; // grams
        addScoredEnergy(energy, mass, calibrationfactor);
    }

    void addScoredEnergy(const EnergyScored<T>& energy, T mass, T calibrationfactor = 1)
    {
        const auto mass = volume * density; // grams
        const auto factor = calibrationfactor / mass;
        m_dose += energy.energyImparted() * factor;
        m_doseVariance += energy.variance() * factor * factor;
        m_nEvents += energy.numberOfEvents();
    }

    T dose() const
    {
        return m_dose;
    }

    T variance() const
    {
        return m_doseVariance;
    }

    T standardDeviation() const
    {
        return std::sqrt(variance());
    }

    T relativeUncertainty() const
    {
        return T { 1.96 } * standardDeviation() / dose();
    }

    std::uint64_t numberOfEvents() const
    {
        return m_nEvents;
    }

    void clear()
    {
        m_nEvents = 0;
        m_dose = 0;
        m_doseSquared = 0;
    }

private:
    T m_dose = 0;
    T m_doseVariance = 0;
    std::uint64_t m_nEvents = 0;
};
}
