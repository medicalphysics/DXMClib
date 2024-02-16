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
#include <atomic>

namespace dxmc {

class EnergyScore {
public:
    EnergyScore() { }
    void scoreEnergy(double energy)
    {
        if (energy <= 0.0)
            return;

        {
            auto aref = std::atomic_ref(m_energyImparted);
            aref.fetch_add(energy, std::memory_order_relaxed);
        }
        {
            auto aref = std::atomic_ref(m_energyImpartedSquared);
            aref.fetch_add(energy * energy, std::memory_order_relaxed);
        }
        {
            auto aref = std::atomic_ref(m_nEvents);
            ++aref;
        }
    }
    double energyImparted() const
    {
        return m_energyImparted;
    }

    double energyImpartedSquared() const
    {
        return m_energyImpartedSquared;
    }

    double variance() const
    {
        if (numberOfEvents() > 0) {
            const auto e_exp = energyImparted() / numberOfEvents();
            const auto e2_exp = energyImpartedSquared() / numberOfEvents();
            return (e2_exp - e_exp * e_exp) * numberOfEvents();
        }

        return 0;
    }

    double standardDeviation() const
    {
        return std::sqrt(variance());
    }

    double relativeUncertainty() const
    {
        return 1.96 * standardDeviation() / energyImparted();
    }

    std::uint64_t numberOfEvents() const
    {
        return m_nEvents;
    }

    void clear()
    {
        m_nEvents = 0;
        m_energyImparted = 0;
        m_energyImpartedSquared = 0;
    }

private:
    std::uint64_t m_nEvents = 0;
    double m_energyImparted = 0;
    double m_energyImpartedSquared = 0;
};
}
