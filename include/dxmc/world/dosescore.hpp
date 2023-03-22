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
#include "dxmc/floating.hpp"
#include <atomic>

namespace dxmc {

template <Floating T>
class DoseScore {
public:
    DoseScore() { }
    void scoreEnergy(T energy)
    {
        if (energy <= T { 0 })
            return;
        // Thread safe update of values with Kahn summation in case of large sums
        {
            auto aref = std::atomic_ref(m_energyImparted);
            const auto oldacc = aref.fetch_add(energy, std::memory_order_relaxed);
            const auto newacc = oldacc + energy;
            const auto r = energy - (newacc - oldacc);
            auto restref = std::atomic_ref(m_energyImparted_rest);
            restref.fetch_add(r, std::memory_order_relaxed);
            if (restref.load(std::memory_order_relaxed) > T { 1e5 }) {
                aref.fetch_add(restref.exchange(T { 0 }, std::memory_order_relaxed), std::memory_order_relaxed);
            }
        }
        {
            auto aref = std::atomic_ref(m_energyImpartedSquared);
            const auto esq = energy * energy;
            const auto oldacc = aref.fetch_add(esq, std::memory_order_relaxed);
            const auto newacc = oldacc + esq;
            const auto r = esq - (newacc - oldacc);
            auto restref = std::atomic_ref(m_energyImpartedSquared_rest);
            restref.fetch_add(r, std::memory_order_relaxed);
            if (restref.load(std::memory_order_relaxed) > T { 1e5 }) {
                aref.fetch_add(restref.exchange(T { 0 }, std::memory_order_relaxed), std::memory_order_relaxed);
            }
        }
        {
            auto aref = std::atomic_ref(m_nEvents);
            ++aref;
        }
    }
    T energyImparted() const
    {
        return m_energyImparted + m_energyImparted_rest;
    }

    T energyImpartedSquared() const
    {
        return m_energyImpartedSquared + m_energyImpartedSquared_rest;
    }

    T varianceEnergyImparted() const
    {
        if (numberOfEvents() > 0) {
            const auto e_exp = energyImparted() / numberOfEvents();
            const auto e2_exp = energyImpartedSquared() / numberOfEvents();
            return (e2_exp - e_exp * e_exp) * numberOfEvents();
        }
        return std::numeric_limits<T>::infinity();
    }

    T stdEnergyImparted() const
    {
        return std::sqrt(varianceEnergyImparted());
    }

    T relativeUncertainty() const
    {
        return T { 1.96 } * stdEnergyImparted() / energyImparted();
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
        m_energyImparted_rest = 0;
        m_energyImpartedSquared_rest = 0;
    }

private:
    std::uint64_t m_nEvents = 0;
    T m_energyImparted = 0;
    T m_energyImparted_rest = 0;
    T m_energyImpartedSquared = 0;
    T m_energyImpartedSquared_rest = 0;
};

template <>
class DoseScore<double> {
public:
    DoseScore() { }
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

    double varianceEnergyImparted() const
    {
        if (numberOfEvents() > 0) {
            const auto e_exp = energyImparted() / numberOfEvents();
            const auto e2_exp = energyImpartedSquared() / numberOfEvents();
            return (e2_exp - e_exp * e_exp) * numberOfEvents();
        }
        return std::numeric_limits<double>::infinity();
    }

    double stdEnergyImparted() const
    {
        return std::sqrt(varianceEnergyImparted());
    }

    double relativeUncertainty() const
    {
        return 1.96 * stdEnergyImparted() / energyImparted();
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
