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
#include "dxmc/floating.hpp"
#include "dxmc/material/material.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/world/worldintersectionresult.hpp"

#include <algorithm>
#include <atomic>
#include <optional>

namespace dxmc {

template <Floating T, bool KAHN_SUMMATION = true>
class DoseScore {
public:
    DoseScore() { }
    void scoreEnergy(T energy)
    {
        // threadsafe update
        if constexpr (KAHN_SUMMATION) {
            auto aref = std::atomic_ref(m_energyImparted);
            volatile const auto oldacc = aref.fetch_add(energy);
            volatile const auto newacc = oldacc + energy;
            volatile const auto r = energy - (newacc - oldacc); // Recover rounding error using Fast2Sum, assumes |oldacc|>=|val|
            auto restref = std::atomic_ref(m_energyImparted_rest);
            restref.fetch_add(r);
        } else {
            auto aref = std::atomic_ref(m_energyImparted);
            aref.fetch_add(energy);
        }

        if constexpr (KAHN_SUMMATION) {
            auto aref = std::atomic_ref(m_energyImpartedSquared);
            const auto esq = energy * energy;
            const auto oldacc = aref.fetch_add(esq);
            const auto newacc = oldacc + esq;
            const auto r = esq - (newacc - oldacc); // Recover rounding error using Fast2Sum, assumes |oldacc|>=|val|
            auto restref = std::atomic_ref(m_energyImpartedSquared_rest);
            restref.fetch_add(r);
        } else {
            auto aref = std::atomic_ref(m_energyImpartedSquared);
            aref.fetch_add(energy * energy);
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
    std::uint64_t numberOfEvents() const
    {
        return m_nEvents;
    }

private:
    std::uint64_t m_nEvents = 0;
    T m_energyImparted = 0;
    T m_energyImpartedSquared = 0;
    T m_energyImparted_rest = 0;
    T m_energyImpartedSquared_rest = 0;
};

template <Floating T>
class WorldItemBase {
public:
    virtual void translate(const std::array<T, 3>& dist) = 0;
    virtual std::array<T, 3> center() const = 0;
    virtual std::array<T, 6> AABB() const = 0;
    virtual WorldIntersectionResult<T> intersect(const Particle<T>& p) const = 0;

    virtual const DoseScore<T>& dose(std::size_t index = 0) const = 0;
    virtual void transport(Particle<T>& p, RandomState& state) = 0;

protected:
private:
};
}