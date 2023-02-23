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

#include <algorithm>
#include <atomic>
#include <optional>

namespace dxmc {

template <Floating T>
class DoseScore {
public:
    DoseScore() { }
    void scoreEnergy(T energy)
    {
        // threadsafe update
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
    T energyImparted() const
    {
        return m_energyImparted;
    }

private:
    std::uint64_t m_nEvents = 0;
    T m_energyImparted = 0;
    T m_energyImpartedSquared = 0;
};

template <Floating T>
class WorldItemBase {
public:
    virtual void translate(const std::array<T, 3>& dist) = 0;
    virtual std::array<T, 3> center() const = 0;
    virtual std::array<T, 6> AABB() const = 0;
    virtual std::optional<T> intersectForward(const Particle<T>& p) const = 0;
    virtual const DoseScore<T>& dose(std::size_t index = 0) const = 0;
    virtual void transport(Particle<T>& p, RandomState& state) = 0;

protected:
private:
};
}