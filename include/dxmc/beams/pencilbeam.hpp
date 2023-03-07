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
#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"

#include <array>

namespace dxmc {

template <Floating T>
class PencilBeamExposure {
public:
    PencilBeamExposure(const std::array<T, 3>& pos, const std::array<T, 3>& dir, T energy, std::uint64_t N)
        : m_pos(pos)
        , m_dir(dir)
        , m_energy(energy)
        , m_NParticles(N)
    {
    }

    std::uint64_t numberOfParticles() const { return m_NParticles; }

    Particle<T> sampleParticle(RandomState& state) const noexcept
    {
        Particle<T> p = { .pos = m_pos,
            .dir = m_dir,
            .energy = m_energy,
            .weight = T { 1 } };
        return p;
    }

private:
    T m_energy = 60;
    std::array<T, 3> m_pos = { 0, 0, 0 };
    std::array<T, 3> m_dir = { 0, 0, 1 };
    std::uint64_t m_NParticles = 100;
};

template <Floating T>
class PencilBeam {
public:
    PencilBeam(const std::array<T, 3>& pos = { 0, 0, 0 }, const std::array<T, 3>& dir = { 0, 0, 1 }, T energy = 60)
        : m_pos(pos)
        , m_dir(dir)
        , m_energy(energy)
    {
    }

    std::uint64_t numberOfExposures() const { return m_Nexposures; }

    std::uint64_t numberOfParticles() const { return m_Nexposures * m_particlesPerExposure; }
    void setNumberOfParticlesPerExposure(std::uint64_t n) { m_particlesPerExposure = n; }

    void setPosition(const std::array<T, 3>&& pos)
    {
        m_pos = pos;
    }

    void setDirection(const std::array<T, 3>& dir)
    {
        m_dir = dir;
        vectormath::normalize(m_dir);
    }

    PencilBeamExposure<T> exposure(std::size_t i) const noexcept
    {
        PencilBeamExposure<T> exp(m_pos, m_dir, m_energy, m_particlesPerExposure);
        return exp;
    }

private:
    T m_energy = 60;
    std::array<T, 3> m_pos = { 0, 0, 0 };
    std::array<T, 3> m_dir = { 0, 0, 1 };
    std::uint64_t m_Nexposures = 100;
    std::uint64_t m_particlesPerExposure = 100;
};

}