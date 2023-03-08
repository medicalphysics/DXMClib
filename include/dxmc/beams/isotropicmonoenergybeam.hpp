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
class IsotropicMonoEnergyBeamExposure {
public:
    IsotropicMonoEnergyBeamExposure(const std::array<T, 3>& pos, const std::array<std::array<T, 3>, 2>& dircosines, T energy, std::uint64_t N)
        : m_pos(pos)
        , m_dirCosines(dircosines)
        , m_energy(energy)
        , m_NParticles(N)
    {
    }

    void setCollimationAngles(const std::array<T, 2>& angles) { m_collimationAngles = angles; }

    std::uint64_t numberOfParticles() const { return m_NParticles; }

    Particle<T> sampleParticle(RandomState& state) const noexcept
    {
        auto dir = vectormath::cross(m_dirCosines[0], m_dirCosines[1]);
        const auto theta_x = state.randomUniform(-m_collimationAngles[0], m_collimationAngles[0]);
        dir = vectormath::rotate(dir, m_dirCosines[1], theta_x);
        const auto theta_y = state.randomUniform(-m_collimationAngles[1], m_collimationAngles[1]);
        dir = vectormath::rotate(dir, m_dirCosines[0], theta_y);

        Particle<T> p = { .pos = m_pos,
            .dir = dir,
            .energy = m_energy,
            .weight = T { 1 } };
        return p;
    }

private:
    T m_energy = 60;
    std::array<T, 3> m_pos = { 0, 0, 0 };
    std::array<std::array<T, 3>, 2> m_dirCosines = { 1, 0, 0, 0, 1, 0 };
    std::array<T, 2> m_collimationAngles = { 0, 0 };
    std::uint64_t m_NParticles = 100;
};

template <Floating T>
class IsotropicMonoEnergyBeam {
public:
    IsotropicMonoEnergyBeam(const std::array<T, 3>& pos = { 0, 0, 0 }, const std::array<std::array<T, 3>, 2>& dircosines = { 1, 0, 0, 0, 1, 0 }, T energy = 60)
        : m_pos(pos)
        , m_energy(energy)
    {
        setDirectionCosines(dircosines);
    }

    std::uint64_t numberOfExposures() const { return m_Nexposures; }
    void setNumberOfExposures(std::uint64_t n) { m_Nexposures = std::max(n, std::uint64_t { 1 }); }

    std::uint64_t numberOfParticles() const { return m_Nexposures * m_particlesPerExposure; }
    void setNumberOfParticlesPerExposure(std::uint64_t n) { m_particlesPerExposure = n; }

    void setPosition(const std::array<T, 3>& pos)
    {
        m_pos = pos;
    }

    void setDirectionCosines(const std::array<std::array<T, 3>, 2>& dir)
    {
        m_dirCosines = dir;
        vectormath::normalize(m_dirCosines[0]);
        vectormath::normalize(m_dirCosines[1]);
    }

    void setCollimationAngles(const std::array<T, 2>& angles) { m_collimationAngles = angles; }

    IsotropicMonoEnergyBeamExposure<T> exposure(std::size_t i) const noexcept
    {
        
        IsotropicMonoEnergyBeamExposure<T> exp(m_pos, m_dirCosines, m_energy, m_particlesPerExposure);
        exp.setCollimationAngles(m_collimationAngles);
        return exp;
    }

private:
    T m_energy = 60;
    std::array<T, 3> m_pos = { 0, 0, 0 };
    std::array<std::array<T, 3>, 2> m_dirCosines = { 1, 0, 0, 0, 1, 0 };
    std::array<T, 2> m_collimationAngles = { 0, 0 };
    std::uint64_t m_Nexposures = 100;
    std::uint64_t m_particlesPerExposure = 100;
};

}