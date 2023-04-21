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

    void setCollimationAngles(const std::array<T, 4>& angles)
    {
        for (std::size_t i = 0; i < angles.size(); ++i)
            m_collimationAngles[i] = angles[i];
    }

    std::uint64_t numberOfParticles() const { return m_NParticles; }

    Particle<T> sampleParticle(RandomState& state) const noexcept
    {
        auto dir = vectormath::cross(m_dirCosines[0], m_dirCosines[1]);
        const auto angx = state.randomUniform(m_collimationAngles[0], m_collimationAngles[2]);
        const auto angy = state.randomUniform(m_collimationAngles[1], m_collimationAngles[3]);
        
        const auto sinx = std::sin(angx);
        const auto siny = std::sin(angy);
        const auto sinz = std::sqrt(1 - sinx * sinx - siny * siny);
        std::array pdir = {
            m_dirCosines[0][0] * sinx + m_dirCosines[1][0] * siny + dir[0] * sinz,
            m_dirCosines[0][1] * sinx + m_dirCosines[1][1] * siny + dir[1] * sinz,
            m_dirCosines[0][2] * sinx + m_dirCosines[1][2] * siny + dir[2] * sinz
        };

        //auto pdir = vectormath::rotate(vectormath::rotate(dir, m_dirCosines[1], angx), m_dirCosines[0], angy);
        
        Particle<T> p = { .pos = m_pos,
            .dir = pdir,
            .energy = m_energy,
            .weight = T { 1 } };
        return p;
    }

private:
    T m_energy = 60;
    std::array<T, 3> m_pos = { 0, 0, 0 };
    std::array<std::array<T, 3>, 2> m_dirCosines = { 1, 0, 0, 0, 1, 0 };
    std::array<T, 4> m_collimationAngles = { 0, 0, 0, 0 };
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

    void setPosition(const std::array<T, 3>& pos) { m_pos = pos; }

    void setDirectionCosines(const std::array<std::array<T, 3>, 2>& dir)
    {
        m_dirCosines = dir;
        vectormath::normalize(m_dirCosines[0]);
        vectormath::normalize(m_dirCosines[1]);
    }
    void setDirectionCosines(const std::array<T, 3>& xdir, const std::array<T, 3>& ydir)
    {
        m_dirCosines[0] = xdir;
        m_dirCosines[1] = ydir;
        vectormath::normalize(m_dirCosines[0]);
        vectormath::normalize(m_dirCosines[1]);
    }

    void setEnergy(T energy) { m_energy = std::min(std::max(MIN_ENERGY<T>(), energy), MAX_ENERGY<T>()); }

    T energy() { return m_energy; }

    void setCollimationAngles(const std::array<T, 4>& angles) { m_collimationAngles = angles; }
    void setCollimationAngles(T minX, T minY, T maxX, T maxY)
    {
        m_collimationAngles[0] = minX;
        m_collimationAngles[1] = minY;
        m_collimationAngles[2] = maxX;
        m_collimationAngles[3] = maxY;
    }

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
    std::array<T, 4> m_collimationAngles = { 0, 0, 0, 0 };
    std::uint64_t m_Nexposures = 100;
    std::uint64_t m_particlesPerExposure = 100;
};

}