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

#include "dxmc/constants.hpp"
#include "dxmc/dxmcrandom.hpp"
#include "dxmc/floating.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/transportprogress.hpp"
#include "dxmc/vectormath.hpp"

#include <array>

namespace dxmc {

class IsotropicMonoEnergyBeamExposure {
public:
    IsotropicMonoEnergyBeamExposure(const std::array<double, 3>& pos, const std::array<std::array<double, 3>, 2>& dircosines, double energy, std::uint64_t N)
        : m_energy(energy)
        , m_pos(pos)
        , m_dirCosines(dircosines)
        , m_NParticles(N)
    {
        m_dir = vectormath::cross(m_dirCosines);
    }

    const std::array<double, 3>& position() const { return m_pos; }

    void setCollimationAngles(const std::array<double, 4>& angles)
    {
        for (std::size_t i = 0; i < angles.size(); ++i)
            m_collimationAngles[i] = angles[i];
    }

    std::uint64_t numberOfParticles() const { return m_NParticles; }

    Particle sampleParticle(RandomState& state) const noexcept
    {

        const auto angx = state.randomUniform(m_collimationAngles[0], m_collimationAngles[2]);
        const auto angy = state.randomUniform(m_collimationAngles[1], m_collimationAngles[3]);

        Particle p = { .pos = m_pos,
            .dir = particleDirection(angx, angy),
            .energy = m_energy,
            .weight = 1 };
        return p;
    }

protected:
    std::array<double, 3> particleDirection(double anglex, double angley) const
    {
        const auto dx = std::tan(anglex);
        const auto dy = std::tan(angley);
        const std::array pdir = {
            m_dirCosines[0][0] * dx + m_dirCosines[1][0] * dy + m_dir[0],
            m_dirCosines[0][1] * dx + m_dirCosines[1][1] * dy + m_dir[1],
            m_dirCosines[0][2] * dx + m_dirCosines[1][2] * dy + m_dir[2]
        };
        return vectormath::normalized(pdir);
    }

private:
    double m_energy = 60;
    std::array<double, 3> m_pos = { 0, 0, 0 };
    std::array<double, 3> m_dir = { 0, 0, 1 };
    std::array<std::array<double, 3>, 2> m_dirCosines = { { { 1, 0, 0 }, { 0, 1, 0 } } };
    std::array<double, 4> m_collimationAngles = { 0, 0, 0, 0 };
    std::uint64_t m_NParticles = 100;
};

class IsotropicMonoEnergyBeam {
public:
    IsotropicMonoEnergyBeam(const std::array<double, 3>& pos = { 0, 0, 0 }, const std::array<std::array<double, 3>, 2>& dircosines = { { { 1, 0, 0 }, { 0, 1, 0 } } }, double energy = 60)
        : m_energy(energy)
        , m_pos(pos)
    {
        setDirectionCosines(dircosines);
    }

    std::uint64_t numberOfExposures() const { return m_Nexposures; }
    void setNumberOfExposures(std::uint64_t n) { m_Nexposures = std::max(n, std::uint64_t { 1 }); }

    std::uint64_t numberOfParticles() const { return m_Nexposures * m_particlesPerExposure; }
    void setNumberOfParticlesPerExposure(std::uint64_t n) { m_particlesPerExposure = n; }

    void setPosition(const std::array<double, 3>& pos) { m_pos = pos; }
    const std::array<double, 3>& position() const { return m_pos; }

    const std::array<std::array<double, 3>, 2>& directionCosines() const
    {
        return m_dirCosines;
    }

    void setDirectionCosines(const std::array<std::array<double, 3>, 2>& dir)
    {
        m_dirCosines = dir;
        vectormath::normalize(m_dirCosines[0]);
        vectormath::normalize(m_dirCosines[1]);
    }
    void setDirectionCosines(const std::array<double, 3>& xdir, const std::array<double, 3>& ydir)
    {
        m_dirCosines[0] = xdir;
        m_dirCosines[1] = ydir;
        vectormath::normalize(m_dirCosines[0]);
        vectormath::normalize(m_dirCosines[1]);
    }

    void setEnergy(double energy) { m_energy = std::min(std::max(MIN_ENERGY(), energy), MAX_ENERGY()); }

    double energy() const { return m_energy; }

    const std::array<double, 4>& collimationAngles() const { return m_collimationAngles; }

    void setCollimationAngles(const std::array<double, 4>& angles) { m_collimationAngles = angles; }
    void setCollimationAngles(double minX, double minY, double maxX, double maxY)
    {
        m_collimationAngles[0] = minX;
        m_collimationAngles[1] = minY;
        m_collimationAngles[2] = maxX;
        m_collimationAngles[3] = maxY;
    }

    IsotropicMonoEnergyBeamExposure exposure(std::size_t i) const noexcept
    {
        IsotropicMonoEnergyBeamExposure exp(m_pos, m_dirCosines, m_energy, m_particlesPerExposure);
        exp.setCollimationAngles(m_collimationAngles);
        return exp;
    }

    double calibrationFactor(TransportProgress* progress = nullptr) const noexcept
    {
        return 1;
    }

private:
    double m_energy = 60;
    std::array<double, 3> m_pos = { 0, 0, 0 };
    std::array<std::array<double, 3>, 2> m_dirCosines = { { { 1, 0, 0 }, { 0, 1, 0 } } };
    std::array<double, 4> m_collimationAngles = { 0, 0, 0, 0 };
    std::uint64_t m_Nexposures = 100;
    std::uint64_t m_particlesPerExposure = 100;
};

}