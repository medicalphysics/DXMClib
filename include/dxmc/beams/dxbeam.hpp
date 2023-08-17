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

#include "dxmc/beams/tube/tube.hpp"
#include "dxmc/dxmcrandom.hpp"
#include "dxmc/floating.hpp"
#include "dxmc/material/material.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"

#include <array>

namespace dxmc {

template <Floating T>
class DXBeamExposure {
public:
    DXBeamExposure(const std::array<T, 3>& pos, const std::array<std::array<T, 3>, 2>& dircosines, std::uint64_t N, T weight,
        const std::array<T, 4>& collimationAngles, const SpecterDistribution<T> specter)
        : m_pos(pos)
        , m_dirCosines(dircosines)
        , m_NParticles(N)
        , m_weight(weight)
        , m_collimationAngles(collimationAngles)
        , m_specter(specter)
    {
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

        Particle<T> p = {
            .pos = m_pos,
            .dir = pdir,
            .energy = m_specter.sampleValue(state),
            .weight = m_weight
        };
        return p;
    }

private:
    std::array<T, 3> m_pos = { 0, 0, 0 };
    std::array<std::array<T, 3>, 2> m_dirCosines = { 1, 0, 0, 0, 1, 0 };
    std::array<T, 4> m_collimationAngles = { 0, 0, 0, 0 };
    std::uint64_t m_NParticles = 100;
    T m_weight = 1;
    SpecterDistribution<T> m_specter;
};

template <Floating T>
class DXBeam {
public:
    DXBeam(const std::array<T, 3>& pos = { 0, 0, 0 }, const std::array<std::array<T, 3>, 2>& dircosines = { 1, 0, 0, 0, 1, 0 })
        : m_pos(pos)
    {
        setDirectionCosines(dircosines);
        tubeChanged();
    }

    std::uint64_t numberOfExposures() const { return m_Nexposures; }
    void setNumberOfExposures(std::uint64_t n) { m_Nexposures = std::max(n, std::uint64_t { 1 }); }

    std::uint64_t numberOfParticles() const { return m_Nexposures * m_particlesPerExposure; }
    void setNumberOfParticlesPerExposure(std::uint64_t n) { m_particlesPerExposure = n; }

    const std::array<T, 3>& position() const { return m_pos; }
    void setPosition(const std::array<T, 3>& pos) { m_pos = pos; }

    const std::array<std::array<T, 3>, 2>& directionCosines() const
    {
        return m_dirCosines;
    }

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

    T DAPvalue() const { return m_measuredDAP; }
    void setDAPvalue(T dap) { m_measuredDAP = std::abs(dap); }

    const Tube<T>& tube() const { return m_tube; }
    void setTube(const Tube<T>&& tube)
    {
        m_tube = tube;
        tubeChanged();
    }
    void setTubeVoltage(T voltage)
    {
        m_tube.setVoltage(voltage);
        tubeChanged();
    }
    void setTubeAnodeAngle(T ang)
    {
        m_tube.setAnodeAngle(ang);
        tubeChanged();
    }
    void setTubeAnodeAngleDeg(T ang)
    {
        m_tube.setAnodeAngleDeg(ang);
        tubeChanged();
    }
    void addTubeFiltrationMaterial(std::size_t Z, T mm)
    {
        m_tube.addFiltrationMaterial(Z, mm);
        tubeChanged();
    }
    void clearTubeFiltrationMaterials(std::size_t Z, T mm)
    {
        m_tube.clearFiltrationMaterials();
        tubeChanged();
    }
    void setTubeEnergyResolution(T energyResolution)
    {
        m_tube.setEnergyResolution(energyResolution);
        tubeChanged();
    }

    const std::array<T, 4>& collimationAngles() const { return m_collimationAngles; }
    void setCollimationAngles(const std::array<T, 4>& angles) { m_collimationAngles = angles; }
    void setCollimationAngles(T minX, T minY, T maxX, T maxY)
    {
        m_collimationAngles[0] = minX;
        m_collimationAngles[1] = minY;
        m_collimationAngles[2] = maxX;
        m_collimationAngles[3] = maxY;
    }

    DXBeamExposure<T> exposure(std::size_t i) const noexcept
    {
        DXBeamExposure<T> exp(m_pos, m_dirCosines, m_particlesPerExposure, m_weight, m_collimationAngles, m_specter);
        return exp;
    }

    T calibrationFactor() const
    {
        const auto energies = m_tube.getEnergy();
        const auto weights = m_tube.getSpecter(energies, true);
        auto air_cand = Material<T, 5, 0>::byNistName("Air, Dry (near sea level)");
        if (!air_cand)
            return 0;
        const auto& air = air_cand.value();

        const auto kerma_per_history = std::transform_reduce(std::execution::par_unseq, energies.cbegin(), energies.cend(), weights.cbegin(), T { 0 }, std::plus<>, [&air](const auto e, const auto w) {
            const auto uen = air.massEnergyTransferAttenuation(e);
            return w * e * uen;
        });

        const auto kerma_total = kerma_per_history * numberOfParticles();
        return m_measuredDAP / kerma_total;
    }

protected:
    void tubeChanged()
    {
        const auto energies = m_tube.getEnergy();
        const auto weights = m_tube.getSpecter(energies, false);
        m_specter = SpecterDistribution(energies, weights);
    }

private:
    std::array<T, 3> m_pos = { 0, 0, 0 };
    std::array<std::array<T, 3>, 2> m_dirCosines = { 1, 0, 0, 0, 1, 0 };
    std::array<T, 4> m_collimationAngles = { 0, 0, 0, 0 };
    std::uint64_t m_Nexposures = 100;
    std::uint64_t m_particlesPerExposure = 100;
    T m_weight = 1;
    T m_measuredDAP = 1;
    Tube<T> m_tube;
    SpecterDistribution<T> m_specter;
};

}