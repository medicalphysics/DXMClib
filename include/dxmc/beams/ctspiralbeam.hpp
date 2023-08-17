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
#include <cmath>

namespace dxmc {

template <Floating T>
class CTSpiralBeamExposure {
public:
    CTSpiralBeamExposure(const std::array<T, 3>& pos, const std::array<std::array<T, 3>, 2>& dircosines, std::uint64_t N, T weight,
        const std::array<T, 2>& collimationAngles, const SpecterDistribution<T> specter)
        : m_pos(pos)
        , m_dirCosines(dircosines)
        , m_NParticles(N)
        , m_weight(weight)
        , m_specter(specter)
        , m_collimationAngles(collimationAngles)
    {
        m_dir = vectormath::cross(m_dirCosines[0], m_dirCosines[1]);
    }

    std::uint64_t numberOfParticles() const
    {
        return m_NParticles;
    }

    Particle<T> sampleParticle(RandomState& state) const noexcept
    {

        const auto angx = state.randomUniform(-m_collimationAngles[0], m_collimationAngles[0]);
        const auto angy = state.randomUniform(-m_collimationAngles[1], m_collimationAngles[1]);

        const auto sinx = std::sin(angx);
        const auto siny = std::sin(angy);
        const auto sinz = std::sqrt(1 - sinx * sinx - siny * siny);
        std::array pdir = {
            m_dirCosines[0][0] * sinx + m_dirCosines[1][0] * siny + m_dir[0] * sinz,
            m_dirCosines[0][1] * sinx + m_dirCosines[1][1] * siny + m_dir[1] * sinz,
            m_dirCosines[0][2] * sinx + m_dirCosines[1][2] * siny + m_dir[2] * sinz
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
    std::array<T, 3> m_dir = { 0, 0, 1 };
    std::array<std::array<T, 3>, 2> m_dirCosines = { 1, 0, 0, 0, 1, 0 };
    std::array<T, 2> m_collimationAngles = { 0, 0 };
    std::uint64_t m_NParticles = 100;
    T m_weight = 1;
    SpecterDistribution<T> m_specter;
};

template <Floating T>
class CTSpiralBeam {
public:
    CTSpiralBeam(const std::array<T, 3>& start_pos = { 0, 0, 0 }, const std::array<T, 3>& stop_pos = { 0, 0, 0 })
        : m_start(start_pos)
        , m_stop(stop_pos)
    {
        tubeChanged();
    }

    std::uint64_t numberOfExposures() const
    {
        const auto direction = vectormath::subtract(m_stop, m_start);
        const auto dz = m_pitch * m_collimation;
        const auto total_rot_angle = vectormath::lenght(direction) * (PI_VAL<T>() * 2) / dz;
        auto N_angles = static_cast<std::uint64_t>(std::ceil(total_rot_angle / m_stepAngle));
        return N_angles;
    }

    std::uint64_t numberOfParticles() const { return numberOfExposures() * m_particlesPerExposure; }
    void setNumberOfParticlesPerExposure(std::uint64_t n) { m_particlesPerExposure = n; }

    const std::array<T, 3>& start() const { return m_start; }
    const std::array<T, 3>& stop() const { return m_stop; }
    void setStartStop(const std::array<T, 3>& start, const std::array<T, 3>& stop)
    {
        m_start = start;
        m_stop = stop;
    }

    T collimation() const { return m_collimation; }
    void setCollimation(T coll_cm)
    {
        // collimation must be larger than 1 mm (0.1 cm)
        m_collimation = std::max(std::abs(coll_cm), T { 0.1 });
    }

    T sourceDetectorDistance() const { return m_SDD; }
    T setSourceDetectorDistance(T SDD_cm)
    {
        m_SDD = std::max(std::abs(SDD_cm), T { 1 });
    }

    T scanFieldOfView() const { return m_FOV; }
    void setScanFieldOfView(T fov_cm)
    {
        m_FOV = std::max(std::abs(fov_cm), T { 1 });
    }

    T pitch() const { return m_pitch; }
    void setPitch(T p)
    {
        m_pitch = std::max(std::abs(p), T { 0.1 });
    }

    T startAngle() const { return m_startAngle; }
    void setStartAngle(T angle) { m_startAngle = angle; }
    T startAngleDeg() const { return m_startAngle * RAD_TO_DEG<T>(); }
    void setStartAngleDeg(T angle) { m_startAngle = angle * DEG_TO_RAD<T>(); }

    T stepAngle() const { return m_stepAngle; }
    void setStepAngle(T angle)
    {
        m_stepAngle = std::max(std::abs(angle), DEG_TO_RAD<T>() / 10;)
    }
    T stepAngleDeg() const { return m_stepAngle * RAD_TO_DEG<T>(); }
    void setStepAngleDeg(T angle) { setStepAngle(angle * DEG_TO_RAD<T>()); }

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

    CTSpiralBeamExposure<T> exposure(std::size_t i) const noexcept
    {
        constexpr auto pi2 = PI_VAL<T>() * 2;
        const T angle = i * m_stepAngle;
        const auto dz = m_pitch * m_collimation * angle / pi2;
        const auto directionZ = vectormath::subtract(m_stop, m_start);
        const auto direction = vectormath::normalized(directionZ);

        // finding normal vector to direction
        const auto normal_ind = vectormath::argmin3(direction);
        std::array<T, 3> normal = { 0, 0, 0 };
        normal[normal_ind] = 1;
        normal = vectormath::normalized(vectormath::cross(normal, direction));
        normal = vectormath::rotate(normal, direction, m_startAngle + angle);

        const auto beamdir = vectormath::cross(normal, direction);

        const std::array<std::array<T, 3>, 2> cosines = { normal, direction };

        auto pos = vectormath::add(vectormath::add(m_start, vectormath::scale(direction, dz)), vectormath::scale(beamdir, -m_SDD / 2));

        // position along cylinder axis
        const auto angx = std::atan(m_FOV / m_SDD);
        const auto angy = std::atan(T { 0.5 } * m_collimation / m_SDD);

        std::array<T, 2> angles = { angx, angy };

        CTSpiralBeamExposure<T> exp(pos, cosines, m_particlesPerExposure, m_weight, angles, m_specter);
        return exp;
    }

    T calibrationFactor() const
    {
        // todo
        return 1;
    }

protected:
    void tubeChanged()
    {
        const auto energies = m_tube.getEnergy();
        const auto weights = m_tube.getSpecter(energies, false);
        m_specter = SpecterDistribution(energies, weights);
    }

private:
    std::array<T, 3> m_start = { 0, 0, 0 };
    std::array<T, 3> m_stop = { 0, 0, 0 };
    T m_FOV = 50;
    T m_SDD = 100;
    T m_collimation = 1; // cm
    T m_pitch = 1;
    T m_startAngle = 0;
    T m_stepAngle = T { 0.018 }; // about a degree;
    T m_weight = 1;
    std::uint64_t m_particlesPerExposure = 100;
    Tube<T> m_tube;
    SpecterDistribution<T> m_specter;
};

}