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

#include "dxmc/beams/ctdibeam.hpp"
#include "dxmc/beams/ctspiralbeam.hpp"
#include "dxmc/beams/tube/tube.hpp"
#include "dxmc/dxmcrandom.hpp"
#include "dxmc/floating.hpp"
#include "dxmc/material/material.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/transport.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/worlditems/ctdiphantom.hpp"

#include <array>
#include <cmath>

namespace dxmc {

template <Floating T>
class CTSpiralDualEnergyBeam {

public:
    CTSpiralDualEnergyBeam(const std::array<T, 3>& start_pos = { 0, 0, 0 }, const std::array<T, 3>& stop_pos = { 0, 0, 0 })
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
        auto N_angles = static_cast<std::uint64_t>(total_rot_angle / m_stepAngle);
        return N_angles * 2;
    }

    std::uint64_t numberOfParticles() const { return numberOfExposures() * m_particlesPerExposure; }
    std::uint64_t numberOfParticlesPerExposure() const { return m_particlesPerExposure; }
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
    void setSourceDetectorDistance(T SDD_cm)
    {
        m_SDD = std::max(std::abs(SDD_cm), T { 1 });
    }

    T scanFieldOfViewB() const { return m_FOVB; }
    void setScanFieldOfViewB(T fov_cm)
    {
        m_FOVB = std::max(std::abs(fov_cm), T { 1 });
    }
    T scanFieldOfViewA() const { return m_FOVA; }
    void setScanFieldOfViewA(T fov_cm)
    {
        m_FOVA = std::max(std::abs(fov_cm), T { 1 });
    }

    T pitch() const { return m_pitch; }
    void setPitch(T p)
    {
        m_pitch = std::max(std::abs(p), T { 0.1 });
    }

    void setCTDIvol(T ctdi) { m_CTDIvol = ctdi; }
    T CTDIvol() const { return m_CTDIvol; }
    void setCTDIdiameter(T d) { m_CTDIdiameter = d; }
    T CTDIdiameter() const { return m_CTDIdiameter; }

    T startAngle() const { return m_startAngle; }
    void setStartAngle(T angle) { m_startAngle = angle; }
    T startAngleDeg() const { return m_startAngle * RAD_TO_DEG<T>(); }
    void setStartAngleDeg(T angle) { m_startAngle = angle * DEG_TO_RAD<T>(); }

    T tubeBoffsetAngle() const { return m_tubeBoffsetAngle; }
    void setTubeBoffsetAngle(T angle) { m_tubeBoffsetAngle = angle; }
    T tubeBoffsetAngleDeg() const { return m_tubeBoffsetAngle * RAD_TO_DEG<T>(); }
    void setTubeBoffsetAngleDeg(T angle) { m_tubeBoffsetAngle = angle * DEG_TO_RAD<T>(); }

    T stepAngle() const { return m_stepAngle; }
    void setStepAngle(T angle)
    {
        m_stepAngle = std::max(std::abs(angle), DEG_TO_RAD<T>() / 10);
    }
    T stepAngleDeg() const { return m_stepAngle * RAD_TO_DEG<T>(); }
    void setStepAngleDeg(T angle) { setStepAngle(angle * DEG_TO_RAD<T>()); }

    const Tube<T>& tubeA() const { return m_tubeA; }
    void setTubeA(const Tube<T>&& tube)
    {
        m_tubeA = tube;
        tubeChanged();
    }
    const Tube<T>& tubeB() const { return m_tubeB; }
    void setTubeB(const Tube<T>&& tube)
    {
        m_tubeB = tube;
        tubeChanged();
    }
    void setTubeAVoltage(T voltage)
    {
        m_tubeA.setVoltage(voltage);
        tubeChanged();
    }
    void setTubeBVoltage(T voltage)
    {
        m_tubeB.setVoltage(voltage);
        tubeChanged();
    }
    void setTubesAnodeAngle(T ang)
    {
        m_tubeA.setAnodeAngle(ang);
        m_tubeB.setAnodeAngle(ang);
        tubeChanged();
    }
    void setTubesAnodeAngleDeg(T ang)
    {
        m_tubeA.setAnodeAngleDeg(ang);
        m_tubeB.setAnodeAngleDeg(ang);
        tubeChanged();
    }
    void addTubeAFiltrationMaterial(std::size_t Z, T mm)
    {
        auto success = m_tubeA.addFiltrationMaterial(Z, mm);
        if (success)
            tubeChanged();
    }
    void addTubeBFiltrationMaterial(std::size_t Z, T mm)
    {
        auto success = m_tubeB.addFiltrationMaterial(Z, mm);
        if (success)
            tubeChanged();
    }
    void clearTubesFiltrationMaterials()
    {
        m_tubeA.clearFiltrationMaterials();
        m_tubeB.clearFiltrationMaterials();
        tubeChanged();
    }
    void setTubesEnergyResolution(T energyResolution)
    {
        m_tubeA.setEnergyResolution(energyResolution);
        m_tubeB.setEnergyResolution(energyResolution);
        tubeChanged();
    }

    CTSpiralBeamExposure<T> exposure(std::size_t dualExposureIndex) const noexcept
    {
        const std::size_t i = dualIndex / 2;
        const bool isTubeB = dualExposureIndex % 2 == 1;

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
        const auto rot_angle = isTubeB ? m_startAngle + angle + m_tubeBoffsetAngle : m_startAngle + angle;
        normal = vectormath::rotate(normal, direction, rot_angle);

        const auto beamdir = vectormath::cross(normal, direction);

        const std::array<std::array<T, 3>, 2> cosines = { normal, direction };

        auto pos = vectormath::add(vectormath::add(m_start, vectormath::scale(direction, dz)), vectormath::scale(beamdir, -m_SDD / 2));

        // position along cylinder axis
        const auto angx = isTubeB ? std::atan(m_FOVB / m_SDD) : std::atan(m_FOVA / m_SDD);
        const auto angy = std::atan(T { 0.5 } * m_collimation / m_SDD);

        std::array<T, 2> angles = { angx, angy };

        const auto weight = isTubeB ? m_weightB : m_weightA;
        const SpecterDstribution<T>* const specter = isTubeB ? &m_specterB : &m_specterA;
        CTSpiralBeamExposure<T> exp(pos, cosines, m_particlesPerExposure, weight, angles, specter);
        return exp;
    }

    T calibrationFactor(TransportProgress* progress = nullptr) const
    {
        // generating scoring world
        using Phantom = CTDIPhantom<T, 5, 1>;
        World<T, Phantom> world;
        world.reserveNumberOfItems(1);
        const auto& ctdi = world.template addItem<Phantom>({ m_CTDIdiameter });
        world.build();

        // generating CTDIbeam
        const auto angxA = std::atan(m_FOVA / m_SDD);
        const auto angy = std::atan(T { 0.5 } * m_collimation / m_SDD);
        const std::array<T, 2> collimationAnglesA = { angxA, angy };
        CTDIBeam<T> beamA(m_stepAngle, m_SDD, collimationAnglesA, m_particlesPerExposure, m_specterA, m_weightA);
        Transport transport;
        transport(world, beamA, progress, false);

        const auto angxB = std::atan(m_FOVB / m_SDD);
        const std::array<T, 2> collimationAnglesB = { angxB, angy };
        CTDIBeam<T> beamB(m_stepAngle, m_SDD, collimationAnglesB, m_particlesPerExposure, m_specterB, m_weightB);
        transport(world, beamB, progress, false);

        const T ctdiw_calc = (ctdi.centerDoseScored() + 2 * ctdi.pheriferyDoseScored()) * T { 10 } / (3 * m_collimation);

        const T ctdiw_beam = m_CTDIvol * m_pitch;
        return ctdiw_beam / ctdiw_calc;
    }

protected:
    void tubeChanged()
    {
        auto energiesA = m_tubeA.getEnergy();
        auto weightsA = m_tubeA.getSpecter(energiesA, false);
        const auto weightA = std::reduce(std::execution::par_unseq, weightsA.cbegin(), weightsA.cend(), T { 0 });
        m_specterA = SpecterDistribution(energiesA, weightsA);

        auto energiesB = m_tubeB.getEnergy();
        auto weightsB = m_tubeB.getSpecter(energiesB, false);
        const auto weightB = std::reduce(std::execution::par_unseq, weightsB.cbegin(), weightsB.cend(), T { 0 });
        m_specterB = SpecterDistribution(energiesB, weightsB);

        m_weightA = weightA / (weightA + weightB);
        m_weightB = weightB / (weightA + weightB);
    }

private:
    std::array<T, 3> m_start = { 0, 0, 0 };
    std::array<T, 3> m_stop = { 0, 0, 0 };
    T m_FOVA = 50;
    T m_FOVB = 50;
    T m_SDD = 100;
    T m_collimation = 1; // cm
    T m_pitch = 1;
    T m_startAngle = 0;
    T m_tubeBoffsetAngle = DEG_TO_RAD<T>() * 90;
    T m_stepAngle = T { 0.018 }; // about a degree;
    T m_weightA = 1;
    T m_weightB = 1;
    T m_CTDIvol = 1;
    T m_CTDIdiameter = 32;
    std::uint64_t m_particlesPerExposure = 100;
    Tube<T> m_tubeA;
    Tube<T> m_tubeB;
    SpecterDistribution<T> m_specterA;
    SpecterDistribution<T> m_specterB;
};

}
