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

#include "dxmc/beams/filters/bowtiefilter.hpp"
#include "dxmc/constants.hpp"
#include "dxmc/dxmcrandom.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/transportprogress.hpp"
#include "dxmc/vectormath.hpp"

#include <array>

namespace dxmc {

template <Floating T>
class CTDIBeamExposure {
public:
    CTDIBeamExposure(T angle, T SDD, std::uint64_t historiesPerExposure,
        const std::array<T, 2>& collimationAngles, const SpecterDistribution<T>* specter, const BowtieFilter<T>* bowtie, T weight = 1)
        : m_Nparticles(historiesPerExposure)
        , m_collimationAngles(collimationAngles)
        , m_specter(specter)
        , m_bowtieFilter(bowtie)
        , m_weight(weight)
    {
        m_dirCosineX = vectormath::rotate(m_dirCosineX, m_dirCosineY, angle);
        const auto beamdir = vectormath::cross(m_dirCosineX, m_dirCosineY);
        m_pos = vectormath::scale(beamdir, -SDD / 2);
        m_dir = vectormath::cross(m_dirCosineX, m_dirCosineY);
    }

    CTDIBeamExposure() = delete;

    const std::array<T, 3>& position() const { return m_pos; }

    std::uint64_t numberOfParticles() const
    {
        return m_Nparticles;
    }

    Particle<T> sampleParticle(RandomState& state) const noexcept
    {
        const auto angx = state.randomUniform(-m_collimationAngles[0], m_collimationAngles[0]);
        const auto angy = state.randomUniform(-m_collimationAngles[1], m_collimationAngles[1]);

        const auto sinx = std::sin(angx);
        const auto siny = std::sin(angy);
        const auto sinz = std::sqrt(1 - sinx * sinx - siny * siny);
        std::array pdir = {
            m_dirCosineX[0] * sinx + m_dirCosineY[0] * siny + m_dir[0] * sinz,
            m_dirCosineX[1] * sinx + m_dirCosineY[1] * siny + m_dir[1] * sinz,
            m_dirCosineX[2] * sinx + m_dirCosineY[2] * siny + m_dir[2] * sinz
        };

        const auto bowtie_weight = m_bowtieFilter->operator()(angx);

        Particle<T> p = {
            .pos = m_pos,
            .dir = pdir,
            .energy = m_specter->sampleValue(state),
            .weight = m_weight * bowtie_weight
        };
        return p;
    }

private:
    std::uint64_t m_Nparticles = 1;
    std::array<T, 3> m_pos = { 0, 0, 0 };
    std::array<T, 3> m_dir = { 0, 1, 0 };
    std::array<T, 3> m_dirCosineX = { 1, 0, 0 };
    std::array<T, 3> m_dirCosineY = { 0, 0, 1 };
    std::array<T, 2> m_collimationAngles = { 0, 0 };
    T m_weight = 1;
    const SpecterDistribution<T>* m_specter = nullptr;
    const BowtieFilter<T>* m_bowtieFilter = nullptr;
};

template <Floating T>
class CTDIBeam {
public:
    CTDIBeam(T angleStep, T SDD, const std::array<T, 2>& collimationAngles, std::uint64_t particlesPerExposure, const SpecterDistribution<T>& specter, const BowtieFilter<T>& bowtie, T weight = 1)
        : m_angleStep(angleStep)
        , m_sdd(SDD)
        , m_collimationAngles(collimationAngles)
        , m_particlesPerExposure(particlesPerExposure)
        , m_specter(specter)
        , m_bowtieFilter(bowtie)
        , m_weight(weight)
    {
    }

    std::uint64_t numberOfExposures() const
    {
        constexpr auto pi2 = PI_VAL<T>() * 2;
        return static_cast<std::uint64_t>(pi2 / m_angleStep);
    }

    std::uint64_t numberOfParticles() const { return numberOfExposures() * m_particlesPerExposure; }

    CTDIBeamExposure<T> exposure(std::size_t i) const noexcept
    {
        const auto angle = i * m_angleStep;
        return CTDIBeamExposure(angle, m_sdd, m_particlesPerExposure, m_collimationAngles, &m_specter, &m_bowtieFilter, m_weight);
    }

    T calibrationFactor(TransportProgress* progress = nullptr) const
    {
        return 1;
    }

private:
    T m_angleStep = 0;
    T m_sdd = 1;
    T m_weight = 1;
    std::array<T, 2> m_collimationAngles = { 0, 0 };
    std::uint64_t m_particlesPerExposure = 1;
    SpecterDistribution<T> m_specter;
    BowtieFilter<T> m_bowtieFilter;
};
}