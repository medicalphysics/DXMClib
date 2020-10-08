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

Copyright 2019 Erlend Andersen
*/

#include "dxmc/exposure.h"
#include "dxmc/vectormath.h"

namespace dxmc {
template <Floating T>
Exposure<T>::Exposure(const BeamFilter* filter, const SpecterDistribution* specter, const HeelFilter* heelFilter)
{
    for (std::size_t i = 0; i < 3; ++i) {
        m_position[i] = T { 0.0 };
        m_directionCosines[i] = T { 0.0 };
        m_directionCosines[i + 3] = T { 0.0 };
    }
    m_directionCosines[0] = T { 1.0 };
    m_directionCosines[5] = T { 1.0 };
    calculateBeamDirection();
    m_collimationAngles[0] = T { 0.35 }; // about 20 deg
    m_collimationAngles[1] = T { 0.35 };
    m_beamIntensityWeight = T { 1.0 };
    m_specterDistribution = specter;
    m_beamFilter = filter;
    m_heelFilter = heelFilter;
}
template <Floating T>
void Exposure<T>::setPosition(T x, T y, T z)
{
    m_position[0] = x;
    m_position[1] = y;
    m_position[2] = z;
}
template <Floating T>
void Exposure<T>::setPosition(const T pos[3])
{
    for (std::size_t i = 0; i < 3; i++) {
        m_position[i] = pos[i];
    }
}
template <Floating T>
const std::array<T, 3>& Exposure<T>::position(void) const
{
    return m_position;
}
template <Floating T>
void Exposure<T>::setDirectionCosines(T x1, T x2, T x3, T y1, T y2, T y3)
{
    m_directionCosines[0] = x1;
    m_directionCosines[1] = x2;
    m_directionCosines[2] = x3;
    m_directionCosines[3] = y1;
    m_directionCosines[4] = y2;
    m_directionCosines[5] = y3;
    normalizeDirectionCosines();
}

template <Floating T>
void Exposure<T>::setDirectionCosines(const T cosines[6])
{
    for (std::size_t i = 0; i < 6; i++) {
        m_directionCosines[i] = cosines[i];
    }
    normalizeDirectionCosines();
}

template <Floating T>
void Exposure<T>::setDirectionCosines(const std::array<T, 6>& cosines)
{
    m_directionCosines = cosines;
    normalizeDirectionCosines();
}

template <Floating T>
void Exposure<T>::setDirectionCosines(const std::array<T, 3>& cosinesX, const std::array<T, 3>& cosinesY)
{
    for (std::size_t i = 0; i < 3; i++) {
        m_directionCosines[i] = cosinesX[i];
        m_directionCosines[i + 3] = cosinesY[i];
    }
    normalizeDirectionCosines();
}

template <Floating T>
void Exposure<T>::normalizeDirectionCosines()
{
    T sumx { 0 };
    T sumy { 0 };
    for (std::size_t i = 0; i < 3; i++) {
        sumx += m_directionCosines[i] * m_directionCosines[i];
        sumy += m_directionCosines[i + 3] * m_directionCosines[i + 3];
    }

    sumx = T { 1.0 } / std::sqrt(sumx);
    sumy = T { 1.0 } / std::sqrt(sumy);

    for (std::size_t i = 0; i < 3; i++) {
        m_directionCosines[i] *= sumx;
        m_directionCosines[i + 3] *= sumy;
    }

    calculateBeamDirection();
}

template <Floating T>
void Exposure<T>::calculateBeamDirection()
{
    vectormath::cross(m_directionCosines.data(), m_beamDirection.data());
}

template <Floating T>
void Exposure<T>::setCollimationAngles(const T angles[2])
{
    m_collimationAngles[0] = angles[0];
    m_collimationAngles[1] = angles[1];
}

template <Floating T>
void Exposure<T>::setCollimationAngles(const std::array<T, 2>& angles)
{
    m_collimationAngles = angles;
}

template <Floating T>
void Exposure<T>::setCollimationAngles(const T angleX, const T angleY)
{
    m_collimationAngles[0] = angleX;
    m_collimationAngles[1] = angleY;
}

template <Floating T>
T Exposure<T>::collimationAngleX() const
{
    return m_collimationAngles[0];
}

template <Floating T>
T Exposure<T>::collimationAngleY() const
{
    return m_collimationAngles[1];
}

template <Floating T>
void Exposure<T>::setBeamFilter(const BeamFilter* filter)
{
    m_beamFilter = filter;
}

template <Floating T>
void Exposure<T>::setSpecterDistribution(const SpecterDistribution* specter)
{
    m_specterDistribution = specter;
}

template <Floating T>
void Exposure<T>::setHeelFilter(const HeelFilter* filter)
{
    m_heelFilter = filter;
}

template <Floating T>
void Exposure<T>::setMonoenergeticPhotonEnergy(T energy)
{
    if (energy > T { 500.0 })
        energy = T { 500.0 };
    if (energy < T { 0.0 })
        energy = T { 0.0 };
    m_monoenergeticPhotonEnergy = energy;
}

template <Floating T>
void Exposure<T>::alignToDirectionCosines(const std::array<T, 6>& directionCosines)
// aligning coordinates to new basis given by cosines
{
    const T* b1 = directionCosines.data();
    const T* b2 = &b1[3];
    T b3[3];
    vectormath::cross(b1, b2, b3);
    vectormath::changeBasisInverse(b1, b2, b3, m_position.data());
    vectormath::changeBasisInverse(b1, b2, b3, m_directionCosines.data());
    vectormath::changeBasisInverse(b1, b2, b3, &m_directionCosines[3]);
    vectormath::changeBasisInverse(b1, b2, b3, m_beamDirection.data());
}

template <Floating T>
void Exposure<T>::sampleParticle(Particle<T>& p, RandomState& state) const
{

    p.pos[0] = m_position[0];
    p.pos[1] = m_position[1];
    p.pos[2] = m_position[2];

    // particle direction
    const T theta = state.randomUniform(-m_collimationAngles[0] / 2.0, m_collimationAngles[0] / 2.0);
    const T phi = state.randomUniform(-m_collimationAngles[1] / 2.0, m_collimationAngles[1] / 2.0);
    const T sintheta = std::sin(theta);
    const T sinphi = std::sin(phi);
    const T sin2theta = sintheta * sintheta;
    const T sin2phi = sinphi * sinphi;
    const T norm = 1.0 / std::sqrt(1.0 + sin2phi + sin2theta);
    for (std::size_t i = 0; i < 3; i++) {
        p.dir[i] = norm * (m_beamDirection[i] + sintheta * m_directionCosines[i] + sinphi * m_directionCosines[i + 3]);
    }

    if (m_specterDistribution) {
        p.energy = m_specterDistribution->sampleValue(state);
    } else {
        p.energy = m_monoenergeticPhotonEnergy;
    }

    p.weight = m_beamIntensityWeight;
    if (m_beamFilter) {
        p.weight *= m_beamFilter->sampleIntensityWeight(theta);
    }
    if (m_heelFilter) {
        p.weight *= m_heelFilter->sampleIntensityWeight(phi, p.energy);
    }
}
}