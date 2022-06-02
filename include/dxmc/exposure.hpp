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

#pragma once // include guard

#include "dxmc/beamfilters.hpp"
#include "dxmc/dxmcrandom.hpp"
#include "dxmc/floating.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"

#include <array>
#include <vector>

namespace dxmc {

/**
 * @brief Class to describe an exposure.
 * An exposure is used to describe a particle emitter where random distributions, position and, direction and collimation do not change for particles created.
 */
template <Floating T = double>
class Exposure {

private:
    std::array<T, 3> m_position;
    std::array<T, 6> m_directionCosines;
    std::array<T, 3> m_beamDirection;
    std::array<T, 4> m_collimationAngles; // x0 x1 y0 y1
    T m_beamIntensityWeight;
    const BeamFilter<T>* m_beamFilter = nullptr;
    const SpecterDistribution<T>* m_specterDistribution = nullptr;
    const HeelFilter<T>* m_heelFilter = nullptr;
    T m_monoenergeticPhotonEnergy { 0 };
    std::uint64_t m_nHistories;

protected:
    void normalizeDirectionCosines(void)
    {
        vectormath::normalize(&m_directionCosines[0]);
        vectormath::normalize(&m_directionCosines[3]);
        calculateBeamDirection();
    }

    void calculateBeamDirection(void)
    {
        vectormath::cross(m_directionCosines.data(), m_beamDirection.data());
    }

public:
    /**
     * @brief Initialize the exposure class
     * @param filter A BeamFilter to modify photon or particle weight according to the BeamFilter weight distribution. An example is a CT bowtie filter.
     * @param specter A SpecterDistribution to sample photon energies, if null all particles will be emittet with same energy. (see set setMonoenergeticPhotonEnergy)
     * @param heelFilter A HeelFilter to model the Heel effect of an x-ray tube. If null no Heel effect is modelled.
     */
    Exposure(const std::array<T, 3>& position,
        const std::array<T, 6>& directionCosines,
        const std::array<T, 4>& collimationAngles,
        std::uint64_t nHistories = 1000,
        T beamIntensityWeight = 1,
        const SpecterDistribution<T>* specterDistribution = nullptr,
        const HeelFilter<T>* heelFilter = nullptr,
        const BeamFilter<T>* filter = nullptr)
        : m_position(position)
        , m_directionCosines(directionCosines)
        , m_collimationAngles(collimationAngles)
        , m_nHistories(nHistories)
        , m_beamIntensityWeight(beamIntensityWeight)
        , m_specterDistribution(specterDistribution)
        , m_heelFilter(heelFilter)
        , m_beamFilter(filter)
    {
        normalizeDirectionCosines();
    }

    Exposure(const std::array<T, 3>& position,
        const std::array<T, 6>& directionCosines,
        const std::array<T, 2>& collimationAngles,
        std::uint64_t nHistories = 1000,
        T beamIntensityWeight = 1,
        const SpecterDistribution<T>* specterDistribution = nullptr,
        const HeelFilter<T>* heelFilter = nullptr,
        const BeamFilter<T>* filter = nullptr)
        : m_position(position)
        , m_directionCosines(directionCosines)
        , m_nHistories(nHistories)
        , m_beamIntensityWeight(beamIntensityWeight)
        , m_specterDistribution(specterDistribution)
        , m_heelFilter(heelFilter)
        , m_beamFilter(filter)
    {
        setCollimationAngles(collimationAngles);
        normalizeDirectionCosines();
    }

    /**
     * @brief Set position of the exposure and also each particle starting point. Initial value is zero in all dimensions.
     * @param x
     * @param y
     * @param z
     */
    void setPosition(T x, T y, T z)
    {
        m_position[0] = x;
        m_position[1] = y;
        m_position[2] = z;
    }
    /**
     * @brief Set position of the exposure and also each particle starting point. Initial value is zero in all dimensions.
     * @param pos
     */
    void setPosition(const T pos[3])
    {
        for (std::size_t i = 0; i < 3; i++) {
            m_position[i] = pos[i];
        }
    }
    /**
     * @brief Set position of the exposure and also each particle starting point. Initial value is zero in all dimensions.
     * @param pos
     */
    void setPosition(const std::array<T, 3>& pos) { m_position = pos; }
    /**
     * @brief Set position of the exposure in the Z or third dimension.
     * @param posZ
     */
    void setPositionZ(const T posZ) { m_position[2] = posZ; }
    /**
     * @brief Position of the exposure and also each particle starting point.
     * @param
     * @return
     */
    const std::array<T, 3>& position(void) const { return m_position; }

    /**
     * @brief Add to current exposure position
     * @param pos
     */
    void addPosition(const std::array<T, 3>& pos)
    {
        for (std::size_t i = 0; i < 3; ++i)
            m_position[i] += pos[i];
    }
    /**
     * @brief Subtract from the current position
     * @param pos
     */
    void subtractPosition(const std::array<T, 3>& pos)
    {
        for (std::size_t i = 0; i < 3; ++i)
            m_position[i] -= pos[i];
    }

    /**
     * @brief Set the direction cosines for this exposure.
     * Direction cosines describe a plane intersecting the exposure position. The plane is spanned by the orthonormal vectors X and Y. The cross product of X and Y gives the beam direction. The Beamfilter will modify the weights of a photon according to the X direction and the HeelFilter will modify photon weights by the Y direction.
     * @param x1
     * @param x2
     * @param x3
     * @param y1
     * @param y2
     * @param y3
     */
    void setDirectionCosines(T x1, T x2, T x3, T y1, T y2, T y3)
    {
        m_directionCosines[0] = x1;
        m_directionCosines[1] = x2;
        m_directionCosines[2] = x3;
        m_directionCosines[3] = y1;
        m_directionCosines[4] = y2;
        m_directionCosines[5] = y3;
        normalizeDirectionCosines();
    }
    /**
     * @brief Direction cosines describe a plane intersecting the exposure position. The plane is spanned by the orthonormal vectors X and Y. The cross product of X and Y gives the beam direction. The Beamfilter will modify the weights of a photon according to the X direction and the HeelFilter will modify photon weights by the Y direction.
     * @param cosines
     */
    void setDirectionCosines(const T cosines[6])
    {
        for (std::size_t i = 0; i < 6; i++) {
            m_directionCosines[i] = cosines[i];
        }
        normalizeDirectionCosines();
    }

    /**
     * @brief Direction cosines describe a plane intersecting the exposure position. The plane is spanned by the orthonormal vectors X and Y. The cross product of X and Y gives the beam direction. The Beamfilter will modify the weights of a photon according to the X direction and the HeelFilter will modify photon weights by the Y direction.
     * @param cosines
     */
    void setDirectionCosines(const std::array<T, 6>& cosines)
    {
        m_directionCosines = cosines;
        normalizeDirectionCosines();
    }

    /**
     * @brief Direction cosines describe a plane intersecting the exposure position. The plane is spanned by the orthonormal vectors X and Y. The cross product of X and Y gives the beam direction. The Beamfilter will modify the weights of a photon according to the X direction and the HeelFilter will modify photon weights by the Y direction.
     * @param cosinesX
     * @param cosinesY
     */
    void setDirectionCosines(const std::array<T, 3>& cosinesX, const std::array<T, 3>& cosinesY)
    {
        for (std::size_t i = 0; i < 3; i++) {
            m_directionCosines[i] = cosinesX[i];
            m_directionCosines[i + 3] = cosinesY[i];
        }
        normalizeDirectionCosines();
    }
    /**
     * @brief Direction cosines for the exposure
     * @return An array of six values where the first three is the X directional cosine and the last three is the Y directional cosine
     */
    const std::array<T, 6>& directionCosines(void) const { return m_directionCosines; }

    /**
     * @brief Beam direction, given by X direction cosine cross Y direction cosine
     * @return
     */
    const std::array<T, 3>& beamDirection(void) const { return m_beamDirection; }

    void setCollimationAngles(const std::array<T, 4>& angles) { m_collimationAngles = angles; }
    void setCollimationAngles(const std::array<T, 2>& angles) { setCollimationAngles(angles[0], angles[1]); }
    void setCollimationAngles(const T angleX, const T angleY)
    {
        m_collimationAngles[0] = -angleX / 2;
        m_collimationAngles[1] = angleX / 2;
        m_collimationAngles[2] = -angleY / 2;
        m_collimationAngles[3] = angleY / 2;
    }

    const std::array<T, 4>& collimationAngles(void) const { return m_collimationAngles; }
    T collimationAngleX() const { return m_collimationAngles[1] - m_collimationAngles[0]; }
    T collimationAngleY() const { return m_collimationAngles[3] - m_collimationAngles[2]; }

    void setBeamIntensityWeight(T weight) { m_beamIntensityWeight = weight; }
    T beamIntensityWeight(void) const { return m_beamIntensityWeight; }

    void setBeamFilter(const BeamFilter<T>* filter) { m_beamFilter = filter; }
    void setSpecterDistribution(const SpecterDistribution<T>* specter) { m_specterDistribution = specter; }
    void setHeelFilter(const HeelFilter<T>* filter) { m_heelFilter = filter; }
    void setMonoenergeticPhotonEnergy(T energy)
    {
        if (energy > T { 500.0 })
            energy = T { 500.0 };
        if (energy < T { 0.0 })
            energy = T { 0.0 };
        m_monoenergeticPhotonEnergy = energy;
    }

    void setNumberOfHistories(std::size_t nHistories) { m_nHistories = nHistories; }
    std::size_t numberOfHistories(void) const { return m_nHistories; }

    void alignToDirectionCosines(const std::array<T, 6>& directionCosines) noexcept
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

    Particle<T> sampleParticle(RandomState& state) const noexcept // thread safe
    {
        // particle direction
        const T theta = state.randomUniform(m_collimationAngles[0], m_collimationAngles[1]);
        const T phi = state.randomUniform(m_collimationAngles[2], m_collimationAngles[3]);

        Particle<T> p { .pos = m_position, .dir = m_beamDirection, .weight = m_beamIntensityWeight };

        vectormath::rotate(p.dir.data(), &m_directionCosines[3], theta);
        vectormath::rotate(p.dir.data(), &m_directionCosines[0], phi);

        if (m_specterDistribution) {
            p.energy = m_specterDistribution->sampleValue(state);
        } else {
            p.energy = m_monoenergeticPhotonEnergy;
        }

        if (m_beamFilter) {
            p.weight *= m_beamFilter->sampleIntensityWeight(theta);
        }
        if (m_heelFilter) {
            p.weight *= m_heelFilter->sampleIntensityWeight(phi, p.energy);
        }
        return p;
    }
};
}
