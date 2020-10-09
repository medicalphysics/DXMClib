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

#include "dxmc/floating.h"
#include "dxmc/beamfilters.h"
#include "dxmc/dxmcrandom.h"
#include <array>
#include <vector>

namespace dxmc {

/**
 * @brief Simple struct to describe a photon 
*/
template <Floating T>
struct Particle {
    /**
     * @brief Position vector in three dimensions.
    */
    T pos[3];
    /**
     * @brief Direction vector in three dimension. This vector is threated as a normal vector.
    */
    T dir[3];
    /**
     * @brief Photon energy in keV.
    */
    T energy;
    /**
     * @brief Photon relative weight.
    */
    T weight;
};
/**
 * @brief Class to describe an exposure.
 * An exposure is used to describe a particle emitter where random distributions, position and, direction and collimation do not change for particles created.  
*/
template <Floating T>
class Exposure {
public:
    /**
     * @brief Initialize the exposure class
     * @param filter A BeamFilter to modify photon or particle weight according to the BeamFilter weight distribution. An example is a CT bowtie filter.
     * @param specter A SpecterDistribution to sample photon energies, if null all particles will be emittet with same energy. (see set setMonoenergeticPhotonEnergy)
     * @param heelFilter A HeelFilter to model the Heel effect of an x-ray tube. If null no Heel effect is modelled.
    */
    Exposure(const BeamFilter* filter = nullptr, const SpecterDistribution* specter = nullptr, const HeelFilter* heelFilter = nullptr);

    /**
     * @brief Set position of the exposure and also each particle starting point. Initial value is zero in all dimensions.
     * @param x 
     * @param y 
     * @param z 
    */
    void setPosition(T x, T y, T z);
    /**
     * @brief Set position of the exposure and also each particle starting point. Initial value is zero in all dimensions.
     * @param pos 
    */
    void setPosition(const T pos[3]);
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
    const std::array<T, 3>& position(void) const;

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
    void setDirectionCosines(T x1, T x2, T x3, T y1, T y2, T y3);
    /**
     * @brief Direction cosines describe a plane intersecting the exposure position. The plane is spanned by the orthonormal vectors X and Y. The cross product of X and Y gives the beam direction. The Beamfilter will modify the weights of a photon according to the X direction and the HeelFilter will modify photon weights by the Y direction.
     * @param cosines 
    */
    void setDirectionCosines(const T cosines[6]);
    /**
     * @brief Direction cosines describe a plane intersecting the exposure position. The plane is spanned by the orthonormal vectors X and Y. The cross product of X and Y gives the beam direction. The Beamfilter will modify the weights of a photon according to the X direction and the HeelFilter will modify photon weights by the Y direction.
     * @param cosines 
    */
    void setDirectionCosines(const std::array<T, 6>& cosines);
    /**
     * @brief Direction cosines describe a plane intersecting the exposure position. The plane is spanned by the orthonormal vectors X and Y. The cross product of X and Y gives the beam direction. The Beamfilter will modify the weights of a photon according to the X direction and the HeelFilter will modify photon weights by the Y direction.
     * @param cosinesX 
     * @param cosinesY 
    */
    void setDirectionCosines(const std::array<T, 3>& cosinesX, const std::array<T, 3>& cosinesY);
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

    void setCollimationAngles(const T angles[2]);
    void setCollimationAngles(const std::array<T, 2>& angles);
    void setCollimationAngles(const T angleX, const T angleY);
    const std::array<T, 2>& collimationAngles(void) const { return m_collimationAngles; }
    T collimationAngleX() const;
    T collimationAngleY() const;

    void setBeamIntensityWeight(T weight) { m_beamIntensityWeight = weight; }
    T beamIntensityWeight(void) const { return m_beamIntensityWeight; }

    void setBeamFilter(const BeamFilter* filter);
    void setSpecterDistribution(const SpecterDistribution* specter);
    void setHeelFilter(const HeelFilter* filter);
    void setMonoenergeticPhotonEnergy(T energy);

    void setNumberOfHistories(std::size_t nHistories) { m_nHistories = nHistories; }
    std::size_t numberOfHistories(void) const { return m_nHistories; }

    void alignToDirectionCosines(const std::array<T, 6>& directionCosines);

    void sampleParticle(Particle<T>& p, RandomState& state) const; // thread safe

protected:
    void normalizeDirectionCosines(void);
    void calculateBeamDirection(void);

private:
    std::array<T, 3> m_position;
    std::array<T, 6> m_directionCosines;
    std::array<T, 3> m_beamDirection;
    std::array<T, 2> m_collimationAngles;
    T m_beamIntensityWeight { 1.0 };
    const BeamFilter* m_beamFilter = nullptr;
    const SpecterDistribution* m_specterDistribution = nullptr;
    const HeelFilter* m_heelFilter = nullptr;
    T m_monoenergeticPhotonEnergy { 0.0 };
    std::size_t m_nHistories;
};
Exposure<double>;
};