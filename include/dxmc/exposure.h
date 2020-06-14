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

#include "dxmc/beamfilters.h"
#include "dxmc/dxmcrandom.h"
#include <array>
#include <vector>

struct Particle {
    double pos[3];
    double dir[3];
    double energy;
    double weight;
};

class Exposure {
public:
    Exposure(const BeamFilter* filter = nullptr, const SpecterDistribution* specter = nullptr, const HeelFilter* heelFilter = nullptr);

    void setPosition(double x, double y, double z);
    void setPosition(const double pos[3]);
    void setPosition(const std::array<double, 3>& pos) { m_position = pos; }
    void setPositionZ(const double posZ) { m_position[2] = posZ; }
    const std::array<double, 3>& position(void) const;

    void addPosition(const std::array<double, 3>& pos)
    {
        for (std::size_t i = 0; i < 3; ++i)
            m_position[i] += pos[i];
    }
    void subtractPosition(const std::array<double, 3>& pos)
    {
        for (std::size_t i = 0; i < 3; ++i)
            m_position[i] -= pos[i];
    }

    void setDirectionCosines(double x1, double x2, double x3, double y1, double y2, double y3);
    void setDirectionCosines(const double cosines[6]);
    void setDirectionCosines(const std::array<double, 6>& cosines);
    void setDirectionCosines(const std::array<double, 3>& cosinesX, const std::array<double, 3>& cosinesY);
    const std::array<double, 6>& directionCosines(void) const { return m_directionCosines; }

    const std::array<double, 3>& beamDirection(void) const { return m_beamDirection; }

    void setCollimationAngles(const double angles[2]);
    void setCollimationAngles(const std::array<double, 2>& angles);
    void setCollimationAngles(const double angleX, const double angleY);
    const std::array<double, 2>& collimationAngles(void) const { return m_collimationAngles; }
    double collimationAngleX() const;
    double collimationAngleY() const;

    void setBeamIntensityWeight(double weight) { m_beamIntensityWeight = weight; }
    double beamIntensityWeight(void) const { return m_beamIntensityWeight; }

    void setBeamFilter(const BeamFilter* filter);
    void setSpecterDistribution(const SpecterDistribution* specter);
    void setHeelFilter(const HeelFilter* filter);
    void setMonoenergeticPhotonEnergy(double energy);

    void setNumberOfHistories(std::size_t nHistories) { m_nHistories = nHistories; }
    std::size_t numberOfHistories(void) const { return m_nHistories; }

    void alignToDirectionCosines(const std::array<double, 6>& directionCosines);

    void sampleParticle(Particle& p, RandomState &state) const; // thread safe

protected:
    void normalizeDirectionCosines(void);
    void calculateBeamDirection(void);

private:
    std::array<double, 3> m_position;
    std::array<double, 6> m_directionCosines;
    std::array<double, 3> m_beamDirection;
    std::array<double, 2> m_collimationAngles;
    double m_beamIntensityWeight = 1.0;
    const BeamFilter* m_beamFilter = nullptr;
    const SpecterDistribution* m_specterDistribution = nullptr;
    const HeelFilter* m_heelFilter = nullptr;
    double m_monoenergeticPhotonEnergy = 0.0;
    std::size_t m_nHistories;
};
