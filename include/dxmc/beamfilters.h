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

#pragma once
#include "dxmc/world.h"
#include "dxmc/tube.h"

#include <vector>
#include <array>

/**
 * @brief Base class for filters on a photon beam
 * 
 * Base class for various filters on a photon beam based on angle of photom from beam direction along the first direction cosine vector of the source 
 */
class BeamFilter
{
public:
	/**
	 * @brief Construct a new Beam Filter object
	 * 
	 */
	BeamFilter() {};
	/**
	 * @brief Sample intensity weigh
	 * Sample the photon weight from this filter. The mean of photon weights sampled is 1.0 for N samples when N->infinity.  
	 * @param angle Angle in radians.
	 * @return double Photon weigh
	 */
	virtual double sampleIntensityWeight(const double angle) const = 0;

};


class BowTieFilter : public BeamFilter
{
public:
	BowTieFilter(const std::vector<double>& angles, const std::vector<double>& weights);
    BowTieFilter(const std::vector<std::pair<double, double>>& angleWeightsPairs);
	BowTieFilter(const BowTieFilter& other);
	virtual ~BowTieFilter() = default;
	double sampleIntensityWeight(const double angle) const override;
	const std::vector<std::pair<double, double>>& data() const { return m_data; }
protected:
	void normalizeData();
private:
    std::vector<std::pair<double, double>> m_data;
};


class XCareFilter :public BeamFilter
{
public:
	XCareFilter();
	virtual ~XCareFilter() = default;
	double filterAngle() const { return m_filterAngle; }
	double filterAngleDeg() const;
	void setFilterAngle(double angle);
	void setFilterAngleDeg(double angle);

	double spanAngle() const { return m_spanAngle; }
	double spanAngleDeg() const;
	void setSpanAngle(double angle);
	void setSpanAngleDeg(double angle);

	double rampAngle() const { return m_rampAngle; }
	double rampAngleDeg() const;
	void setRampAngle(double angle);
	void setRampAngleDeg(double angle);

	double lowWeight() const { return m_lowWeight; }
	void setLowWeight(double weight);
	double highWeight() const;

	double sampleIntensityWeight(const double angle) const override;

private:
	double m_filterAngle;
	double m_spanAngle;
	double m_rampAngle;
	double m_lowWeight;
};

class HeelFilter
{
public:
	HeelFilter(const Tube& tube, const double heel_angle_span = 0.0);
	void update(const Tube& tube, const double heel_angle_span = 0.0);
	double sampleIntensityWeight(const double angle, const double energy) const;

	std::size_t energySize() const { return m_energySize; }
	std::size_t angleSize() const { return m_angleSize; }
	const std::vector<double>& weights() const { return m_weights; }
private:
	double m_energyStep = 2.0;
	double m_energyStart = 20.0;
	std::size_t m_energySize = 65;

	double m_angleStep = 0.07; // about 4 deg step size
	double m_angleStart = 0.07;
	std::size_t m_angleSize = 5;
	std::vector<double> m_energies;
	std::vector<double> m_weights; //vector of m_angleSize*m_energySize weights
};


class AECFilter
{
public:
	AECFilter(const std::vector<double>& densityImage, const std::array<double, 3> spacing, const std::array<std::size_t, 3> dimensions, const std::vector<double>& exposuremapping);
	AECFilter(std::shared_ptr<std::vector<double>>& densityImage, const std::array<double, 3> spacing, const std::array<std::size_t, 3> dimensions, const std::vector<double>& exposuremapping);
	AECFilter(const std::vector<double>& mass, const std::vector<double>& intensity);
	double sampleIntensityWeight(const std::array<double, 3>& position) const;
	void updateFromWorld(const World& world);
	bool isValid() const { return m_valid; }
	const std::vector<double>& mass() const { return m_mass; }
	const std::vector<double>& massIntensity() const { return m_massIntensity; }

protected:
	void generateMassWeightMap(std::vector<double>::const_iterator densBeg, std::vector<double>::const_iterator densEnd, const std::array<double, 3> spacing, const std::array<std::size_t, 3> dimensions, const std::vector<double>& exposuremapping);
	void generatePositionWeightMap(std::vector<double>::const_iterator densBeg, std::vector<double>::const_iterator densEnd, const std::array<double, 3> spacing, const std::array<std::size_t, 3> dimensions, const std::array<double, 3>& origin);
	double interpolateMassIntensity(double mass) const;

private:
	bool m_valid = false;
	std::vector<double> m_mass;
	std::vector<double> m_massIntensity;
	double m_positionStep = 0.0;
	double m_positionMin = 0.0;
	double m_positionMax = 0.0;
	std::vector<double> m_positionIntensity;
};
