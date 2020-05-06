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
	 * @return double Photon weight
	 */
	virtual double sampleIntensityWeight(const double angle) const = 0;

};

/**
 * @brief Class for simple CT bowtie filter modelling. 
 * 
 * Filter to adjust photon weight based on a measured bowtie filter profile. Note that beamhardening of a specter is not modelled
 * and the bowtie filter simply adjust photon fluence according to a profile. 
 */
class BowTieFilter : public BeamFilter
{
public:
	/**
	 * @brief Construct a new Bow Tie Filter object
	 * 
	 * 
	 * @param angles A vector of angles in radians for the fluence profile. The filter is assumed to be symmetrical and only absolute values of angles are used.
	 * @param weights Photon fluence for the corresponding angles. The fluence measurements are normalized such that the average photon weight is 1.0 when number of samples -> infinity. 
	 */
	BowTieFilter(const std::vector<double>& angles, const std::vector<double>& weights);

	/**
	 * @brief Construct a new Bow Tie Filter object
	 * 
	 * @param angleWeightsPairs A vector of angle in radians and photon fluence pairs.
	 */
    BowTieFilter(const std::vector<std::pair<double, double>>& angleWeightsPairs);

	/**
	 * @brief Construct a new Bow Tie Filter object from another BowTieFilter
	 * 
	 * @param other BowTiefilter
	 */
	BowTieFilter(const BowTieFilter& other);
	virtual ~BowTieFilter() = default;

	/**
	 * @brief Sample photon weight modifier for a photon with direction angle in radians along the first Source cosine direction 
	 * 
	 * @param angle Photon angle in radians 
	 * @return double Photon fluence weight modifier 
	 */
	double sampleIntensityWeight(const double angle) const override;

	/**
	 * @brief Filter data
	 * 
	 * @return const std::vector<std::pair<double, double>>& A vector reference to angle in radians and normalized fluence weight pairs  
	 */
	const std::vector<std::pair<double, double>>& data() const { return m_data; }
protected:
	void normalizeData();
private:
    std::vector<std::pair<double, double>> m_data;
};

/**
 * @brief Filter for modelling of organ exposure control for Siemens CT scanners (XCare).
 * 
 * This filter modifies a particle's weight along the rotation angle in the same manner as Simens have iplementet on some CT scanner models.
 * The result is decreased fluence (photon weight) along the filters span angle and increased fluence outside the filter angle. The mean photon weight over all angles in a rotation is unity. 
 */
class XCareFilter :public BeamFilter
{
public:
	/**
	 * @brief Construct a new XCareFilter object
	 * Constructs new XCareFilter with default values
	 */
	XCareFilter();
	virtual ~XCareFilter() = default;

	/**
	 * @brief Center filter angle
	 * 
	 * Center angle for the filter where photon weight is reduced
	 * @return double Angle in radians
	 */
	double filterAngle() const { return m_filterAngle; }
	
	/**
	 * @brief Center filter angle
	 * 
	 * Center angle for the filter where photon weight is reduced
	 * @return double Angle in degrees
	 */
	double filterAngleDeg() const;

	/**
	 * @brief Set the Filter Angle
	 * 
	 * Center angle for the filter where photon weight is reduced
	 * @param angle New filter angle in radians 
	 */
	void setFilterAngle(double angle);
	
	/**
	 * @brief Set the Filter Angle
	 * 
	 * Center angle for the filter where photon weight is reduced
	 * @param angle New filter angle in degrees
	 */
	void setFilterAngleDeg(double angle);

	/**
	 * @brief Span angle of the filter
	 * 
	 * The span where photon weight is reduced, this angle span includes the ramp angle. Valid values are 5 degrees > span angle > PI. 
	 * @return double Span angle in radians
	 */
	double spanAngle() const { return m_spanAngle; }
	
	/**
	 * @brief Span angle of the filter
	 * 
	 * The span where photon weight is reduced, this angle span includes the ramp angle. Valid values are 5 degrees > span angle > PI. 
	 * @return double Span angle in degrees
	 */
	double spanAngleDeg() const;
	
	/**
	 * @brief Set span angle of the filter
	 * 
	 * The span where photon weight is reduced, this angle span includes the ramp angle. Valid values are 5 degrees > span angle > PI. 
	 * @return double Span angle in radians
	 */
	void setSpanAngle(double angle);

	/**
	 * @brief Set span angle of the filter
	 * 
	 * The span where photon weight is reduced, this angle span includes the ramp angle. Valid values are 5 degrees > span angle > PI. 
	 * @return double Span angle in degrees
	 */
	void setSpanAngleDeg(double angle);

	/**
	 * @brief Ramp angle 
	 * 
	 * The ramp angle span for ramping down the photon weight. The ramp angle are included in the span angle and must be less than spanAngle/2. 
	 * @return double Ramp angle in radians
	 */
	double rampAngle() const { return m_rampAngle; }

	/**
	 * @brief Ramp angle 
	 * 
	 * The ramp angle span for ramping down the photon weight. The ramp angle are included in the span angle and must be less than spanAngle/2. 
	 * @return double Ramp angle in degrees
	 */
	double rampAngleDeg() const;

	/**
	 * @brief Set ramp angle
	 * 
	 * The ramp angle span for ramping down the photon weight. The ramp angle are included in the span angle and must be less than spanAngle/2. 
	 * @param angle ramp angle in radians 
	 */
	void setRampAngle(double angle);

	/**
	 * @brief Set ramp angle
	 * 
	 * The ramp angle span for ramping down the photon weight. The ramp angle are included in the span angle and must be less than spanAngle/2. 
	 * @param angle ramp angle in degrees 
	 */
	void setRampAngleDeg(double angle);

	/**
	 * @brief Photon weight in the filters span angle excluding angle ramps 
	 * 
	 * Photon weight in the span angle excluding ramp angles. Must be in 0.0 < lowWeight <= 1.0.
	 * @return double 
	 */
	double lowWeight() const { return m_lowWeight; }

	/**
	 * @brief Set photon weight in the filters span angle excluding angle ramps 
	 * 
	 * Photon weight in the span angle excluding ramp angles. Must be in 0.0 < lowWeight <= 1.0.
	 * @param weight 
	 */
	void setLowWeight(double weight);

	/**
	 * @brief Photon weight for angles outside span angle
	 * 
	 * Highets photon weight outside span angle. This value is calculated such that mean weight from all angles is 1.0
	 * @return double Photon hight weight
	 */
	double highWeight() const;

	/**
	 * @brief Sample photon weight
	 * 
	 * Sample photon weight from rotation angle. Mean weight from all angles is 1.0
	 * 
	 * @param angle Rotation angle in radians
	 * @return double Photon weight modifier
	 */
	double sampleIntensityWeight(const double angle) const override;

private:
	double m_filterAngle;
	double m_spanAngle;
	double m_rampAngle;
	double m_lowWeight;
};

/**
 * @brief Filter to model Heel effect of a tube.
 * 
 * Filter for modelling of Heel effect of a x-ray tube. The filter do not model beam hardening of the Heel effect, only photon fluence effects. 
 * 
 */
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
