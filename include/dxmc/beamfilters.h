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
#include "dxmc/tube.h"
#include "dxmc/world.h"

#include <array>
#include <string>
#include <vector>

/**
 * @brief Base class for filters on a photon beam
 * 
 * Base class for various filters on a photon beam based on angle of photom from beam direction along the first direction cosine vector of the source 
 */
class BeamFilter {
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
    /**
	 * @brief Get the name of the filter
	 * This is provided so we can keep track of many filters if needed. 
	 * @param  name of filter. Defaults to empty string.
	 * @return name of filter. Defaults to empty string.
	*/
    const std::string& filterName(void) const { return m_filterName; }
    /**
	 * @brief Set the name of the filter
	 * This is provided so we can keep track of many filters if needed. 
	 * @param  name of filter. Defaults to empty string.
	 * @return name of filter. Defaults to empty string.
	*/
    void setFilterName(const std::string& name) { m_filterName = name; }

private:
    std::string m_filterName = "";
};

/**
 * @brief Class for simple CT bowtie filter modeling. 
 * 
 * Filter to adjust photon weight based on a measured bowtie filter profile. Note that beamhardening of a specter is not modeled
 * and the bowtie filter simply adjust photon fluence according to a profile. 
 */
class BowTieFilter : public BeamFilter {
public:
    /**
	 * @brief Construct a new Bow Tie Filter object
	 * Construct a new Bow Tie Filter object
	 * 
	 * @param angles A vector of angles in radians for the fluence profile. The filter is assumed to be symmetrical and only absolute values of angles are used.
	 * @param weights Photon fluence for the corresponding angles. The fluence measurements are normalized such that the average photon weight is 1.0 when number of samples -> infinity. 
	 */
    BowTieFilter(const std::vector<double>& angles, const std::vector<double>& weights);

    /**
	 * @brief Construct a new Bow Tie Filter object
	 * Construct a new Bow Tie Filter object
	 * @param angleWeightsPairs A vector of angle in radians and photon fluence pairs.
	 */
    BowTieFilter(const std::vector<std::pair<double, double>>& angleWeightsPairs);

    /**
	 * @brief Construct a new Bow Tie Filter object from another BowTieFilter
	 * Construct a new Bow Tie Filter object
	 * @param other BowTiefilter
	 */
    BowTieFilter(const BowTieFilter& other);
    virtual ~BowTieFilter() = default;

    /**
	 * @brief Sample photon weight
	 * Sample photon weight modifier for a photon with direction angle in radians along the first Source cosine direction 
	 * @param angle Photon angle in radians 
	 * @return double Photon fluence weight modifier 
	 */
    double sampleIntensityWeight(const double angle) const override;

    /**
	 * @brief Filter data
	 * Returns a vector reference to angle in radians and normalized fluence weight pairs  
	 * @return const std::vector<std::pair<double, double>>& Vector of <angle, weight> pairs
	 */
    const std::vector<std::pair<double, double>>& data() const { return m_data; }

protected:
    /**
     * @brief Normalize filter
	 * Normalize filter such that expectation value of sampleIntensityWeight approches unity for large number of samples.
    */
    void normalizeData();

private:
    std::vector<std::pair<double, double>> m_data;
};

/**
 * @brief Filter for modeling of organ exposure control for Siemens CT scanners (XCare).
 * 
 * This filter modifies a particle's weight along the rotation angle in the same manner as Simens have iplementet on some CT scanner models.
 * The result is decreased fluence (photon weight) along the filters span angle and increased fluence outside the filter angle. The mean photon weight over all angles in a rotation is unity. 
 */
class XCareFilter : public BeamFilter {
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
	 * Highest photon weight outside span angle. This value is calculated such that expectation weight across all angles is 1.0
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
 * Filter for modeling of Heel effect of a x-ray tube. The filter do not model beam hardening of the Heel effect, only photon fluence effects. 
 * 
 */
class HeelFilter {
public:
    /**
     * @brief Constructs a new Heel filter
     * @param tube Tube to use when calculating Heel effect
     * @param heel_angle_span Span angle in radians to calculate Heel effect. This should be atleast as large as the beam collimation from the radiation source.
    */
    HeelFilter(const Tube& tube, const double heel_angle_span = 0.0);
    /**
     * @brief Update this filter with a new tube.
     * @param tube Tube to use when calculating Heel effect
     * @param heel_angle_span Span angle in radians to calculate Heel effect. This should be atleast as large as the beam collimation from the radiation source.
    */
    void update(const Tube& tube, const double heel_angle_span = 0.0);
    /**
     * @brief Sample weight from angle and photon energy
	 * The expectation value of weight from a large number of samples for a specific pfoton energy is normalized to unity. 
     * @param angle Angle of photon direction in anode cathode direction of a tube.
     * @param energy Ebergy in keV for photon. 
     * @return Weight of photon.
    */
    double sampleIntensityWeight(const double angle, const double energy) const;
    /**
     * @brief Size of energy bins used for interpolation
     * @return 
    */
    std::size_t energySize() const { return m_energySize; }
    /**
     * @brief Size of angle bins used for interpolation
     * @return 
    */
    std::size_t angleSize() const { return m_angleSize; }
    /**
     * @brief Buffer for angle and energy weights
	 * A buffer of angle size * energy size of weights used to interpolate across photon energies and angles.  
     * @return 
    */
    const std::vector<double>& weights() const { return m_weights; }
    /**
	 * @brief Get the name of the filter
	 * This is provided so we can keep track of many filters if needed. 
	 * @param  name of filter. Defaults to empty string.
	 * @return name of filter. Defaults to empty string.
	*/
    const std::string& filterName(void) const { return m_filterName; }
    /**
	 * @brief Set the name of the filter
	 * This is provided so we can keep track of many filters if needed. 
	 * @param  name of filter. Defaults to empty string.
	 * @return name of filter. Defaults to empty string.
	*/
    void setFilterName(const std::string& name) { m_filterName = name; }

private:
    double m_energyStep = 2.0;
    double m_energyStart = 20.0;
    std::size_t m_energySize = 65;

    double m_angleStep = 0.07; // about 4 deg step size
    double m_angleStart = 0.07;
    std::size_t m_angleSize = 5;
    std::vector<double> m_energies;
    std::vector<double> m_weights; //vector of m_angleSize*m_energySize weights
    std::string m_filterName = "";
};

/**
 * @brief Filter to adjust photon weights according to a tube current profile for CT examinations.
 * Filter to simulate automatic exposure control for CT examinations. This filter will match a given tube current profile to a denisty image to generate a lookup table of photon weights according to a ddensity profile. 
*/
class AECFilter {
public:
    /**
     * @brief Constructs a new AECFilter
     * @param densityImage Density array of size dimension[0]*dimension[1]*dimension[2]. 
     * @param spacing Spacing of each voxel for the density array.
     * @param dimensions Dimensions of the densityImage.
     * @param exposuremapping Exposure along the z (3rd) dimension, must be of size dimension[2].
    */
    AECFilter(const std::vector<double>& densityImage, const std::array<double, 3> spacing, const std::array<std::size_t, 3> dimensions, const std::vector<double>& exposuremapping);
    /**
     * @brief Constructs a new AECFilter
     * @param densityImage Density array of size dimension[0]*dimension[1]*dimension[2]. 
     * @param spacing Spacing of each voxel for the density array.
     * @param dimensions Dimensions of the densityImage.
     * @param exposuremapping Exposure along the z (3rd) dimension, must be of size dimension[2].
    */
    AECFilter(std::shared_ptr<std::vector<double>>& densityImage, const std::array<double, 3> spacing, const std::array<std::size_t, 3> dimensions, const std::vector<double>& exposuremapping);
    /**
     * @brief Constructs a new AECFilter directly from mass and weight arrays
     * @param mass Array of masses. Must be of same size as intensity array.
     * @param intensity Array of intensity or exposure corresponding to each mass in mass array. Must be of same size as mass array
    */
    AECFilter(const std::vector<double>& mass, const std::vector<double>& intensity);
    /**
     * @brief Sample weight of photon from a position in a world object.
     * Sample AEC weight of photon originating from position. This function assumes that updateFromWorld is called with the same World object the position is generated from. The expected value of weight from a large number of samples across all valid positions (the extent of the world object) is unity.
     * @param position Position to interpolate photon weight.
     * @return photon weight
    */
    double sampleIntensityWeight(const std::array<double, 3>& position) const;
    /**
     * @brief Update AECFilter from a World. 
     * This function will generate position lookup tables based on world, such as that photon weights can be sampled from a position within the world object.
     * @param world 
    */
    void updateFromWorld(const World& world);
    /**
     * @brief Returns if the AECFilter is valid and has been updated from a world.
     * @return 
    */
    bool isValid() const { return m_valid; }
    /**
     * @brief Masses used for interpolation of photon weights
     * @return Array of mass values
    */
    const std::vector<double>& mass() const { return m_mass; }
    /**
     * @brief Intensities used for interpolation of photon weights
     * @return Array of intensity values.
    */
    const std::vector<double>& massIntensity() const { return m_massIntensity; }
    /**
	 * @brief Get the name of the filter
	 * This is provided so we can keep track of many filters if needed. 
	 * @param  name of filter. Defaults to empty string.
	 * @return name of filter. Defaults to empty string.
	*/
    const std::string& filterName(void) const { return m_filterName; }
    /**
	 * @brief Set the name of the filter
	 * This is provided so we can keep track of many filters if needed. 
	 * @param  name of filter. Defaults to empty string.
	 * @return name of filter. Defaults to empty string.
	*/
    void setFilterName(const std::string& name) { m_filterName = name; }

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
    std::string m_filterName = "";
};
