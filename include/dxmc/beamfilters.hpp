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
#include "dxmc/constants.hpp"
#include "dxmc/floating.hpp"
#include "dxmc/interpolation.hpp"
#include "dxmc/tube.hpp"
#include "dxmc/world.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <execution>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

namespace dxmc {

/**
 * @brief Base class for filters on a photon beam
 * 
 * Base class for various filters on a photon beam based on angle of photom from beam direction along the first direction cosine vector of the source 
 */
template <Floating T = double>
class BeamFilter {

private:
    std::string m_filterName = "";

public:
    /**
	 * @brief Sample intensity weigh
	 * Sample the photon weight from this filter. The mean of photon weights sampled is 1.0 for N samples when N->infinity.  
	 * @param angle Angle in radians.
	 * @return double Photon weight
	 */
    virtual T sampleIntensityWeight(const T angle) const = 0;
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
};

/**
 * @brief Class for simple CT bowtie filter modeling. 
 * 
 * Filter to adjust photon weight based on a measured bowtie filter profile. Note that beamhardening of a specter is not modeled
 * and the bowtie filter simply adjust photon fluence according to a profile. 
 */
template <Floating T = double>
class BowTieFilter : public BeamFilter<T> {
private:
    std::vector<std::pair<T, T>> m_data;

protected:
    /**
     * @brief Normalize filter
	 * Normalize filter such that expectation value of sampleIntensityWeight approches unity for large number of samples.
    */
    void normalizeData()
    {
        const auto mean = std::reduce(m_data.begin(), m_data.end(), 0.0, [](auto a, auto el) { return a + el.second; }) / m_data.size();
        std::transform(m_data.begin(), m_data.end(), m_data.begin(), [=](const auto& el) { return std::make_pair(el.first, el.second / mean); });
    }

public:
    /**
	 * @brief Construct a new Bow Tie Filter object
	 * Construct a new Bow Tie Filter object
	 * 
	 * @param angles A vector of angles in radians for the fluence profile. The filter is assumed to be symmetrical and only absolute values of angles are used.
	 * @param weights Photon fluence for the corresponding angles. The fluence measurements are normalized such that the average photon weight is 1.0 when number of samples -> infinity. 
	 */
    BowTieFilter(const std::vector<T>& angles, const std::vector<T>& weights)
        : BeamFilter<T>()
    {
        if (angles.size() == weights.size()) {
            m_data.resize(angles.size());
            for (std::size_t i = 0; i < angles.size(); i++) {
                m_data[i].first = std::abs(angles[i]);
                m_data[i].second = weights[i];
            }
        }
        std::sort(m_data.begin(), m_data.end());
        normalizeData();
    }

    /**
	 * @brief Construct a new Bow Tie Filter object
	 * Construct a new Bow Tie Filter object
	 * @param angleWeightsPairs A vector of angle in radians and photon fluence pairs.
	 */
    BowTieFilter(const std::vector<std::pair<T, T>>& angleWeightsPairs)
        : BeamFilter<T>()
        , m_data(angleWeightsPairs)
    {
        for (auto& p : m_data) {
            p.first = std::abs(p.first);
        }
        std::sort(m_data.begin(), m_data.end());
        normalizeData();
    }

    /**
	 * @brief Construct a new Bow Tie Filter object from another BowTieFilter
	 * Construct a new Bow Tie Filter object
	 * @param other BowTiefilter
	 */
    BowTieFilter(const BowTieFilter& other)
        : BeamFilter<T>()
    {
        m_data = other.data();
    }
    virtual ~BowTieFilter() = default;

    /**
	 * @brief Sample photon weight
	 * Sample photon weight modifier for a photon with direction angle in radians along the first Source cosine direction 
	 * @param angle Photon angle in radians 
	 * @return double Photon fluence weight modifier 
	 */
    T sampleIntensityWeight(T anglePlusAndMinus) const override
    {

        const auto angle = std::abs(anglePlusAndMinus);

        auto first = m_data.begin();
        auto last = m_data.end();
        std::advance(last, -1);

        //binary search for angle
        auto mid = std::distance(first, last) / 2;
        auto it = first;
        std::advance(it, mid);
        while (it != first) {
            if (angle < it->first) {
                last = it;
                it = first;
            } else {
                first = it;
            }
            mid = std::distance(first, last) / 2;
            std::advance(it, mid);
        }
        // end search

        if (angle < first->first) {
            return m_data[0].second;
        }
        if (angle > last->first) {
            return m_data.back().second;
        }

        //linear interpolation

        const auto x0 = first->first;
        const auto x1 = last->first;
        const auto y0 = first->second;
        const auto y1 = last->second;

        return y0 + (angle - x0) * (y1 - y0) / (x1 - x0);
    }

    /**
	 * @brief Filter data
	 * Returns a vector reference to angle in radians and normalized fluence weight pairs  
	 * @return const std::vector<std::pair<double, double>>& Vector of <angle, weight> pairs
	 */
    const std::vector<std::pair<T, T>>& data() const { return m_data; }
};

/**
 * @brief Filter for modeling of organ exposure control for Siemens CT scanners (XCare).
 * 
 * This filter modifies a particle's weight along the rotation angle in the same manner as Simens have iplementet on some CT scanner models.
 * The result is decreased fluence (photon weight) along the filters span angle and increased fluence outside the filter angle. The mean photon weight over all angles in a rotation is unity. 
 */
template <Floating T = double>
class XCareFilter : public BeamFilter<T> {

private:
    T m_filterAngle;
    T m_spanAngle;
    T m_rampAngle;
    T m_lowWeight;

public:
    /**
	 * @brief Construct a new XCareFilter object
	 * Constructs new XCareFilter with default values
	 */
    XCareFilter()
        : BeamFilter<T>()
    {
        m_filterAngle = 0;
        m_spanAngle = 120 * DEG_TO_RAD<T>();
        m_rampAngle = 20 * DEG_TO_RAD<T>();
        m_lowWeight = T { 0.6 };
    }
    virtual ~XCareFilter() = default;

    /**
	 * @brief Center filter angle
	 * 
	 * Center angle for the filter where photon weight is reduced
	 * @return double Angle in radians
	 */
    T filterAngle() const { return m_filterAngle; }

    /**
	 * @brief Center filter angle
	 * 
	 * Center angle for the filter where photon weight is reduced
	 * @return double Angle in degrees
	 */
    T filterAngleDeg() const
    {
        return m_filterAngle * RAD_TO_DEG<T>();
    }

    /**
	 * @brief Set the Filter Angle
	 * 
	 * Center angle for the filter where photon weight is reduced
	 * @param angle New filter angle in radians 
	 */
    void setFilterAngle(T angle)
    {
        constexpr T pi_2 = T { 2 } * PI_VAL<T>();
        m_filterAngle = std::fmod(angle, pi_2);
        if (m_filterAngle < 0.0)
            m_filterAngle += pi_2;
    }

    /**
	 * @brief Set the Filter Angle
	 * 
	 * Center angle for the filter where photon weight is reduced
	 * @param angle New filter angle in degrees
	 */
    void setFilterAngleDeg(T angle)
    {
        setFilterAngle(angle * DEG_TO_RAD<T>());
    }

    /**
	 * @brief Span angle of the filter
	 * 
	 * The span where photon weight is reduced, this angle span includes the ramp angle. Valid values are 5 degrees > span angle > PI. 
	 * @return double Span angle in radians
	 */
    T spanAngle() const { return m_spanAngle; }

    /**
	 * @brief Span angle of the filter
	 * 
	 * The span where photon weight is reduced, this angle span includes the ramp angle. Valid values are 5 degrees > span angle > PI. 
	 * @return double Span angle in degrees
	 */
    T spanAngleDeg() const
    {
        return m_spanAngle * RAD_TO_DEG<T>();
    }

    /**
	 * @brief Set span angle of the filter
	 * 
	 * The span where photon weight is reduced, this angle span includes the ramp angle. Valid values are 5 degrees > span angle > PI. 
	 * @return T Span angle in radians
	 */
    void setSpanAngle(T angle)
    {
        constexpr T smallestDegree = T { 5.0 } * DEG_TO_RAD<T>();
        if ((angle > smallestDegree) && (angle < PI_VAL<T>())) {
            m_spanAngle = angle;
        }
    }
    /**
	 * @brief Set span angle of the filter
	 * 
	 * The span where photon weight is reduced, this angle span includes the ramp angle. Valid values are 5 degrees > span angle > PI. 
	 * @return double Span angle in degrees
	 */
    void setSpanAngleDeg(T angle)
    {
        setSpanAngle(angle * DEG_TO_RAD<T>());
    }

    /**
	 * @brief Ramp angle 
	 * 
	 * The ramp angle span for ramping down the photon weight. The ramp angle are included in the span angle and must be less than spanAngle/2. 
	 * @return double Ramp angle in radians
	 */
    T rampAngle() const { return m_rampAngle; }

    /**
	 * @brief Ramp angle 
	 * 
	 * The ramp angle span for ramping down the photon weight. The ramp angle are included in the span angle and must be less than spanAngle/2. 
	 * @return double Ramp angle in degrees
	 */
    T rampAngleDeg() const
    {
        return m_rampAngle * RAD_TO_DEG<T>();
    }

    /**
	 * @brief Set ramp angle
	 * 
	 * The ramp angle span for ramping down the photon weight. The ramp angle are included in the span angle and must be less than spanAngle/2. 
	 * @param angle ramp angle in radians 
	 */
    void setRampAngle(T angle)
    {
        if ((angle >= 0.0) && (angle <= 0.5 * m_spanAngle)) {
            m_rampAngle = angle;
        }
    }

    /**
	 * @brief Set ramp angle
	 * 
	 * The ramp angle span for ramping down the photon weight. The ramp angle are included in the span angle and must be less than spanAngle/2. 
	 * @param angle ramp angle in degrees 
	 */
    void setRampAngleDeg(T angle)
    {
        setRampAngle(angle * DEG_TO_RAD<T>());
    }

    /**
	 * @brief Photon weight in the filters span angle excluding angle ramps 
	 * 
	 * Photon weight in the span angle excluding ramp angles. Must be in 0.0 < lowWeight <= 1.0.
	 * @return double 
	 */
    T lowWeight() const { return m_lowWeight; }

    /**
	 * @brief Set photon weight in the filters span angle excluding angle ramps 
	 * 
	 * Photon weight in the span angle excluding ramp angles. Must be in 0.0 < lowWeight <= 1.0.
	 * @param weight
	 */
    void setLowWeight(T weight)
    {
        if ((weight > 0.0) && (weight <= 1.0)) {
            m_lowWeight = weight;
        }
    }
    /**
	 * @brief Photon weight for angles outside span angle
	 * 
	 * Highest photon weight outside span angle. This value is calculated such that expectation weight across all angles is 1.0
	 * @return double Photon hight weight
	 */
    T highWeight() const
    {
        constexpr T pi_2 = T { 2 } * PI_VAL<T>();
        return (pi_2 - m_spanAngle * m_lowWeight + m_lowWeight * m_rampAngle) / (pi_2 - m_spanAngle + m_rampAngle);
    }

    /**
	 * @brief Sample photon weight
	 * 
	 * Sample photon weight from rotation angle. Mean weight from all angles is 1.0
	 * 
	 * @param angle Rotation angle in radians
	 * @return double Photon weight modifier
	 */
    T sampleIntensityWeight(const T angle) const override
    {
        constexpr T pi_2 = T { 2 } * PI_VAL<T>();
        auto angleMod = std::fmod(angle - m_filterAngle + PI_VAL<T>(), pi_2); // centering angle on 180 degrees
        if (angleMod < 0)
            angleMod += pi_2;

        const auto high = highWeight();

        const auto startFilter = PI_VAL<T>() - m_spanAngle * T { 0.5 };
        if (angleMod < startFilter)
            return high;

        const auto endRamp1 = startFilter + m_rampAngle;
        if (angleMod < endRamp1)
            return interp<T>(startFilter, endRamp1, high, m_lowWeight, angleMod);

        const auto startRamp2 = endRamp1 + m_spanAngle - m_rampAngle;
        if (angleMod < startRamp2)
            return m_lowWeight;

        const auto endramp2 = startFilter + m_spanAngle;
        if (angleMod < endramp2)
            return interp<T>(startRamp2, endramp2, m_lowWeight, high, angleMod);

        return high;
    }
};

/**
 * @brief Filter to model Heel effect of a tube.
 * 
 * Filter for modeling of Heel effect of a x-ray tube. The filter do not model beam hardening of the Heel effect, only photon fluence effects. 
 * 
 */
template <Floating T = double>
class HeelFilter {
private:
    T m_energyStep = 2.0;
    T m_energyStart = 10.0;
    std::size_t m_energySize = 65;

    T m_angleStep = 0.07; // about 4 deg step size
    T m_angleStart = 0.07;
    std::size_t m_angleSize = 5;
    std::vector<T> m_energies;
    std::vector<T> m_weights; //vector of m_angleSize*m_energySize weights
    std::string m_filterName = "";

public:
    /**
     * @brief Constructs a new Heel filter
     * @param tube Tube to use when calculating Heel effect
     * @param heel_angle_span Span angle in radians to calculate Heel effect. This should be atleast as large as the beam collimation from the radiation source.
    */
    HeelFilter(const Tube<T>& tube, const T heel_angle_span = 0.0)
    {
        update(tube, heel_angle_span);
    }

    /**
     * @brief Update this filter with a new tube.
     * @param tube Tube to use when calculating Heel effect
     * @param heel_angle_span Span angle in radians to calculate Heel effect. This should be atleast as large as the beam collimation from the radiation source.
    */
    void update(const Tube<T>& tube, const T heel_angle_span = 0.0)
    {
        // recalculating energy range
        m_energySize = static_cast<std::size_t>((tube.voltage() - m_energyStart) / m_energyStep);
        if (m_energySize < 2)
            m_energySize = 2;
        m_energies.clear();
        m_energies.reserve(m_energySize);
        for (std::size_t i = 0; i < m_energySize; ++i)
            m_energies.push_back(m_energyStart + i * m_energyStep);

        // recalculating weights
        const auto abs_span_angle = std::abs(heel_angle_span);
        if (tube.anodeAngle() < abs_span_angle / 2) // minimum span angle is greater than anode angle
            m_angleStart = -tube.anodeAngle();
        else
            m_angleStart = -abs_span_angle / 2;
        m_angleStep = (abs_span_angle / 2 - m_angleStart) / m_angleSize;

        m_weights.clear();
        m_weights.resize(m_energySize * m_angleSize);
        for (std::size_t i = 0; i < m_angleSize; ++i) {
            const T angle = m_angleStart + i * m_angleStep + tube.anodeAngle();
            auto specter = tube.getSpecter(m_energies, angle, false);
            for (std::size_t j = 0; j < m_energySize; ++j)
                m_weights[j * m_angleSize + i] = specter[j];
        }
        //normalizing weights
        auto w_start = m_weights.begin();
        for (std::size_t i = 0; i < m_energySize; ++i) {
            auto w_end = w_start + m_angleSize;
            //normalize to mean
            const T sum = std::reduce(w_start, w_end, 0.0) / m_angleSize;
            if (sum > T { 0 }) {
                w_start = std::transform(w_start, w_end, w_start, [=](auto w) -> T { return w / sum; });
            } else {
                std::fill(w_start, w_end, T { 1 });
                w_start = w_end;
            }
        }
    }

    /**
     * @brief Sample weight from angle and photon energy
	 * The expectation value of weight from a large number of samples for a specific pfoton energy is normalized to unity. 
     * @param angle Angle of photon direction in anode cathode direction of a tube.
     * @param energy Ebergy in keV for photon. 
     * @return Weight of photon.
    */
    T sampleIntensityWeight(const T angle, const T energy) const
    {
        std::size_t e_index = static_cast<std::size_t>((energy - m_energyStart + T { 0.5 } * m_energyStep) / m_energyStep);
        if (e_index >= m_energySize)
            e_index = m_energySize - 1;
        if (energy < m_energyStart)
            e_index = 0;

        std::size_t a_index = static_cast<std::size_t>((angle - m_angleStart) / m_angleStep); // we want lowest angle index since we are interpolating
        if (a_index >= m_angleSize)
            a_index = m_angleSize - 1;
        if (angle < m_angleStart)
            a_index = 0;

        std::size_t w_index = e_index * m_angleSize + a_index;
        if (a_index < m_angleSize - 1) // we want to interpolate if not on edge
        {
            const auto a0 = m_angleStart + m_angleStep * a_index;
            const auto a1 = m_angleStart + m_angleStep * (a_index + 1);
            const auto w0 = m_weights[w_index];
            const auto w1 = m_weights[w_index + 1];
            return interp(a0, a1, w0, w1, angle);
        }
        return m_weights[w_index];
    }
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
    const std::vector<T>& weights() const { return m_weights; }
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
};

/**
 * @brief Filter to adjust photon weights according to a tube current profile for CT examinations.
 * Filter to simulate automatic exposure control for CT examinations. This filter will match a given tube current profile to a denisty image to generate a lookup table of photon weights according to a ddensity profile. 
*/
template <Floating T = double>
class AECFilter {

private:
    std::vector<T> m_mass;
    std::vector<T> m_massIntensity;
    T m_positionStep = 0.0;
    T m_positionMin = 0.0;
    T m_positionMax = 0.0;
    std::vector<T> m_positionIntensity;
    std::string m_filterName = "";
    bool m_valid = false;

protected:
    void generateMassWeightMap(typename std::vector<T>::const_iterator densBeg, typename std::vector<T>::const_iterator densEnd, const std::array<T, 3>& spacing, const std::array<std::size_t, 3>& dimensions, const std::vector<T>& exposuremapping)
    {
        m_valid = false;
        if (std::distance(densBeg, densEnd) != dimensions[0] * dimensions[1] * dimensions[2]) {
            m_mass.resize(1);
            m_massIntensity.resize(1);
            m_massIntensity[0] = 1;
            return;
        }
        if (exposuremapping.size() != dimensions[2]) {
            m_mass.resize(1);
            m_massIntensity.resize(1);
            m_massIntensity[0] = 1;
            return;
        }

        std::vector<std::pair<T, T>> posExp(dimensions[2]);
        const auto sliceStep = dimensions[0] * dimensions[1];
        const auto voxelArea = spacing[0] * spacing[1];
        for (std::size_t i = 0; i < dimensions[2]; ++i) {
            const auto start = densBeg + sliceStep * i;
            const auto sliceMass = std::reduce(std::execution::par_unseq, start, start + sliceStep, 0.0) * voxelArea;
            posExp[i] = std::make_pair(sliceMass, exposuremapping[i]);
        }
        std::sort(posExp.begin(), posExp.end());

        m_mass.resize(dimensions[2]); // this should now be sorted
        std::transform(std::execution::par_unseq, posExp.cbegin(), posExp.cend(), m_mass.begin(), [](auto pair) -> T { return pair.first; });

        const auto meanExposure = std::reduce(std::execution::par_unseq, exposuremapping.cbegin(), exposuremapping.cend(), T { 0.0 }) / dimensions[2];
        m_massIntensity.resize(dimensions[2]);
        std::transform(std::execution::par_unseq, posExp.cbegin(), posExp.cend(), m_massIntensity.begin(), [=](auto pair) -> T { return pair.second / meanExposure; });

        std::array<T, 3> origin({ 0, 0, 0 });
        generatePositionWeightMap(densBeg, densEnd, spacing, dimensions, origin);
        return;
    }

    void generatePositionWeightMap(typename std::vector<T>::const_iterator densBeg, typename std::vector<T>::const_iterator densEnd, const std::array<T, 3> spacing, const std::array<std::size_t, 3>& dimensions, const std::array<T, 3>& origin)
    {
        m_positionMin = origin[2] - spacing[2] * dimensions[2] * T { 0.5 };
        m_positionMax = m_positionMin + spacing[2] * dimensions[2];
        m_positionStep = (m_positionMax - m_positionMin) / dimensions[2]; // this should be equal to spacing[2]
        m_positionIntensity.resize(dimensions[2]);

        const auto sliceStep = dimensions[0] * dimensions[1];
        const auto voxelArea = spacing[0] * spacing[1];
        for (std::size_t i = 0; i < dimensions[2]; ++i) {
            auto start = densBeg + sliceStep * i;
            const auto sliceMass = std::reduce(std::execution::par_unseq, start, start + sliceStep, T { 0.0 }) * voxelArea;
            m_positionIntensity[i] = interpolateMassIntensity(sliceMass);
        }
        m_valid = true;
    }
    T interpolateMassIntensity(T mass) const
    {
        auto beg = m_mass.cbegin();
        auto end = m_mass.cend();
        auto pos = std::upper_bound(m_mass.cbegin(), m_mass.cend(), mass);
        if (pos == beg)
            return m_massIntensity[0];
        if (mass > m_mass.back())
            return m_massIntensity.back();

        if (pos == end)
            return m_massIntensity.back();
        const auto x0 = *(pos - 1);
        const auto x1 = *pos;
        const auto y0 = *(m_massIntensity.cbegin() + std::distance(beg, pos) - 1);
        const auto y1 = *(m_massIntensity.cbegin() + std::distance(beg, pos));
        return interp(x0, x1, y0, y1, mass);
    }

public:
    /**
     * @brief Constructs a new AECFilter
     * @param densityImage Density array of size dimension[0]*dimension[1]*dimension[2]. 
     * @param spacing Spacing of each voxel for the density array.
     * @param dimensions Dimensions of the densityImage.
     * @param exposuremapping Exposure along the z (3rd) dimension, must be of size dimension[2].
    */
    AECFilter(const std::vector<T>& densityImage, const std::array<T, 3> spacing, const std::array<std::size_t, 3> dimensions, const std::vector<T>& exposuremapping)
    {
        generateMassWeightMap(densityImage.cbegin(), densityImage.cend(), spacing, dimensions, exposuremapping);
    }

    /**
     * @brief Constructs a new AECFilter
     * @param densityImage Density array of size dimension[0]*dimension[1]*dimension[2]. 
     * @param spacing Spacing of each voxel for the density array.
     * @param dimensions Dimensions of the densityImage.
     * @param exposuremapping Exposure along the z (3rd) dimension, must be of size dimension[2].
    */
    AECFilter(std::shared_ptr<std::vector<T>>& densityImage, const std::array<T, 3> spacing, const std::array<std::size_t, 3> dimensions, const std::vector<T>& exposuremapping)
    {
        generateMassWeightMap(densityImage->cbegin(), densityImage->cend(), spacing, dimensions, exposuremapping);
    }
    /**
     * @brief Constructs a new AECFilter directly from mass and weight arrays
     * @param mass Array of masses. Must be of same size as intensity array.
     * @param intensity Array of intensity or exposure corresponding to each mass in mass array. Must be of same size as mass array
    */
    AECFilter(const std::vector<T>& mass, const std::vector<T>& intensity)
    {
        m_mass = mass;
        m_massIntensity = intensity;
        m_positionIntensity.resize(1);
        m_positionIntensity[0] = 1.0;
        m_positionStep = 1.0;
        m_positionMax = 1.0;
        m_valid = false;
    }

    /**
     * @brief Sample weight of photon from a position in a world object.
     * Sample AEC weight of photon originating from position. This function assumes that updateFromWorld is called with the same World object the position is generated from. The expected value of weight from a large number of samples across all valid positions (the extent of the world object) is unity.
     * @param position Position to interpolate photon weight.
     * @return photon weight
    */
    T sampleIntensityWeight(const std::array<T, 3>& position) const
    {
        const T p = position[2];
        if (p < m_positionMin + m_positionStep)
            return m_positionIntensity[0];
        if (p >= m_positionMax - m_positionStep)
            return m_positionIntensity.back();
        const std::size_t ind = static_cast<std::size_t>((p - m_positionMin) / m_positionStep);
        const auto x0 = ind * m_positionStep + m_positionMin;
        const auto x1 = x0 + m_positionStep;
        const auto y0 = m_positionIntensity[ind];
        const auto y1 = m_positionIntensity[ind + 1];
        return interp(x0, x1, y0, y1, p);
    }

    /**
     * @brief Update AECFilter from a World. 
     * This function will generate position lookup tables based on world, such as that photon weights can be sampled from a position within the world object.
     * @param world 
    */

    void updateFromWorld(const World<T>& world)
    {
        const auto& dens = world.densityArray();
        const auto& spacing = world.spacing();
        const auto& dim = world.dimensions();
        const auto& origin = world.origin();
        generatePositionWeightMap(dens->cbegin(), dens->cend(), spacing, dim, origin);
    }
    /**
     * @brief Returns if the AECFilter is valid and has been updated from a world.
     * @return 
    */
    bool isValid() const { return m_valid; }
    /**
     * @brief Masses used for interpolation of photon weights
     * @return Array of mass values
    */
    const std::vector<T>& mass() const { return m_mass; }
    /**
     * @brief Intensities used for interpolation of photon weights
     * @return Array of intensity values.
    */
    const std::vector<T>& massIntensity() const { return m_massIntensity; }
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
};
}