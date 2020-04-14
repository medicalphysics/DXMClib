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

#include "dxmc/beamfilters.h"
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <execution>

constexpr double PI = 3.14159265359;
constexpr double PI_2 = PI + PI;
constexpr double DEG_TO_RAD = PI / 180.0;
constexpr double RAD_TO_DEG = 180.0 / PI;

BowTieFilter::BowTieFilter(const std::vector<double>& angles, const std::vector<double>& weights)
	:BeamFilter()
{
    if (angles.size() == weights.size())
    {
        m_data.resize(angles.size());
        for (std::size_t i = 0; i < angles.size(); i++)
        {
            m_data[i].first = std::abs(angles[i]);
            m_data[i].second = weights[i];
        }
    }
    std::sort(m_data.begin(), m_data.end());
	normalizeData();
}

BowTieFilter::BowTieFilter(const std::vector<std::pair<double, double> > & angleWeightsPairs)
    :BeamFilter(), m_data(angleWeightsPairs)
{
    for (auto& p : m_data)
    {
        p.first = std::abs(p.first);
    }
    std::sort(m_data.begin(), m_data.end());
	normalizeData();
}

BowTieFilter::BowTieFilter(const BowTieFilter& other)
	:BeamFilter()
{
	m_data = other.data();
}

double BowTieFilter::sampleIntensityWeight(const double anglePlusAndMinus) const
{

	const double angle = std::abs(anglePlusAndMinus);

    auto first = m_data.begin();
    auto last = m_data.end();
    std::advance(last, -1);

    //binary search for angle
    auto mid = std::distance(first, last) / 2;
    auto it = first;
    std::advance(it, mid);
    while (it != first)
    {
        if (angle < it->first)
        {
            last = it;
            it = first;
        }
        else
        {
            first = it;
        }
        mid = std::distance(first, last) / 2;
        std::advance(it, mid);
    }
    // end search


    if (angle < first->first)
    {
        return m_data[0].second;
    }
    if (angle > last->first)
    {
        return m_data.back().second;
    }

    //linear interpolation

    const double x0 = first->first;
    const double x1 = last->first;
    const double y0 = first->second;
    const double y1 = last->second;

    return y0 + (angle - x0) * (y1 - y0) / (x1 - x0);
}

void BowTieFilter::normalizeData()
{
	const auto mean = std::reduce(m_data.begin(), m_data.end(), 0.0, [](auto a, auto el) {return a + el.second; }) / m_data.size();
	std::transform(m_data.begin(), m_data.end(), m_data.begin(), [=](const auto &el) {return std::make_pair(el.first, el.second / mean); });	
}


XCareFilter::XCareFilter()
	:BeamFilter()
{
	m_filterAngle = 0.0;
	m_spanAngle = 120.0 * DEG_TO_RAD;
	m_rampAngle = 20.0 * DEG_TO_RAD;
	m_lowWeight = 0.6;
}

double XCareFilter::filterAngleDeg() const
{
	return m_filterAngle * RAD_TO_DEG;
}

void XCareFilter::setFilterAngle(double angle)
{
	m_filterAngle = std::fmod(angle, PI_2);
	if (m_filterAngle < 0.0)
		m_filterAngle += PI_2;
}

void XCareFilter::setFilterAngleDeg(double angle)
{
	setFilterAngle(angle*DEG_TO_RAD);
}

double XCareFilter::spanAngleDeg() const
{
	return m_spanAngle * RAD_TO_DEG;
}
void XCareFilter::setSpanAngle(double angle)
{
	constexpr double smallestDegree = 5.0 * DEG_TO_RAD;
	if ((angle > smallestDegree) && (angle < PI))
	{
		m_spanAngle = angle;
	}
}

void XCareFilter::setSpanAngleDeg(double angle)
{
	setSpanAngle(angle * DEG_TO_RAD);
}

double XCareFilter::rampAngleDeg() const
{
	return m_rampAngle * RAD_TO_DEG;
}

void XCareFilter::setRampAngle(double angle)
{
	if ((angle >= 0.0) && (angle <= 0.5 * m_spanAngle))
	{
		m_rampAngle = angle;
	}
}

void XCareFilter::setRampAngleDeg(double angle)
{
	setRampAngle(angle * DEG_TO_RAD);
}


void XCareFilter::setLowWeight(double weight)
{
	if ((weight > 0.0) && (weight <= 1.0))
	{
		m_lowWeight = weight;
	}
}

double XCareFilter::highWeight() const
{
	return (PI_2 - m_spanAngle * m_lowWeight + m_lowWeight * m_rampAngle) / (PI_2 - m_spanAngle + m_rampAngle);
}

template<typename T>
inline T interp(T x0, T x1, T y0, T y1, T x)
{
	return y0 + (y1 - y0) * (x - x0) / (x1 - x0);
}
template<typename T>
inline T interp(T x[2], T y[2], T xi)
{
	return y[0] + (y[1] - y[0]) * (xi - x[0]) / (x[1] - x[0]);
}


double XCareFilter::sampleIntensityWeight(const double angle) const
{
	double angleMod = std::fmod(angle - m_filterAngle + PI, PI_2); // centering angle on 180 degrees
	if (angleMod < 0.0)
		angleMod += PI_2;

	const double high = highWeight();
	
	const double startFilter = PI - m_spanAngle * 0.5;
	if (angleMod < startFilter)
		return high;

	const double endRamp1 = startFilter + m_rampAngle;
	if (angleMod < endRamp1)
		return interp(startFilter, endRamp1, high, m_lowWeight, angleMod);

	const double startRamp2 = endRamp1 + m_spanAngle - m_rampAngle;
	if (angleMod < startRamp2)
		return m_lowWeight;

	const double endramp2 = startFilter + m_spanAngle;
	if (angleMod < endramp2)
		return interp(startRamp2, endramp2, m_lowWeight, high, angleMod);
	
	return high;
}



AECFilter::AECFilter(const std::vector<double>& densityImage, const std::array<double, 3> spacing, const std::array<std::size_t, 3> dimensions, const std::vector<double>& exposure)
{
	generateMassWeightMap(densityImage.cbegin(), densityImage.cend(), spacing, dimensions, exposure);
}

AECFilter::AECFilter(std::shared_ptr<std::vector<double>>& densityImage, const std::array<double, 3> spacing, const std::array<std::size_t, 3> dimensions, const std::vector<double>& exposure)
{
	generateMassWeightMap(densityImage->cbegin(), densityImage->cend(), spacing, dimensions, exposure);
}

AECFilter::AECFilter(const std::vector<double>& mass, const std::vector<double>& intensity)
{
	m_mass = mass;
	m_massIntensity = intensity;
	m_positionIntensity.resize(1);
	m_positionIntensity[0] = 1.0;
	m_positionStep = 1.0;
	m_positionMax = 1.0;
	m_valid = false;
}

double AECFilter::sampleIntensityWeight(const std::array<double, 3>& position) const
{
	const double p = position[2];
	if (p < m_positionMin + m_positionStep)
		return m_positionIntensity[0];
	if (p >= m_positionMax-m_positionStep)
		return m_positionIntensity.back();
	const std::size_t ind = static_cast<std::size_t>((p - m_positionMin) / m_positionStep);
	const double x0 = ind * m_positionStep + m_positionMin;
	const double x1 = x0 + m_positionStep;
	const double y0 = m_positionIntensity[ind];
	const double y1 = m_positionIntensity[ind+1];
	return interp(x0, x1, y0, y1, p);
}

void AECFilter::updateFromWorld(const World& world)
{
	auto dens = world.densityArray();
	auto spacing = world.spacing();
	auto dim = world.dimensions();
	auto origin = world.origin();
	generatePositionWeightMap(dens->cbegin(), dens->cend(), spacing, dim, origin);
}

void AECFilter::generateMassWeightMap(std::vector<double>::const_iterator densBeg, std::vector<double>::const_iterator densEnd, const std::array<double, 3> spacing, const std::array<std::size_t, 3> dimensions, const std::vector<double>& exposuremapping)
{
	m_valid = false;
	if (std::distance(densBeg, densEnd) != dimensions[0] * dimensions[1] * dimensions[2])
	{
		m_mass.resize(1);
		m_massIntensity.resize(1);
		m_massIntensity[0] = 1;
		return;
	}
	if (exposuremapping.size() != dimensions[2])
	{
		m_mass.resize(1);
		m_massIntensity.resize(1);
		m_massIntensity[0] = 1;
		return;
	}
	
	std::vector<std::pair<double, double>> posExp(dimensions[2]);
	const auto sliceStep = dimensions[0] * dimensions[1];
	const double voxelArea = spacing[0] * spacing[1];
	for(std::size_t i=0; i < dimensions[2]; ++i)
	{
		const auto start = densBeg + sliceStep * i;
		const double sliceMass = std::reduce(std::execution::par_unseq, start, start + sliceStep, 0.0) * voxelArea;
		posExp[i] = std::make_pair(sliceMass, exposuremapping[i]);
	}
	std::sort(posExp.begin(), posExp.end());

	m_mass.resize(dimensions[2]); // this should now be sorted
	std::transform(std::execution::par_unseq, posExp.cbegin(), posExp.cend(), m_mass.begin(), [](auto pair)->double {return pair.first; });


	const double meanExposure = std::reduce(std::execution::par_unseq, exposuremapping.cbegin(), exposuremapping.cend(), 0.0) / dimensions[2];
	m_massIntensity.resize(dimensions[2]);
	std::transform(std::execution::par_unseq, posExp.cbegin(), posExp.cend(), m_massIntensity.begin(), [=](auto pair)->double {return pair.second / meanExposure; });
	
	std::array<double, 3> origin({ 0,0,0 });
	generatePositionWeightMap(densBeg, densEnd, spacing, dimensions, origin);
	return;
}

void AECFilter::generatePositionWeightMap(std::vector<double>::const_iterator densBeg, std::vector<double>::const_iterator densEnd, const std::array<double, 3> spacing, const std::array<std::size_t, 3> dimensions, const std::array<double, 3>& origin)
{
	m_positionMin = origin[2] - spacing[2] * dimensions[2] * 0.5;
	m_positionMax = m_positionMin + spacing[2] * dimensions[2];
	m_positionStep = (m_positionMax - m_positionMin) / dimensions[2]; // this should be equal to spacing[2]
	m_positionIntensity.resize(dimensions[2]);

	const auto sliceStep = dimensions[0] * dimensions[1];
	const double voxelArea = spacing[0] * spacing[1];
	for (std::size_t i = 0; i < dimensions[2]; ++i)
	{
		auto start = densBeg + sliceStep * i;
		const double sliceMass = std::reduce(std::execution::par_unseq, start, start + sliceStep, 0.0) * voxelArea;
		m_positionIntensity[i] = interpolateMassIntensity(sliceMass);
	}
	m_valid = true;
}

double AECFilter::interpolateMassIntensity(double mass) const
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
	const double x0 = *(pos - 1);
	const double x1 = *pos;
	const double y0 = *(m_massIntensity.cbegin() + std::distance(beg, pos) - 1);
	const double y1 = *(m_massIntensity.cbegin() + std::distance(beg, pos));
	return interp(x0, x1, y0, y1, mass);
}

HeelFilter::HeelFilter(const Tube& tube, const double heel_angle_span)
{
	update(tube, heel_angle_span);
}

void HeelFilter::update(const Tube& tube, const double heel_angle_span)
{
	// recalculating energy range
	m_energySize = static_cast<std::size_t>((tube.voltage() - m_energyStart) / m_energyStep);
	if (m_energySize < 1)
		m_energySize = 1;
	m_energies.clear();
	m_energies.reserve(m_energySize);
	for (std::size_t i = 0; i < m_energySize; ++i)
		m_energies.push_back(m_energyStart + i * m_energyStep);

	// recalculating weights
	const double abs_span_angle = std::abs(heel_angle_span);
	if (tube.anodeAngle() < abs_span_angle * 0.5) // minimum span angle is greater than anode angle
		m_angleStart = -tube.anodeAngle();
	else
		m_angleStart = -abs_span_angle * 0.5;
	m_angleStep = (abs_span_angle * 0.5 - m_angleStart) / m_angleSize;

	m_weights.clear();
	m_weights.resize(m_energySize * m_angleSize);
	for (std::size_t i = 0; i < m_angleSize; ++i)
	{
		const double angle = m_angleStart + i * m_angleStep + tube.anodeAngle();
		auto specter = tube.getSpecter(m_energies, angle, false);
		for (std::size_t j = 0; j < m_energySize; ++j)
			m_weights[j * m_angleSize + i] = specter[j];
	}
	//normalizing weights
	auto w_start = m_weights.begin();
	for (std::size_t i = 0; i < m_energySize; ++i)
	{
		auto w_end = w_start + m_angleSize;
		//normalize to mean
		const double sum = std::reduce(w_start, w_end, 0.0) / m_angleSize;
		w_start = std::transform(w_start, w_end, w_start,  [=](double w) ->double {return w / sum; });
	}
}

double HeelFilter::sampleIntensityWeight(const double angle, const double energy) const
{
	std::size_t e_index = (energy - m_energyStart + 0.5 * m_energyStep) / m_energyStep;
	if (e_index >= m_energySize)
		e_index = m_energySize - 1;
	if (energy < m_energyStart)
		e_index = 0;

	std::size_t a_index = (angle - m_angleStart ) / m_angleStep; // we want lowest angle index since we are interpolating
	if (a_index >= m_angleSize)
		a_index = m_angleSize - 1;
	if (angle < m_angleStart)
		a_index = 0;

	std::size_t w_index = e_index * m_angleSize + a_index;
	if (a_index < m_angleSize - 1) // we want to interpolate if not on edge
	{
		const double a0 = m_angleStart + m_angleStep * a_index;
		const double a1 = m_angleStart + m_angleStep * (a_index + 1);
		const double w0 = m_weights[w_index];
		const double w1 = m_weights[w_index + 1];
		return interp(a0, a1, w0, w1, angle);
	}
	return m_weights[w_index];
}
