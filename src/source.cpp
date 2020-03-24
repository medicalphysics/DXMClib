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

#include "dxmc/source.h"
#include "dxmc/vectormath.h"
#include "dxmc/transport.h"

#include <future>
#include <numeric>

constexpr double PI = 3.14159265358979323846;  /* pi */
constexpr double PI_2 = 2.0 * PI;
constexpr double DEG_TO_RAD = PI / 180.0;
constexpr double RAD_TO_DEG = 1.0 / DEG_TO_RAD;
constexpr double KEV_TO_MJ = 1.6021773e-13; // milli Joules
constexpr double ANGLE_ERRF = 1E-6;

constexpr std::uint64_t CTDI_MIN_HISTORIES = static_cast<std::uint64_t>(500E6);

template<typename T>
T vectorLenght(const std::array<T, 3>& v)
{
	return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}


Source::Source()
{
	m_type = Type::None;
	for (std::size_t i = 0; i < 3; i++)
	{
		m_position[i] = 0.0;
		m_directionCosines[i] = 0.0;
		m_directionCosines[i + 3] = 0.0;
	}
	m_directionCosines[0] = 1.0;
	m_directionCosines[4] = 1.0;
    m_historiesPerExposure = 1000000;
}

void Source::setDirectionCosines(const std::array<double, 6>& cosines)
{
	m_directionCosines = cosines;
	normalizeDirectionCosines();
}

void Source::normalizeDirectionCosines()
{
	vectormath::normalize(&m_directionCosines[0]);
	vectormath::normalize(&m_directionCosines[3]);
}

void Source::setHistoriesPerExposure(std::uint64_t histories)
{
	m_historiesPerExposure = histories;
}

std::uint64_t Source::historiesPerExposure(void) const
{
	return m_historiesPerExposure;
}


IsotropicSource::IsotropicSource()
	:Source()
{
	m_type = Type::Isotropic;
	
	// initializing a specterdistribution
	const std::vector<double> energies = { 60.0 };
	const std::vector<double> weights = { 1.0 };
	m_maxPhotonEnergy = 60.0;
	m_specterDistribution = std::make_unique<SpecterDistribution>(weights, energies);
}

bool IsotropicSource::getExposure(Exposure& exposure, std::uint64_t exposureNumber) const
{
	exposure.setNumberOfHistories(m_historiesPerExposure);
	exposure.setPosition(m_position);
	exposure.setDirectionCosines(m_directionCosines);
	exposure.setCollimationAngles(m_collimationAngles[0], m_collimationAngles[1]);

	exposure.setSpecterDistribution(m_specterDistribution.get());
	return exposureNumber < m_totalExposures;

}

bool IsotropicSource::isValid() const
{
	return m_specterDistribution != nullptr;
}

bool IsotropicSource::validate()
{
	return m_specterDistribution != nullptr;
}

void IsotropicSource::setSpecter(const std::vector<double> weights, const std::vector<double> energies)
{
	const auto maxEnergyElement = std::max_element(energies.cbegin(), energies.cend());
	m_maxPhotonEnergy = *maxEnergyElement;
	m_specterDistribution = std::make_unique<SpecterDistribution>(weights, energies);


}

void IsotropicSource::setCollimationAngles(double xRad, double yRad)
{
	if (xRad < -PI)
		m_collimationAngles[0] = -PI;
	else if (xRad > PI)
		m_collimationAngles[0] = PI;
	else
		m_collimationAngles[0] = xRad;

	constexpr double PI_H = PI * 0.5;

	if (yRad < -PI_H)
		m_collimationAngles[1] = -PI_H;
	else if (yRad > PI)
		m_collimationAngles[1] = PI_H;
	else
		m_collimationAngles[1] = yRad;
}

std::pair<double, double> IsotropicSource::collimationAngles() const
{
	return std::pair<double, double>(m_collimationAngles[0], m_collimationAngles[1]);
}



PencilSource::PencilSource()
	:Source()
{
	m_type = Type::Pencil;
}


bool PencilSource::getExposure(Exposure& exposure, std::uint64_t i) const
{
	exposure.setNumberOfHistories(m_historiesPerExposure);
	exposure.setPosition(m_position);
	exposure.setDirectionCosines(m_directionCosines);
	exposure.setCollimationAngles(0.0, 0.0);
	exposure.setMonoenergeticPhotonEnergy(m_photonEnergy);
	return i < m_totalExposures;
}
void PencilSource::setPhotonEnergy(double energy)
{
	if (energy > 500.0)
		m_photonEnergy = 500;
	else if (energy < 0.0)
		m_photonEnergy = 0.0;
	else
		m_photonEnergy = energy;
}

double PencilSource::getCalibrationValue(ProgressBar* progressBar) const
{
	Material airMaterial("Air, Dry (near sea level)");
	const double nHistories = totalExposures() * historiesPerExposure();
	const double calcOutput = nHistories * m_photonEnergy * airMaterial.getMassEnergyAbsorbtion(m_photonEnergy) * KEV_TO_MJ;
	const double factor = m_airDose / calcOutput;
	return factor;
}
DXSource::DXSource()
	:Source()
{
	m_type = Type::DX;
	m_sdd = 1000.0;
	m_fieldSize[0] = 100.0;
	m_fieldSize[1] = 100.0;
	updateFieldSize(m_fieldSize);
	m_tube.setAlFiltration(2.0);
	setDirectionCosines(zeroDirectionCosines());
}
bool DXSource::getExposure(Exposure& exposure, std::uint64_t i) const
{
	exposure.setNumberOfHistories(m_historiesPerExposure);
	exposure.setPosition(tubePosition());
	exposure.setDirectionCosines(m_directionCosines);
	exposure.setCollimationAngles(m_collimationAngles.data());
	exposure.setSpecterDistribution(m_specterDistribution.get());
	exposure.setHeelFilter(m_heelFilter.get());
	return i < m_totalExposures;
}

std::uint64_t DXSource::totalExposures() const
{
	return m_totalExposures;
}
void DXSource::setTotalExposures(std::uint64_t nExposures)
{
	if (nExposures > 1)
		m_totalExposures = nExposures;
	else
		m_totalExposures = 1;
}

const std::array<double, 2>& DXSource::collimationAngles() const
{
	return m_collimationAngles;
}
void DXSource::setCollimationAngles(const std::array<double, 2>& angles)
{
	std::array<double, 2> ang = { std::abs(angles[0]), std::abs(angles[1]) };
	updateCollimationAngles(ang);
}
const std::array<double, 2> DXSource::collimationAnglesDeg() const
{
	std::array<double, 2> a{ m_collimationAngles[0] * RAD_TO_DEG,
		m_collimationAngles[1] * RAD_TO_DEG };
	return a;
}
void DXSource::setCollimationAnglesDeg(const std::array<double, 2>& angles)
{
	std::array<double, 2> ang = { 
		std::abs(angles[0]) * DEG_TO_RAD, 
		std::abs(angles[1]) * DEG_TO_RAD 
	};
	updateCollimationAngles(ang);
}

void DXSource::setFieldSize(const std::array<double, 2>& mm)
{
	std::array<double, 2> amm = { std::abs(mm[0]), std::abs(mm[1]) };
	updateFieldSize(amm);
}

const std::array<double, 2>& DXSource::fieldSize(void) const
{
	return m_fieldSize;
}

void DXSource::updateSpecterDistribution()
{
	if (!m_specterValid)
	{
		auto energies = m_tube.getEnergy();
		auto n_obs = m_tube.getSpecter(energies);
		m_specterDistribution = std::make_unique<SpecterDistribution>(n_obs, energies);
		if (m_modelHeelEffect)
			m_heelFilter = std::make_unique<HeelFilter>(m_tube, m_collimationAngles[1]);
		else
			m_heelFilter = nullptr;
		m_specterValid = true;
	}
}

std::array<double, 6> DXSource::zeroDirectionCosines() const
{
	// changing this will break setSourceAngles
	std::array<double, 6> cos = { -1.0, .0, .0, .0, .0, 1.0 };
	return cos;
}


void DXSource::setSourceDetectorDistance(double mm)
{
	m_sdd = std::abs(mm);
	updateFieldSize(m_fieldSize);
}

double DXSource::sourceDetectorDistance() const
{
	return m_sdd;
}

void DXSource::setSourceAngles(double primaryAngle, double secondaryAngle)
{
	if (secondaryAngle > PI * 0.5 - ANGLE_ERRF)
		secondaryAngle = PI * 0.5 - ANGLE_ERRF;
	if (secondaryAngle < -PI * 0.5 + ANGLE_ERRF)
		secondaryAngle = -PI * 0.5 + ANGLE_ERRF;
	while (primaryAngle > PI)
		primaryAngle -= PI;
	while (primaryAngle < -PI)
		primaryAngle += PI;

	std::array<double, 6> cos = zeroDirectionCosines();
	std::array<double, 3> z = { .0, .0, 1.0 };
	vectormath::rotate(cos.data(), z.data(), primaryAngle);
	vectormath::rotate(&cos[3], z.data(), primaryAngle);
	std::array<double, 3> x = { 1.0, .0, .0 };
	vectormath::rotate(cos.data(), x.data(), -secondaryAngle);
	vectormath::rotate(&cos[3], x.data(), -secondaryAngle);
	//handling tube rotation
	std::array<double, 3> beam_dir;
	vectormath::cross(cos.data(), beam_dir.data());
	vectormath::rotate(cos.data(), beam_dir.data(), m_tubeRotationAngle);
	vectormath::rotate(&cos[3], beam_dir.data(), m_tubeRotationAngle);

	setDirectionCosines(cos);
}

void DXSource::setSourceAngles(const std::array<double, 2>& angles)
{
	setSourceAngles(angles[0], angles[1]);
}

std::array<double, 2> DXSource::sourceAngles() const
{
	auto cos = directionCosines();
	std::array<double, 3> beam_direction;
	vectormath::cross(cos.data(), beam_direction.data());
	vectormath::rotate(cos.data(), beam_direction.data(), -m_tubeRotationAngle);
	vectormath::rotate(&cos[3], beam_direction.data(), -m_tubeRotationAngle);
	vectormath::cross(cos.data(), beam_direction.data());

	//handling floating point
	double primAng, secAng;
	if (beam_direction[0] < -1.0 + ANGLE_ERRF)
	{
		primAng = PI / 2.0;
		secAng = 0.0;
	}
	else if (beam_direction[0] > 1.0 - ANGLE_ERRF)
	{
		primAng = -PI / 2.0;
		secAng = 0.0;
	}
	else
	{
		primAng = std::asin(-beam_direction[0]);
		secAng = -std::atan(beam_direction[2] / beam_direction[1]);
	}
	
	std::array<double, 2> angles = { primAng, secAng };
	
	double y[3] = { 0,1,0 };
	const auto dot_dir = vectormath::dot(y, beam_direction.data());
	if (dot_dir < 0.0) {
		if (angles[0] < 0) {
			angles[0] = -PI - angles[0];
		}
		else {
			angles[0] = PI - angles[0];
		}
	}
	
	return angles;
}

void DXSource::setSourceAnglesDeg(double primaryAngle, double secondaryAngle)
{
	setSourceAngles(primaryAngle * DEG_TO_RAD, secondaryAngle * DEG_TO_RAD);
}

void DXSource::setSourceAnglesDeg(const std::array<double, 2>& angles)
{
	setSourceAnglesDeg(angles[0], angles[1]);
}

std::array<double, 2> DXSource::sourceAnglesDeg() const
{
	auto angles = sourceAngles();
	angles[0] *= RAD_TO_DEG;
	angles[1] *= RAD_TO_DEG;
	return angles;
}

void DXSource::setTubeRotation(double angle)
{
	const double diff_angle = angle - m_tubeRotationAngle;
	std::array<double, 3> beam_direction;
	auto cos = directionCosines();
	vectormath::cross(cos.data(), beam_direction.data());
	vectormath::rotate(cos.data(), beam_direction.data(), diff_angle);
	vectormath::rotate(&cos[3], beam_direction.data(), diff_angle);
	setDirectionCosines(cos);
	m_tubeRotationAngle = angle;
}

double DXSource::tubeRotation() const
{
	return m_tubeRotationAngle;
}

void DXSource::setTubeRotationDeg(double angle)
{
	setTubeRotation(angle * DEG_TO_RAD);
}

double DXSource::tubeRotationDeg() const
{
	return tubeRotation() * RAD_TO_DEG;
}

double DXSource::getCalibrationValue(ProgressBar* progressBar) const
{
	auto specter = tube().getSpecter();
	std::vector<double> massAbsorb(specter.size(), 0.0);
	Material airMaterial("Air, Dry (near sea level)");
	std::transform(specter.begin(), specter.end(), massAbsorb.begin(), [&](const auto &el)->double {return airMaterial.getMassEnergyAbsorbtion(el.first); });
	
	const double nHist = totalExposures() * historiesPerExposure();
	
	double calcOutput = 0.0; // Air KERMA [keV/g] 
	for (std::size_t i = 0; i < specter.size(); ++i)
	{
		auto &[keV, weight] = specter[i];
		calcOutput += keV * weight * nHist * massAbsorb[i];
	}
	calcOutput *= KEV_TO_MJ * 1000.0; // Air KERMA [mJ / kg] = [mGy]

	const double output = m_dap / (m_fieldSize[0] * m_fieldSize[1] * 0.01); //mm->cm
	const double factor = output / calcOutput; // mGy/mGy
	return factor; 
}

void DXSource::updateCollimationAngles(const std::array<double, 2>& collimationAngles)
{
	for (std::size_t i = 0; i < 2; ++i)
	{
		m_collimationAngles[i] = collimationAngles[i];
		m_fieldSize[i] = std::tan(m_collimationAngles[i] * 0.5) * m_sdd * 2.0;
	}
	m_specterValid = false;
}
void DXSource::updateFieldSize(const std::array<double, 2>& fieldSize)
{
	for (std::size_t i = 0; i < 2; ++i)
	{
		m_fieldSize[i] = fieldSize[i];
		m_collimationAngles[i] = std::atan(m_fieldSize[i] * 0.5 / m_sdd) * 2.0;
	}
	m_specterValid = false;
}

const std::array<double, 3> DXSource::tubePosition(void) const
{
	std::array<double, 3> beamDirection;
	vectormath::cross(m_directionCosines.data(), beamDirection.data());
	std::array<double, 3> pos;
	for (std::size_t i = 0; i < 3; ++i)
	{
		pos[i] = m_position[i] - beamDirection[i] * m_sdd;
	}
	return pos;
}

CTSource::CTSource()
	:Source()
{
	m_type = Type::None;
	m_sdd = 1190.0;
	m_collimation = 38.4;
	m_fov = 500.0;
    m_startAngle = 0.0;
	m_exposureAngleStep = DEG_TO_RAD;
    m_scanLenght = 100.0;
	auto &t = tube();
	t.setAlFiltration(7.0);
	const std::array<double, 6> ct_cosines{ -1,0,0,0,0,1 };
	setDirectionCosines(ct_cosines);

}


void CTSource::setSourceDetectorDistance(double sdd)
{
	m_sdd = std::abs(sdd);
	m_specterValid = false;
}

double CTSource::sourceDetectorDistance(void) const
{
	return m_sdd;
}

void CTSource::setCollimation(double collimation)
{
	m_collimation = std::abs(collimation);
	m_specterValid = false;
}

double CTSource::collimation(void) const
{
	return m_collimation;
}

void CTSource::setFieldOfView(double fov)
{
	m_fov = std::abs(fov);
}

double CTSource::fieldOfView(void) const
{
	return m_fov;
}

void CTSource::setGantryTiltAngle(double angle)
{
	if (angle < -PI)
		m_gantryTiltAngle = -PI;
	else if (angle > PI)
		m_gantryTiltAngle = PI;
	else
		m_gantryTiltAngle = angle;
}

double CTSource::gantryTiltAngle() const
{
	return m_gantryTiltAngle;
}

void CTSource::setGantryTiltAngleDeg(double angle)
{
	setGantryTiltAngle(angle * DEG_TO_RAD);
}

double CTSource::gantryTiltAngleDeg() const
{
	return m_gantryTiltAngle * RAD_TO_DEG;
}

void CTSource::setStartAngle(double angle)
{
    m_startAngle = angle;
}

double CTSource::startAngle(void) const
{
    return m_startAngle;
}
void CTSource::setStartAngleDeg(double angle)
{
	m_startAngle = DEG_TO_RAD * angle;
}

double CTSource::startAngleDeg(void) const
{
	return RAD_TO_DEG * m_startAngle;
}

void CTSource::setExposureAngleStepDeg(double angleStep)
{
	setExposureAngleStep(angleStep * DEG_TO_RAD);
}

double CTSource::exposureAngleStepDeg(void) const
{
	return m_exposureAngleStep * RAD_TO_DEG;
}

void CTSource::setExposureAngleStep(double angleStep)
{
	const double absAngle = std::abs(angleStep);
	if (absAngle < PI)
		m_exposureAngleStep = absAngle;
	else if (absAngle > 0.1 * DEG_TO_RAD)
		m_exposureAngleStep = absAngle;
}

double CTSource::exposureAngleStep(void) const
{
	return m_exposureAngleStep;
}


void CTSource::setScanLenght(double scanLenght)
{
	m_scanLenght = std::abs(scanLenght);
}
double CTSource::scanLenght(void) const
{
	return m_scanLenght;
}


double ctdiStatIndex(const std::array<double, 5>& measurements)
{
	double mean = 0.0;
	for (std::size_t i = 1; i < 5; ++i)
		mean += measurements[i];
	mean /= 4.0;
	if (mean <= 0.0)
		return 1000.0;
	double stddev = 0;
	for (std::size_t i = 1; i < 5; ++i)
		stddev += (measurements[i] - mean) * (measurements[i] - mean);
	stddev = std::sqrt(stddev / 3.0);

	return stddev / mean;
}

void CTSource::updateSpecterDistribution()
{
	if (!m_specterValid)
	{
		auto energies = m_tube.getEnergy();
		auto n_obs = m_tube.getSpecter(energies);
		m_specterDistribution = std::make_unique<SpecterDistribution>(n_obs, energies);
		const double heel_span_angle = std::atan(m_collimation * 0.5 / m_sdd) * 2.0;
		if (m_modelHeelEffect)
			m_heelFilter = std::make_unique<HeelFilter>(m_tube, heel_span_angle);
		else
			m_heelFilter = nullptr;
		m_specterValid = true;
	}
}

//double CTSource::getCalibrationValue(ProgressBar* progressBar) const

template<typename T>
static double CTSource::ctCalibration(T& sourceCopy, ProgressBar* progressBar)
{
	CTDIPhantom world(sourceCopy.ctdiPhantomDiameter());
	world.setAttenuationLutMaxEnergy(sourceCopy.tube().voltage());
	world.validate();

	sourceCopy.updateFromWorld(world);

	double meanWeight = 0;
	for (std::size_t i = 0; i < sourceCopy.totalExposures(); ++i)
	{
		Exposure dummy;
		sourceCopy.getExposure(dummy, i);
		meanWeight += dummy.beamIntensityWeight();
	}
	meanWeight /= static_cast<double>(sourceCopy.totalExposures());

	const std::array<double, 6> cosines({ -1,0,0,0,0,1 });
	sourceCopy.setDirectionCosines(cosines);

	sourceCopy.setUseXCareFilter(false); // we need to disable organ aec for ctdi statistics, this should be ok 
	std::size_t statCounter = CTDI_MIN_HISTORIES / (sourceCopy.exposuresPerRotatition() * sourceCopy.historiesPerExposure());
	if (statCounter < 1)
		statCounter = 1;
	const auto histories = sourceCopy.historiesPerExposure();
	sourceCopy.setHistoriesPerExposure(histories * statCounter); // ensuring enough histories for ctdi measurement
	sourceCopy.validate();
	
	auto result = transport::run(world, &sourceCopy, progressBar);

	typedef CTDIPhantom::HolePosition holePosition;
	std::array<CTDIPhantom::HolePosition, 5> position = { holePosition::Center, holePosition::West, holePosition::East, holePosition::South, holePosition::North };

	std::array<double, 5> measureDose;
	measureDose.fill(0.0);
	for (std::size_t i = 0; i < 5; ++i)
	{
		auto holeIndices = world.holeIndices(position[i]);
		for (auto idx : holeIndices)
			measureDose[i] += result.dose[idx];
		measureDose[i] /= static_cast<double>(holeIndices.size());
	}
	
	
	//test
	std::array<double, 5> measureTally;
	measureTally.fill(0.0);
	const double totalHistories = sourceCopy.exposuresPerRotatition() * sourceCopy.historiesPerExposure();
	for (std::size_t i = 0; i < 5; ++i)
	{
		auto holeIndices = world.holeIndices(position[i]);
		double tally = 0;
		for (auto idx : holeIndices)
			tally += result.nEvents[idx];
		measureTally[i] = std::sqrt(tally / (totalHistories * totalHistories) - tally * tally / (totalHistories * totalHistories * totalHistories));
	}
	//test end
	


	const double ctdiPher = (measureDose[1] + measureDose[2] + measureDose[3] + measureDose[4]) / 4.0;
	const double ctdiCent = measureDose[0];
	const double ctdiw = (ctdiCent + 2.0 * ctdiPher) / 3.0 / static_cast<double>(statCounter);
	const double factor = sourceCopy.ctdiVol() / ctdiw / meanWeight;
	return factor;
}

std::uint64_t CTSource::exposuresPerRotatition() const 
{
	return static_cast<std::size_t>(PI_2 / m_exposureAngleStep);
}


CTSpiralSource::CTSpiralSource()
	:CTSource()
{
	m_type = Type::CTSpiral;
	m_pitch = 1.0;
}
void CTSpiralSource::setPitch(double pitch)
{
	m_pitch = std::max(0.01, pitch);
}

double CTSpiralSource::pitch(void) const
{
	return m_pitch;
}
bool CTSpiralSource::getExposure(Exposure& exposure, std::uint64_t exposureIndex) const
{
	std::array<double, 3> pos = { 0, -m_sdd / 2.0,0 };

	const double angle = m_startAngle + m_exposureAngleStep * exposureIndex;

	std::array<double, 3> rotationAxis, otherAxis;
	for (std::size_t i = 0; i < 3; ++i)
	{
		rotationAxis[i] = m_directionCosines[i + 3];
		otherAxis[i] = m_directionCosines[i];
	}
	std::array<double, 3> tiltAxis = { 1,0,0 };
	auto tiltCorrection = pos;
	vectormath::rotate(tiltCorrection.data(), tiltAxis.data(), m_gantryTiltAngle);
	vectormath::rotate(rotationAxis.data(), tiltAxis.data(), m_gantryTiltAngle);
	vectormath::rotate(otherAxis.data(), tiltAxis.data(), m_gantryTiltAngle);
	vectormath::rotate(pos.data(), rotationAxis.data(), angle);
	
	//ading transverse step
	pos[2] += (exposureIndex * m_exposureAngleStep) * m_collimation * m_pitch / PI_2 + tiltCorrection[2];

	vectormath::rotate(otherAxis.data(), rotationAxis.data(), angle);
	for (std::size_t i = 0; i < 3; ++i)
		pos[i] += m_position[i];

	exposure.setPosition(pos);
	exposure.setDirectionCosines(otherAxis, rotationAxis);
	exposure.setCollimationAngles(std::atan(m_fov / m_sdd) * 2.0, std::atan(m_collimation / m_sdd) * 2.0);
	exposure.setBeamFilter(m_bowTieFilter.get());
	exposure.setSpecterDistribution(m_specterDistribution.get());
	exposure.setHeelFilter(m_heelFilter.get());
	exposure.setNumberOfHistories(m_historiesPerExposure);
	double weight = 1.0;
	if (m_positionalFilter)
		weight *= m_positionalFilter->sampleIntensityWeight(pos);
	if (m_useXCareFilter)
		weight *= m_xcareFilter.sampleIntensityWeight(angle);
	exposure.setBeamIntensityWeight(weight);
	return true;
}


std::uint64_t CTSpiralSource::totalExposures() const
{
	return static_cast<std::uint64_t>(m_scanLenght * PI_2 / (m_collimation * m_pitch * m_exposureAngleStep));
}

double CTSpiralSource::getCalibrationValue(ProgressBar* progressBar) const
{
	auto copy = *this;
	return ctCalibration(copy, progressBar) * m_pitch;
	//return CTSource::getCalibrationValue(progressBar) / m_pitch;
}


CTAxialSource::CTAxialSource()
	:CTSource()
{
	m_type = Type::CTAxial;
	m_step = m_collimation;
	m_scanLenght = m_step;
}
void CTAxialSource::setStep(double step)
{
	const double absStep = std::abs(step);
	auto n_steps = m_scanLenght / m_step;
	m_step = absStep > 0.01 ? absStep : 0.01;
	setScanLenght(m_step*n_steps);
}

double CTAxialSource::step(void) const
{
	return m_step;
}

void CTAxialSource::setScanLenght(double scanLenght)
{
	m_scanLenght = m_step * std::floor(scanLenght / m_step);
}

bool CTAxialSource::getExposure(Exposure& exposure, std::uint64_t exposureIndex) const
{
	//calculating position
	std::array<double, 3> pos = { 0,-m_sdd / 2.0,0 };
	
	const std::uint64_t anglesPerRotation = static_cast<std::uint64_t>(PI_2 / m_exposureAngleStep);
	const std::uint64_t rotationNumber = exposureIndex / anglesPerRotation;

	const double angle = m_startAngle + m_exposureAngleStep * (exposureIndex - (rotationNumber*anglesPerRotation));

	std::array<double, 3> rotationAxis, otherAxis;
	for (std::size_t i = 0; i < 3; ++i)
	{
		rotationAxis[i] = m_directionCosines[i + 3];
		otherAxis[i] = m_directionCosines[i];
	}
	std::array<double, 3> tiltAxis = { 1,0,0 };
	auto tiltCorrection = pos;
	vectormath::rotate(tiltCorrection.data(), tiltAxis.data(), m_gantryTiltAngle);
	vectormath::rotate(rotationAxis.data(), tiltAxis.data(), m_gantryTiltAngle);
	vectormath::rotate(otherAxis.data(), tiltAxis.data(), m_gantryTiltAngle);
	vectormath::rotate(pos.data(), rotationAxis.data(), angle);
	pos[2] += m_step * rotationNumber + tiltCorrection[2];
	vectormath::rotate(otherAxis.data(), rotationAxis.data(), angle);
	for (std::size_t i = 0; i < 3; ++i)
		pos[i] += m_position[i];

	exposure.setPosition(pos);
	exposure.setDirectionCosines(otherAxis, rotationAxis);
	exposure.setCollimationAngles(std::atan(m_fov / m_sdd) * 2.0, std::atan(m_collimation / m_sdd) * 2.0);
	exposure.setBeamFilter(m_bowTieFilter.get());
	exposure.setHeelFilter(m_heelFilter.get());
	exposure.setSpecterDistribution(m_specterDistribution.get());
	exposure.setNumberOfHistories(m_historiesPerExposure);
	double weight = 1.0;
	if (m_positionalFilter)
		weight *= m_positionalFilter->sampleIntensityWeight(pos);
	if (m_useXCareFilter)
		weight *= m_xcareFilter.sampleIntensityWeight(angle);
	exposure.setBeamIntensityWeight(weight);
	return true;
}


std::uint64_t CTAxialSource::totalExposures() const
{
	const std::uint64_t anglesPerRotation = static_cast<std::uint64_t>(PI_2 / m_exposureAngleStep);
	const std::uint64_t rotationNumbers = static_cast<std::uint64_t>(std::round(m_scanLenght / m_step));
	return anglesPerRotation * (rotationNumbers + 1);
}

double CTAxialSource::getCalibrationValue(ProgressBar* progressBar) const
{
	auto copy = *this;
	return ctCalibration(copy, progressBar);
}

CTDualSource::CTDualSource()
	:CTSource()
{
	m_type = Type::CTDual;
	
	m_sddB = m_sdd;
	m_fovB = m_fov;
	
	m_startAngleB = m_startAngle + PI * 0.5;
	
	m_pitch = 1.0;
	m_tube.setAlFiltration(7.0);
	m_tubeB.setAlFiltration(7.0);
}


bool CTDualSource::getExposure(Exposure& exposure, std::uint64_t exposureIndexTotal) const
{
	double sdd, startAngle, fov;
	uint64_t exposureIndex = exposureIndexTotal / 2;
	BeamFilter* bowTie=nullptr;
	SpecterDistribution* specterDistribution=nullptr;
	HeelFilter* heelFilter = nullptr;
	if (exposureIndexTotal % 2 == 0)
	{
		sdd = m_sdd;
		startAngle = m_startAngle;
		fov = m_fov;
		bowTie = m_bowTieFilter.get();
		specterDistribution = m_specterDistribution.get();
		heelFilter = m_heelFilter.get();
	}
	else
	{
		sdd = m_sddB;
		startAngle = m_startAngleB;
		fov = m_fovB;
		bowTie = m_bowTieFilterB.get();
		specterDistribution = m_specterDistributionB.get();
		heelFilter = m_heelFilterB.get();
	}
	std::array<double, 3> pos = { 0,-m_sdd / 2.0,0 };

	const double angle = startAngle + m_exposureAngleStep * exposureIndex;

	std::array<double, 3> rotationAxis, otherAxis;
	for (std::size_t i = 0; i < 3; ++i)
	{
		rotationAxis[i] = m_directionCosines[i + 3];
		otherAxis[i] = m_directionCosines[i];
	}
	
	std::array<double, 3> tiltAxis = { 1,0,0 };
	auto tiltCorrection = pos;
	vectormath::rotate(tiltCorrection.data(), tiltAxis.data(), m_gantryTiltAngle);
	vectormath::rotate(rotationAxis.data(), tiltAxis.data(), m_gantryTiltAngle);
	vectormath::rotate(otherAxis.data(), tiltAxis.data(), m_gantryTiltAngle);
	vectormath::rotate(pos.data(), rotationAxis.data(), angle);
	pos[2] += (exposureIndex * m_exposureAngleStep) * m_collimation * m_pitch / PI_2 + tiltCorrection[2];
	
	vectormath::rotate(otherAxis.data(), rotationAxis.data(), angle);
	for (std::size_t i = 0; i < 3; ++i)
		pos[i] += m_position[i];

	exposure.setPosition(pos);
	exposure.setDirectionCosines(otherAxis, rotationAxis);
	exposure.setCollimationAngles(std::atan(fov / sdd) * 2.0, std::atan(m_collimation / sdd) * 2.0);
	exposure.setBeamFilter(bowTie);
	exposure.setSpecterDistribution(specterDistribution);
	exposure.setHeelFilter(heelFilter);
	exposure.setNumberOfHistories(m_historiesPerExposure);
	double weight = 1.0;
	if (m_positionalFilter)
		weight *= m_positionalFilter->sampleIntensityWeight(pos);
	if (m_useXCareFilter)
		weight *= m_xcareFilter.sampleIntensityWeight(angle);
	exposure.setBeamIntensityWeight(weight);
	return true;
}

std::uint64_t CTDualSource::totalExposures(void) const
{

	auto singleSourceExposure =  static_cast<std::uint64_t>(m_scanLenght * PI_2 / (m_collimation * m_pitch * m_exposureAngleStep));
	return singleSourceExposure * 2;
}

std::uint64_t CTDualSource::exposuresPerRotatition() const
{
	return 2 * static_cast<std::size_t>(PI_2 / m_exposureAngleStep);
}

double CTDualSource::getCalibrationValue(ProgressBar* progressBar) const
{
	auto copy = *this;
	return ctCalibration(copy, progressBar) * m_pitch;
}

void CTDualSource::setStartAngleDegB(double angle)
{
	m_startAngleB = DEG_TO_RAD * angle;
}

double CTDualSource::startAngleDegB(void) const
{
	return RAD_TO_DEG * m_startAngleB;
}

void CTDualSource::setPitch(double pitch)
{
	m_pitch = std::max(0.01, pitch);
}

double CTDualSource::pitch(void) const
{
	return m_pitch;
}


void CTDualSource::updateSpecterDistribution()
{
	if (!m_specterValid)
	{
		auto energyA = m_tube.getEnergy();
		auto energyB = m_tubeB.getEnergy();
		auto specterA = m_tube.getSpecter(energyA, false);
		auto specterB = m_tubeB.getSpecter(energyB, false);

		double sumA = std::accumulate(specterA.begin(), specterA.end(), 0.0);
		double sumB = std::accumulate(specterB.begin(), specterB.end(), 0.0);
		double weightA = m_tubeAmas * sumA;
		double weightB = m_tubeBmas * sumB;

		std::for_each(specterA.begin(), specterA.end(), [=](double& val) {val /= sumA; });
		std::for_each(specterB.begin(), specterB.end(), [=](double& val) {val /= sumB; });

		m_tubeAweight = 1.0;
		m_tubeBweight = weightB / weightA;

		m_specterDistribution = std::make_unique<SpecterDistribution>(specterA, energyA);
		m_specterDistributionB = std::make_unique<SpecterDistribution>(specterB, energyB);

		const double heel_span_angle = std::atan(m_collimation * 0.5 / m_sdd) * 2.0;
		m_heelFilter = std::make_unique<HeelFilter>(m_tube, heel_span_angle);
		m_heelFilterB = std::make_unique<HeelFilter>(m_tubeB, heel_span_angle);

		m_specterValid = true;
	}
}

