#include "dxmc/source.h"
#include "dxmc/vectormath.h"
#include <iostream>
#include <cassert>

constexpr double RAD2DEG = 180.0 / 3.14159265359;
constexpr double DEG2RAD = 1.0 / RAD2DEG;
constexpr double ERRF = 1e-4;

bool isEqual(double f1, double f2)
{
	return std::abs(f1 - f2) < ERRF;
}


bool testIsotropicSourceSpecter()
{
	std::vector<double> specArr = { 16.25,1.423E-04,16.75,2.157E-04,17.25,3.102E-04,17.75,4.324E-04,18.25,5.840E-04,18.75,7.644E-04,19.25,9.784E-04,19.75,1.222E-03,20.25,1.491E-03,20.75,1.803E-03,21.25,2.129E-03,21.75,2.490E-03,22.25,2.863E-03,22.75,3.263E-03,23.25,3.658E-03,23.75,4.093E-03,24.25,4.504E-03,24.75,4.912E-03,25.25,5.347E-03,25.75,5.769E-03,26.25,6.168E-03,26.75,6.582E-03,27.25,6.965E-03,27.75,7.360E-03,28.25,7.710E-03,28.75,8.067E-03,29.25,8.368E-03,29.75,8.671E-03,30.25,8.975E-03,30.75,9.213E-03,31.25,9.476E-03,31.75,9.694E-03,32.25,9.903E-03,32.75,1.009E-02,33.25,1.025E-02,33.75,1.040E-02,34.25,1.053E-02,34.75,1.063E-02,35.25,1.073E-02,35.75,1.081E-02,36.25,1.087E-02,36.75,1.092E-02,37.25,1.096E-02,37.75,1.099E-02,38.25,1.100E-02,38.75,1.100E-02,39.25,1.099E-02,39.75,1.098E-02,40.25,1.095E-02,40.75,1.091E-02,41.25,1.086E-02,41.75,1.081E-02,42.25,1.076E-02,42.75,1.069E-02,43.25,1.063E-02,43.75,1.055E-02,44.25,1.048E-02,44.75,1.039E-02,45.25,1.031E-02,45.75,1.022E-02,46.25,1.012E-02,46.75,1.003E-02,47.25,9.933E-03,47.75,9.828E-03,48.25,9.732E-03,48.75,9.628E-03,49.25,9.516E-03,49.75,9.412E-03,50.25,9.302E-03,50.75,9.193E-03,51.25,9.084E-03,51.75,8.970E-03,52.25,8.862E-03,52.75,8.749E-03,53.25,8.637E-03,53.75,8.526E-03,54.25,8.409E-03,54.75,8.300E-03,55.25,8.185E-03,55.75,8.072E-03,56.25,7.959E-03,56.75,7.847E-03,57.25,7.737E-03,57.75,2.568E-02,58.25,7.513E-03,58.75,7.405E-03,59.25,3.920E-02,59.75,7.181E-03,60.25,7.071E-03,60.75,6.962E-03,61.25,6.854E-03,61.75,6.746E-03,62.25,6.640E-03,62.75,6.530E-03,63.25,6.425E-03,63.75,6.321E-03,64.25,6.214E-03,64.75,6.107E-03,65.25,6.006E-03,65.75,5.901E-03,66.25,5.797E-03,66.75,1.673E-02,67.25,5.592E-03,67.75,5.491E-03,68.25,5.390E-03,68.75,8.223E-03,69.25,5.055E-03,69.75,4.296E-03,70.25,4.236E-03,70.75,4.171E-03,71.25,4.110E-03,71.75,4.048E-03,72.25,3.982E-03,72.75,3.919E-03,73.25,3.852E-03,73.75,3.787E-03,74.25,3.719E-03,74.75,3.654E-03,75.25,3.585E-03,75.75,3.516E-03,76.25,3.449E-03,76.75,3.379E-03,77.25,3.308E-03,77.75,3.240E-03,78.25,3.169E-03,78.75,3.098E-03,79.25,3.026E-03,79.75,2.954E-03,80.25,2.882E-03,80.75,2.809E-03,81.25,2.736E-03,81.75,2.665E-03,82.25,2.592E-03,82.75,2.519E-03,83.25,2.445E-03,83.75,2.370E-03,84.25,2.296E-03,84.75,2.222E-03,85.25,2.148E-03,85.75,2.073E-03,86.25,1.999E-03,86.75,1.925E-03,87.25,1.850E-03,87.75,1.776E-03,88.25,1.700E-03,88.75,1.625E-03,89.25,1.550E-03,89.75,1.476E-03,90.25,1.400E-03,90.75,1.326E-03,91.25,1.251E-03,91.75,1.177E-03,92.25,1.101E-03,92.75,1.027E-03,93.25,9.529E-04,93.75,8.781E-04,94.25,8.041E-04,94.75,7.302E-04,95.25,6.559E-04,95.75,5.823E-04,96.25,5.089E-04,96.75,4.353E-04,97.25,3.623E-04,97.75,2.892E-04,98.25,2.166E-04,98.75,1.441E-04,99.25,7.193E-05,99.75,5.990E-06 };
	std::vector<double> energies, weighs;
	for (std::size_t i = 0; i < specArr.size(); i = i + 2)
	{
		energies.push_back(specArr[i]);
		weighs.push_back(specArr[i + 1]);
	}

	//normalize weights
	const auto weights_sum = std::accumulate(weighs.cbegin(), weighs.cend(), 0.0);
	std::transform(weighs.begin(), weighs.end(), weighs.begin(), [=](double h)->double {return h / weights_sum; });

	IsotropicSource src;
	src.setSpecter(weighs, energies);
	Exposure exp;
	src.getExposure(exp, 0);
	
	std::vector<double> hist(energies.size(), 0.0);
	const double emin = energies[0];
	const double estep = energies[1] - energies[0];
	Particle particle;
	std::uint64_t seed[2];
	randomSeed(seed);
	for (std::size_t i = 0; i < 1e7; ++i)
	{
		exp.sampleParticle(particle, seed);
		assert(isEqual(vectormath::lenght_sqr(particle.dir), 1.0));
		assert(isEqual(particle.weight, 1.0));
		std::size_t eidx = static_cast<std::size_t>((particle.energy - emin + estep * 0.5) / estep);
		hist[eidx] += particle.weight;
	}

	//normalize hist
	const auto hist_sum = std::accumulate(hist.cbegin(), hist.cend(), 0.0);
	std::transform(hist.begin(), hist.end(), hist.begin(), [=](double h)->double {return h / hist_sum; });
	
	std::cout << "Specter, Histogram\n";
	bool success = true;
	for (int i = 0; i < hist.size(); ++i)
	{
		std::cout << weighs[i] << ", " << hist[i] << "\n";
		success = success && isEqual(weighs[i], hist[i]);
		if (!success)
		{
			auto t_weights = weighs[i];
			auto t_hist = hist[i];
			success = false;
		}
	}

	assert(success);
	return success;
}


bool testSourceAngles(double pang, double sang, double tubeRotation)
{

	DXSource src;

	std::array<double, 2> angles = { pang,sang };
	src.setTubeRotationDeg(tubeRotation);
	src.setSourceAnglesDeg(angles);
	auto anglesres = src.sourceAnglesDeg();
	std::cout << "angles set: " << angles[0] << ", " << angles[1];
	std::cout << " angles res: " << anglesres[0] << ", " << anglesres[1];
	std::cout << " tube rot: " << tubeRotation << '\n';
	return isEqual(angles[0], anglesres[0]) && isEqual(angles[1], anglesres[1]);
}

bool testSourceAnglesMany()
{
	bool success = true;

	testSourceAngles(90, 0, 45);
	testSourceAngles(90, 90, 45);
	testSourceAngles(90, 90, 90);

	std::array<double, 7> angles = { -89, -60,-30,0,30,60, 89 };
	auto tube_rot = angles;
	for (auto ap : angles)
		for (auto as : angles)
			for (auto tr : tube_rot)
			{
				success = success && testSourceAngles(ap, as, tr);
			}

	success = success && testSourceAngles(30, 30, 0);
	success = success && testSourceAngles(150, 30, 0);
	success = success && testSourceAngles(-150, 30, 0);
	success = success && testSourceAngles(-150, 30, 90);
	success = success && testSourceAngles(-150, 30, 180);

	assert(success);
	return success;
}

bool testCTCalibration()
{
	CTSpiralSource src;
	src.setPitch(0.5);
	src.setExposureAngleStepDeg(5.0);
	src.setHistoriesPerExposure(100000);
	auto factor = src.getCalibrationValue();


	return false;
}

int main(int argc, char* argv[])
{
	bool success = true;
	success = success && testCTCalibration();
	success = success && testSourceAnglesMany();
	success = success && testIsotropicSourceSpecter();
	return !success;
}
