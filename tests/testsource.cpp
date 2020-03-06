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



int main(int argc, char* argv[])
{
	bool success = true;
	success = success && testSourceAnglesMany();
	return !success;
}
