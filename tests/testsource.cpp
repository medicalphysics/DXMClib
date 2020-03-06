#include "dxmc/source.h"
#include "dxmc/vectormath.h"
//#include <iostream>
#include <cassert>

constexpr double RAD2DEG = 180.0 / 3.14159265359;
constexpr double DEG2RAD = 1.0 / RAD2DEG;
constexpr double ERRF = 1E-5;

bool isEqual(double f1, double f2)
{
	return std::abs(f1 - f2) < ERRF;
}

bool testSourceAngles(double pang, double sang)
{

	DXSource src;

	std::array<double, 2> angles = { pang,sang };
	src.setSourceAnglesDeg(angles);
	auto anglesres = src.sourceAnglesDeg();
	//std::cout << "angles set: " << angles[0] << ", " << angles[1];
	//std::cout << " angles res: " << anglesres[0] << ", " << anglesres[1] << '\n';

	return isEqual(angles[0], anglesres[0]) && isEqual(angles[1], anglesres[1]);
}

bool testSourceAnglesMany()
{
	bool success = true;

	success = success && testSourceAngles(30, 30);
	success = success && testSourceAngles(150, 30);
	success = success && testSourceAngles(-150, 30);

	std::array<double, 7> angles = { -90, -60,-30,0,30,60, 90 };
	for (auto ap : angles)
		for (auto as : angles)
			success = success && testSourceAngles(ap, as);

	assert(success);
	return success;
}



int main(int argc, char* argv[])
{
	bool success = true;
	success = success && testSourceAnglesMany();
	return !success;
}
