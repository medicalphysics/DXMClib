
#include "dxmc/source.h"
#include "dxmc/vectormath.h"
#include <iostream>
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
	std::cout << "angles set: " << angles[0] << ", " << angles[1];
	std::cout << " angles res: " << anglesres[0] << ", " << anglesres[1] << '\n';

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



bool testVectorAngles()
{
	std::cout << "\ntestVectorAngles:\n";
	double v[3] = { 0,1,0 };
	vectormath::normalize(v);

	double x[3] = { 1,0,0 };
	double y[3] = { 0,1,0 };
	double z[3] = { 0,0,1 };

	std::cout << "vector: " << v[0] << ", " << v[1] << ", " << v[2] << '\n';
	std::array<double, 2> angles = { 30*DEG2RAD,60*DEG2RAD };
	std::cout << "rotation angles: " << angles[0]*RAD2DEG << ", " << angles[1]*RAD2DEG <<'\n';

	
	//double axis[3] = { std::sin(angles[0]), std::cos(angles[0])*std::sin(angles[1]), std::sin(-angles[1]) };
	double axis[3] = { 0,1,0 };
	vectormath::rotate(axis, z, angles[0]);
	vectormath::rotate(axis, x, angles[1]);
	std::cout << "rotation angles calc: " << std::sin(-axis[0])* RAD2DEG << ", " << std::atan(axis[2]/axis[1]) * RAD2DEG << '\n';



	//vectormath::normalize(axis);
	double rotAxis[3];
	vectormath::cross(v, axis, rotAxis);
	vectormath::normalize(rotAxis);

	double angle = std::acos(vectormath::dot(v, axis));

	vectormath::rotate(v, rotAxis, angle);

	
	//vectormath::normalize(v);
	std::cout << "vector rotated: " << v[0] << ", " << v[1] << ", " << v[2] << '\n';
	std::cout << "vector lenght: " << vectormath::lenght(v) << '\n';
	
	std::cout << "angle between v and x: " << vectormath::angleBetween(v, x) * RAD2DEG << '\n';
	std::cout << "angle between v and y: " << vectormath::angleBetween(v, y) * RAD2DEG << '\n';
	std::cout << "angle between v and z: " << vectormath::angleBetween(v, z) * RAD2DEG << '\n';

	std::cout << "angle between v and x on x: " << vectormath::angleBetweenOnPlane(v, x, x) * RAD2DEG << '\n';
	std::cout << "angle between v and y on x: " << vectormath::angleBetweenOnPlane(v, y, x) * RAD2DEG << '\n';
	std::cout << "angle between v and z on x: " << vectormath::angleBetweenOnPlane(v, z, x) * RAD2DEG << '\n';
	std::cout << "angle between v and x on y: " << vectormath::angleBetweenOnPlane(v, x, y) * RAD2DEG << '\n';
	std::cout << "angle between v and y on y: " << vectormath::angleBetweenOnPlane(v, y, y) * RAD2DEG << '\n';
	std::cout << "angle between v and z on y: " << vectormath::angleBetweenOnPlane(v, z, y) * RAD2DEG << '\n'; 
	std::cout << "angle between v and x on z: " << vectormath::angleBetweenOnPlane(v, x, z) * RAD2DEG << '\n';
	std::cout << "angle between v and y on z: " << vectormath::angleBetweenOnPlane(v, y, z) * RAD2DEG << '\n';
	std::cout << "angle between v and z on z: " << vectormath::angleBetweenOnPlane(v, z, z) * RAD2DEG << '\n';
	return true;

}


int main(int argc, char* argv[])
{
	bool success = true;
	//success = success && testVectorAngles();
	std::cout << "\ntestSource:\n";

	success = success && testSourceAnglesMany();
	return success;
}
