
#include "dxmc/source.h"
#include "dxmc/vectormath.h"
#include <iostream>


constexpr double RAD2DEG = 180.0 / 3.14159265359;
constexpr double DEG2RAD = 1.0 / RAD2DEG;

bool testSource(double pang, double sang)
{
	std::cout << "\ntestSource:\n";
	DXSource src;

	std::array<double, 2> angles = { pang,sang };
	src.setSourceAnglesDeg(angles);
	auto anglesres = src.sourceAnglesDeg();
	std::cout << "angles set: " << angles[0] << ", " << angles[1] << '\n';
	std::cout << "angles res: " << anglesres[0] << ", " << anglesres[1] << '\n';



	return true;

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


	double axis[3] = { std::sin(angles[0]), std::cos(angles[0]), std::cos(-angles[1]) };
	vectormath::normalize(axis);
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
	auto success = testVectorAngles();

	std::array<double, 5> angles = { -60, -30, 0, 30, 60 };
	for (auto ap : angles)
		for (auto as : angles)
			testSource(ap, as);

}
