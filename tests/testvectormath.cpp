
#include "dxmc/source.h"
#include "dxmc/vectormath.h"
#include <iostream>


constexpr double RAD2DEG = 180.0 / 3.14159265359;
constexpr double DEG2RAD = 1.0 / RAD2DEG;

bool testSource()
{
	std::cout << "\ntestSource:\n";
	DXSource src;

	std::array<double, 2> angles = { 45,45 };
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
	std::array<double, 2> angles = { 45,45 };
	std::cout << "rotation angles: " << angles[0] << ", " << angles[1] <<'\n';
	vectormath::rotate(v, x, angles[1] * DEG2RAD);
	vectormath::rotate(v, z, angles[0] * DEG2RAD);
	vectormath::normalize(v);
	std::cout << "vector rotated: " << v[0] << ", " << v[1] << ", " << v[2] << '\n';
	
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
	auto ss = testSource();
}
