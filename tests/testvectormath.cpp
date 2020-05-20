
#include "dxmc/source.h"
#include "dxmc/vectormath.h"
#include <iostream>
#include <cassert>


void testArgMax()
{
	std::array<double, 3> arr{ 1,2,3 };
	assert((vectormath::argmax3<std::size_t, double>(arr.data()) == 2));
	assert((vectormath::argmin3<std::size_t, double>(arr.data()) == 0));
	arr[0] = 2; arr[1] = 1; arr[2] = 3;
	assert((vectormath::argmax3<std::size_t, double>(arr.data()) == 2));
	assert((vectormath::argmin3<std::size_t, double>(arr.data()) == 1));
	arr[0] = 3; arr[1] = 2; arr[2] = 1;
	assert((vectormath::argmax3<std::size_t, double>(arr.data()) == 0));
	assert((vectormath::argmin3<std::size_t, double>(arr.data()) == 2));
	return;
}

void testChangeBasis()
{
	std::array<double, 6> cos{ 1,0,0,0,1,0 };
	std::array<double, 3> dir;
	vectormath::cross(cos.data(), dir.data());

	std::array<double, 3> pos{ 1,0,0 };
	auto posFin = pos;

	const double* b1 = &cos[0];
	const double* b2 = &cos[3];
	const double* b3 = &dir[0];

	vectormath::changeBasisInverse(b1, b2, b3, posFin.data());

	for (std::size_t i = 0; i < 3; ++i)
		assert((pos[i] == posFin[i]));
}

int main(int argc, char* argv[])
{
	testArgMax();
	testChangeBasis();
	return EXIT_SUCCESS;
}
