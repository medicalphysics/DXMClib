
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


int main(int argc, char* argv[])
{
	testArgMax();
	return EXIT_SUCCESS;
}
