
#include "dxmc/source.hpp"
#include "dxmc/vectormath.hpp"

#include <cassert>
#include <iostream>

using namespace dxmc;

constexpr double ERRF = 1e-4;

template <typename T>
bool isEqual(T f1, T f2)
{
    return std::abs(f1 - f2) < ERRF;
}

template <typename T>
void testLenght()
{
    
    std::array<T, 3> arr_std = { 1, 1, 1 };

    auto l = vectormath::lenght(arr_std);
    

    assert(isEqual(l, std::sqrt(T { 3 })));
}

void testArgMax()
{
    std::array<double, 3> arr { 1, 2, 3 };
    assert((vectormath::argmax3<std::size_t, double>(arr) == 2));
    assert((vectormath::argmin3<std::size_t, double>(arr) == 0));
    arr[0] = 2;
    arr[1] = 1;
    arr[2] = 3;
    assert((vectormath::argmax3<std::size_t, double>(arr) == 2));
    assert((vectormath::argmin3<std::size_t, double>(arr) == 1));
    arr[0] = 3;
    arr[1] = 2;
    arr[2] = 1;
    assert((vectormath::argmax3<std::size_t, double>(arr) == 0));
    assert((vectormath::argmin3<std::size_t, double>(arr) == 2));
    return;
}

void testChangeBasis()
{
    std::array<double, 6> cos { 1, 0, 0, 0, 0, 1 };
    std::array<double, 3> dir = vectormath::cross(cos);

    std::array<double, 3> pos { 0, 1, 0 };

    const auto [b1, b2] = vectormath::splice(cos);

    auto posFin = vectormath::changeBasisInverse(b1, b2, dir, pos);

    for (std::size_t i = 0; i < 3; ++i)
        assert((pos[i] == posFin[i]));
}

int main(int argc, char* argv[])
{
    testLenght<double>();
    testLenght<float>();
    testArgMax();
    testChangeBasis();
    return EXIT_SUCCESS;
}
