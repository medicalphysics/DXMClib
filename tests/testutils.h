

#pragma once

#include "dxmc/floating.h"

constexpr double ERRF = 1e-6;

template <dxmc::Floating T>
bool isEqual(T f1, T f2)
{
    return std::abs(f1 - f2) < ERRF;
}

template <typename T>
bool isEqual(T f1, T f2)
{
    return f1 == f2;
}
