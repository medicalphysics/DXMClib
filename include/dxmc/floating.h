

#pragma once
#include <type_traits>

namespace dxmc {
template <typename T>
concept Floating = std::is_floating_point<T>::value;
}