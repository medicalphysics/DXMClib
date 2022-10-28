/*This file is part of DXMClib.

DXMClib is free software : you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

DXMClib is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with DXMClib. If not, see < https://www.gnu.org/licenses/>.

Copyright 2022 Erlend Andersen
*/

#pragma once

#include <algorithm>
#include <concepts>
#include <iterator>
#include <vector>

template <typename T>
concept Number = std::is_integral<T>::value || std::is_floating_point<T>::value;

template <Number T>
void serialize(T in, std::vector<char>& buffer)
{
    auto dest = std::back_inserter(buffer);
    auto in_c = reinterpret_cast<char*>(&in);
    std::copy(in_c, in_c + sizeof(T), dest);
}

template <typename T>
void serialize(const std::vector<T>& data, std::vector<char>& buffer)
{
    if constexpr (std::same_as<T, char>) {
        std::copy(data.cbegin(), data.cend(), std::back_inserter(buffer));
    } else {
        std::uint64_t size = data.size() * sizeof(T);
        serialize(size, buffer);
        auto in_c = reinterpret_cast<const char*>(data.data());
        auto dest = std::back_inserter(buffer);
        std::copy(in_c, in_c + size, dest);
    }
}




template <Number T>
char* deserialize(T& value, char* begin)
{
    auto val_ptr = reinterpret_cast<T*>(begin);
    value = *val_ptr;
    return begin + sizeof(T);
}
template<typename T>
char* deserialize(std::vector<T>& val, char* begin, std::size_t size)
{
    auto n_elements = size / (sizeof(T));
    val.clear();
    auto start = reinterpret_cast<T*>(begin);
    std::copy(start, start + n_elements, std::back_inserter(val));
    return begin + n_elements * sizeof(T);
}

template<typename T>
char* deserialize(std::vector<T>& val, char* begin)
{
    std::uint64_t size;
    auto start = deserialize(size, begin);
    return deserialize(val, start, size);
}
