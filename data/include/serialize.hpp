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
void serialize(T* in, std::vector<char>& buffer)
{
    auto dest = std::back_inserter(buffer);
    auto in_c = reinterpret_cast<const char*>(&in);
    std::copy(in_c, in_c + sizeof(T), dest);
}

template <typename T>
void serialize(const std::vector<T>& data, std::vector<char>& buffer)
{
    std::uint64_t size = data.size() * sizeof(T);
    serialize(&size, buffer);

    auto in_c = reinterpret_cast<char const*>(data.data());
    auto dest = std::back_inserter(buffer);
    std::copy(in_c, in_c + size, dest);
}
