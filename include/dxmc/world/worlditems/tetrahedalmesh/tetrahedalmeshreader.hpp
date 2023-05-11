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

Copyright 2023 Erlend Andersen
*/

#pragma once

#include "dxmc/floating.hpp"

#include <fstream>
#include <ranges>
#include <string>
#include <string_view>
#include <vector>

namespace dxmc {
template <Floating T>
class TetrahedalmeshReader {
public:
    TetrahedalmeshReader() { }

    bool readVertices(const std::string& path, int nHeaderLines = 0)
    {
        // Open the file for reading
        std::ifstream file(path);
        if (!file.is_open()) {
            return false;
        }

        std::stringstream buffer;
        buffer << file.rdbuf();
        auto data = buffer.str();

        // auto test = std::views::split(data, '\n') | std::ranges::to<std::vector<std::string>>();
        auto test = std::views::split(data, '\n') | std::ranges::views::transform([](auto str) -> std::optional<T> {
            // std::string_view str(s.begin(), s.end());
            // if (T value; std::from_chars(str.begin(), str.end(), value).ec == std::errc {})
            //   return value;
            // else
            auto start = str.begin();
            auto end = str.end();
            std::string_view s { start, end };

            T value;
            std::from_chars(s.data(), s.data() + s.size(), value);
            auto test = *start;
            return std::nullopt;
        }) | std::ranges::to<std::vector<std::optional<T>>>();

        // Skip the header lines
        for (int i = 0; i < nHeaderLines; ++i) {
            std::string line;
            std::getline(file, line);
        }
        // Read the file and extract the floating point numbers using ranges
        // auto test  = std::ranges::istream_view<std::string>(file) | std::ranges::views::split(' ')
        //    | std::ranges::views::transform([](std::string s) { return std::stof(s); }) | std::ranges::to<std::vector<T>>();

        // constexpr char split = '\n';

        // auto test = std::ranges::istream_view<char>(file) | std::ranges::to<std::vector<char>>();

        // Close the file
        file.close();
        return true;
    }

private:
    std::vector<T> m_vertices;
};
}