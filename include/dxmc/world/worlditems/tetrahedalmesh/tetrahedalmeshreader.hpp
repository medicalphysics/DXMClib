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

#include <execution>
#include <fstream>
#include <limits>
#include <ranges>
#include <string>
#include <string_view>
#include <vector>

namespace dxmc {
template <Floating T>
class TetrahedalmeshReader {
public:
    TetrahedalmeshReader() { }
    bool readTetrahedalIndices(const std::string& path, int nHeaderLines = 1)
    {
        // reads a file formatted as <index n0 n1 n2 n3 matIdx000>, i.e "512 51 80 90 101"

        // Open the file for reading
        std::ifstream file(path);
        if (!file.is_open()) {
            return false;
        }

        std::stringstream buffer;
        buffer << file.rdbuf();
        const auto data = buffer.str();
        file.close();

        // finding line endings
        std::vector<std::pair<std::size_t, std::size_t>> lineIdx;
        std::size_t start = 0;
        for (std::size_t i = 0; i < data.size(); ++i) {
            if (data[i] == '\n') {
                lineIdx.push_back({ start, i });
                start = i;
            }
        }
        auto begin = lineIdx.cbegin() + nHeaderLines;
        auto end = lineIdx.cend();
        if (end - begin < 4)
            return false;

        std::vector<std::array<std::size_t, 6>> nodes(lineIdx.size());

        std::transform(std::execution::par_unseq, begin, end, nodes.begin(), [&data](const auto& idx) {
            const auto start = data.cbegin() + idx.first + 1;
            const auto stop = data.cbegin() + idx.second;
            const std::string_view line { start, stop };

            std::array<std::size_t, 6> v;

            int arrIdx = 0;

            auto word_start = line.data();
            auto line_end = line.data() + line.size();
            while (word_start != line_end) {
                if (*word_start == ' ') {
                    ++word_start;
                } else {
                    if (arrIdx < 6) {
                        auto [ptr, ec] = std::from_chars(word_start, line_end, v[arrIdx]);
                        if (ec == std::errc()) {
                            word_start = ptr;
                            ++arrIdx;
                        } else if (ec == std::errc::invalid_argument)
                            ++word_start;
                        else if (ec == std::errc::result_out_of_range)
                            word_start = ptr;
                    }
                }
            }
            return v;
        });

        // sort nodes
        std::sort(std::execution::par_unseq, nodes.begin(), nodes.end(), [](const auto& lh, const auto& rh) { return lh[0] < rh[0]; });

        // removing nan values and validating;
        m_tetrahedalIdx.clear();
        m_materialIdx.clear();

        m_tetrahedalIdx.reserve(nodes.size());
        m_materialIdx.reserve(nodes.size());

        for (std::size_t i = 0; i < nodes.size(); ++i) {
            if (i == nodes[i][0]) {
                std::array<std::size_t, 4> v;
                for (std::size_t j = 0; j < 4; ++j)
                    v[j] = nodes[i][j + 1];
                m_tetrahedalIdx.push_back(v);
                m_materialIdx.push_back(nodes[i][5]);
            }
        }

        return true;
    }

    bool readVertices(const std::string& path, int nHeaderLines = 1)
    {
        // reads a file formatted as <index v0 v1 v2>, i.e "512 0.2 0.4523 -0.974"

        // Open the file for reading
        std::ifstream file(path);
        if (!file.is_open()) {
            return false;
        }

        std::stringstream buffer;
        buffer << file.rdbuf();
        const auto data = buffer.str();
        file.close();

        // finding line endings
        std::vector<std::pair<std::size_t, std::size_t>> lineIdx;
        std::size_t start = 0;
        for (std::size_t i = 0; i < data.size(); ++i) {
            if (data[i] == '\n') {
                lineIdx.push_back({ start, i });
                start = i;
            }
        }
        auto begin = lineIdx.cbegin() + nHeaderLines;
        auto end = lineIdx.cend();
        if (end - begin < 4)
            return false;

        std::vector<std::pair<std::size_t, std::array<T, 3>>> vertices(lineIdx.size());

        std::transform(std::execution::par_unseq, begin, end, vertices.begin(), [&data](const auto& idx) {
            const auto start = data.cbegin() + idx.first + 1;
            const auto stop = data.cbegin() + idx.second;
            const std::string_view line { start, stop };

            std::array<T, 3> v = {
                std::numeric_limits<T>::quiet_NaN(),
                std::numeric_limits<T>::quiet_NaN(),
                std::numeric_limits<T>::quiet_NaN()
            };
            std::size_t index = std::numeric_limits<std::size_t>::max();
            int arrIdx = -1;

            auto word_start = line.data();
            auto line_end = line.data() + line.size();
            while (word_start != line_end) {
                if (*word_start == ' ') {
                    ++word_start;
                } else {
                    if (arrIdx < 0) {
                        auto [ptr, ec] = std::from_chars(word_start, line_end, index);
                        if (ec == std::errc()) {
                            word_start = ptr;
                            ++arrIdx;
                        } else if (ec == std::errc::invalid_argument)
                            ++word_start;
                        else if (ec == std::errc::result_out_of_range)
                            word_start = ptr;
                    } else if (arrIdx < 3) {
                        auto [ptr, ec] = std::from_chars(word_start, line_end, v[arrIdx]);
                        if (ec == std::errc()) {
                            word_start = ptr;
                            ++arrIdx;
                        } else if (ec == std::errc::invalid_argument)
                            ++word_start;
                        else if (ec == std::errc::result_out_of_range)
                            word_start = ptr;
                    }
                }
            }
            return std::make_pair(index, v);
        });

        // sort vertices
        std::sort(std::execution::par_unseq, vertices.begin(), vertices.end(), [](const auto& lh, const auto& rh) { return lh.first < rh.first; });

        // removing nan values and validating;
        m_vertices.clear();
        m_vertices.reserve(vertices.size());
        for (std::size_t i = 0; i < vertices.size(); ++i) {
            const auto& [ind, v] = vertices[i];
            if (v[0] == v[0] && v[1] == v[1] && v[2] == v[2] && ind == i) // test for NaNs
                m_vertices.push_back(v);
        }

        return true;
    }

private:
    std::vector<std::array<T, 3>> m_vertices;
    std::vector<std::array<std::size_t, 4>> m_tetrahedalIdx;
    std::vector<std::size_t> m_materialIdx;
    std::size_t m_maxNodeIndex = 0;
};
}