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
#include "dxmc/world/worlditems/tetrahedalmesh/tetrahedron.hpp"

#include <execution>
#include <fstream>
#include <limits>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

namespace dxmc {
template <Floating T>
class TetrahedalmeshReader {
public:
    TetrahedalmeshReader() { }

    std::vector<Tetrahedron<T>> readICRP145Phantom(const std::string& nodeFile, const std::string& elementFile)
    {
        auto vertices = readVertices(nodeFile, 1, 80);
        auto nodes = readTetrahedalIndices(elementFile, 1, 80);

        const auto max_ind = std::transform_reduce(
            std::execution::par_unseq, nodes.cbegin(), nodes.cend(), std::size_t { 0 }, [](const auto lh, const auto rh) { return std::max(rh, lh); }, [](const auto& node) {
                const auto& v = std::get<1>(node);
                return std::max(v[0], std::max(v[1], v[2])); });

        std::vector<Tetrahedron<T>> tets;
        if (max_ind < vertices.size()) {
            tets.resize(nodes.size());
            std::transform(std::execution::par_unseq, nodes.cbegin(), nodes.cend(), tets.begin(), [&vertices](const auto& n) {
                const auto& vIdx = std::get<1>(n);
                const auto collection = static_cast<std::uint8_t>(std::get<2>(n) / 1000);

                const auto& v0 = vertices[vIdx[0]].second;
                const auto& v1 = vertices[vIdx[1]].second;
                const auto& v2 = vertices[vIdx[2]].second;
                const auto& v3 = vertices[vIdx[3]].second;

                return Tetrahedron<T> { v0, v1, v2, v3, collection };
            });
        }

        return tets;
    }

    static std::vector<std::tuple<std::size_t, std::array<std::size_t, 4>, std::size_t>> readTetrahedalIndices(const std::string& path, int nHeaderLines = 1, std::size_t colLenght = 80)
    {
        // reads a file formatted as <index n0 n1 n2 n3 matIdx000>, i.e "512 51 80 90 101"

        std::vector<std::tuple<std::size_t, std::array<std::size_t, 4>, std::size_t>> nodes;

        const auto data = readBufferFromFile(path);
        if (data.size() == 0)
            return nodes;

        // finding line endings
        std::vector<std::pair<std::size_t, std::size_t>> lineIdx;
        lineIdx.reserve(data.size() / colLenght);
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
            return nodes;

        nodes.resize(end - begin);

        std::transform(std::execution::par_unseq, begin, end, nodes.begin(), [&data](const auto& idx) {
            const auto start = data.cbegin() + idx.first + 1;
            const auto stop = data.cbegin() + idx.second;
            const std::string_view line { start, stop };

            std::size_t index = std::numeric_limits<std::size_t>::max();

            std::array<std::size_t, 4> v {
                std::numeric_limits<std::size_t>::max(),
                std::numeric_limits<std::size_t>::max(),
                std::numeric_limits<std::size_t>::max(),
                std::numeric_limits<std::size_t>::max()
            };
            std::size_t collection = std::numeric_limits<std::size_t>::max();

            int arrIdx = -1;

            auto word_start = line.data();
            auto line_end = line.data() + line.size();
            while (word_start != line_end) {
                if (*word_start == ' ') {
                    ++word_start;
                } else if (*word_start == '#') {
                    // we exits if we find #
                    word_start = line_end;
                } else {
                    if (arrIdx < 5) {
                        std::size_t val;
                        auto [ptr, ec] = std::from_chars(word_start, line_end, val);
                        if (ec == std::errc()) {
                            word_start = ptr;
                            if (arrIdx < 0) {
                                index = val;
                            } else if (arrIdx < 4) {
                                v[arrIdx] = val;
                            } else {
                                collection = val;
                            }
                            ++arrIdx;
                        } else if (ec == std::errc::invalid_argument)
                            ++word_start;
                        else if (ec == std::errc::result_out_of_range)
                            word_start = ptr;
                    }
                }
            }
            return std::make_tuple(index, v, collection);
        });

        // sort nodes
        std::sort(std::execution::par_unseq, nodes.begin(), nodes.end(), [](const auto& lh, const auto& rh) { return std::get<0>(lh) < std::get<0>(rh); });

        // removing nan values and validating;
        auto delete_from = nodes.cbegin();
        std::size_t index = 0;
        while (delete_from != nodes.cend() && index == std::get<0>(*delete_from)) {
            ++index;
            ++delete_from;
        }

        if (delete_from != nodes.cend())
            nodes.erase(delete_from, nodes.cend());

        return nodes;
    }

    static std::vector<std::pair<std::size_t, std::array<T, 3>>> readVertices(const std::string& path, int nHeaderLines = 1, std::size_t colLenght = 80)
    {
        // reads a file formatted as <index v0 v1 v2>, i.e "512 0.2 0.4523 -0.974"
        // # is treated as comment start

        std::vector<std::pair<std::size_t, std::array<T, 3>>> vertices;

        const auto data = readBufferFromFile(path);
        if (data.size() == 0)
            return vertices;

        // finding line endings
        std::vector<std::pair<std::size_t, std::size_t>> lineIdx;
        lineIdx.reserve(data.size() / colLenght);
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
            return vertices;

        vertices.resize(end - begin);

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
            const auto line_end = line.data() + line.size();
            while (word_start != line_end) {
                if (*word_start == ' ') {
                    ++word_start;
                } else if (*word_start == '#') {
                    // we exits if we find #
                    word_start = line_end;
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
        auto delete_from = vertices.cbegin();
        std::size_t index = 0;
        while (delete_from != vertices.cend() && index == delete_from->first) {
            ++index;
            ++delete_from;
        }

        if (delete_from != vertices.cend())
            vertices.erase(delete_from, vertices.cend());

        return vertices;
    }

    bool validateIndices() const
    {
        std::vector<std::uint8_t> testIdx(m_vertices.size(), 0);
        const auto maxInd = m_vertices.size();
        bool valid = true;
        for (const auto& v : m_tetrahedalIdx) {
            for (const auto& t : v) {
                if (t < maxInd)
                    testIdx[t] = 1;
                else
                    valid = false;
            }
        }
        valid = valid && std::all_of(std::execution::par_unseq, testIdx.cbegin(), testIdx.cend(), [](const auto i) { return i > 0; });
        return valid;
    }

protected:
    static std::string readBufferFromFile(const std::string& path)
    {
        std::string buffer_str;

        // Open the file for reading
        std::ifstream file(path);
        if (file.is_open()) {
            std::stringstream buffer;
            buffer << file.rdbuf();
            buffer_str = buffer.str();
            file.close();
        }
        return buffer_str;
    }

private:
    std::vector<std::array<T, 3>> m_vertices;
    std::vector<std::array<std::size_t, 4>> m_tetrahedalIdx;
    std::vector<std::size_t> m_materialIdx;
    std::size_t m_maxNodeIndex = 0;
};
}