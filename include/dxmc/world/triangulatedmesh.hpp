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

#include "dxmc/floating.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/baseobject.hpp"
#include "dxmc/world/triangle.hpp"

#include <algorithm>
#include <array>
#include <execution>
#include <fstream>
#include <limits>
#include <numeric>
#include <optional>
#include <string>
#include <vector>

namespace dxmc {

template <Floating T>
class STLReader {
public:
    STLReader(const std::string& path)
        : m_filePath(path)
    {
    }
    STLReader() { }

    const std::string& message() const { return m_error; }

    std::vector<Triangle<T>> operator()()
    {
        return operator()(m_filePath);
    }
    std::vector<Triangle<T>> operator()(const std::string& path)
    {
        m_filePath = path;
        std::vector<Triangle<T>> data;
        // First we chech for binary or ascii file
        std::ifstream f(path);
        if (f.is_open()) {
            std::string line;
            std::getline(f, line);
            const auto pos = line.find("solid", 0);
            f.close();
            if (pos != std::string::npos) {
                data = readSTLfileASCII(path);
            } else {
                data = readSTLfileBinary(path);
            }
        } else {
            m_error = "Could not open file: " + path;
            f.close();
        }
        return data;
    }

protected:
    template <typename S>
    static S readFromBytes(const std::uint8_t* bytes)
    {
        S value;
        memcpy(&value, bytes, sizeof(S));
        return value;
    }

    std::vector<float> readTrianglesFromBytes(const std::vector<std::uint8_t>& buffer, const std::size_t header = 84, const std::size_t offset = 50) const
    {
        const auto n_elements = ((buffer.size() - 84) / 50);
        std::vector<float> vec(n_elements * 12);

        std::size_t vec_pos = 0;
        for (std::size_t i = 0; i < n_elements; i++) {
            const auto b = &buffer[84 + 50 * i];
            memcpy(&vec[vec_pos], b, sizeof(float) * 12);
            vec_pos += 12;
        }
        return vec;
    }

    std::vector<Triangle<T>> readSTLfileBinary(const std::string& path)
    {
        constexpr auto MIN_STL_SIZE = 80 + 4 + 50;

        std::vector<Triangle<T>> data;
        std::vector<std::uint8_t> buffer;

        std::ifstream f(path, std::ios::binary | std::ios::in | std::ios::ate);
        if (f.is_open()) {
            const auto filesize = f.tellg();
            f.seekg(0, std::ios::beg);
            buffer.resize(filesize);
            f.read((char*)&buffer[0], filesize);
        } else {
            m_error = "Could not open file: " + path;
        }
        f.close();

        if (buffer.size() < MIN_STL_SIZE) {
            m_error = "File do not contains any triangles";
            return data;
        }

        const std::uint32_t n_triangles = readFromBytes<std::uint32_t>(&buffer[80]);
        const bool valid_file = MIN_STL_SIZE + 50 * (n_triangles - 1) == buffer.size();
        if (!valid_file) {
            m_error = "Error reading triangles";
            return data;
        } else {
            const std::vector<float> data_buffer = readTrianglesFromBytes(buffer);

            std::vector<T> vertices;
            std::vector<T> normals;
            vertices.reserve(data_buffer.size() * 9);
            normals.reserve(data_buffer.size() * 3);
            for (std::size_t i = 0; i < data_buffer.size(); i = i + 12) {
                for (std::size_t j = i; j < i + 3; j++) {
                    normals.push_back(data_buffer[j]);
                }
                for (std::size_t j = i + 3; j < i + 12; j++) {
                    vertices.push_back(data_buffer[j]);
                }
            }

            const bool valid_n_points = vertices.size() % 9 == 0;

            if (valid_n_points) {
                const auto n_triangles = vertices.size() / 9;
                data.clear();
                data.reserve(n_triangles);
                for (std::size_t i = 0; i < vertices.size(); i += 9) {
                    data.emplace_back(&vertices[i]);
                }
            }
        }

        return data;
    }

    static std::vector<std::string> stringSplit(const std::string& text, char sep)
    {
        // this function splits a string into a vector based on a sep character
        // it will skip empty tokens
        std::vector<std::string> tokens;
        std::size_t start = 0, end = 0;
        while ((end = text.find(sep, start)) != std::string::npos) {
            if (end != start) {
                tokens.push_back(text.substr(start, end - start));
            }
            start = end + 1;
        }
        if (end != start) {
            tokens.push_back(text.substr(start));
        }
        return tokens;
    }

    std::vector<Triangle<T>> readSTLfileASCII(const std::string& path)
    {

        auto processLine = [](const std::string& line) -> std::optional<std::array<T, 3>> {
            auto words = stringSplit(line, ' ');
            if (words.size() < 4)
                return std::nullopt;
            const std::string cmp("vertex");
            std::transform(words[0].begin(), words[0].end(), words[0].begin(), [](unsigned char u) { return std::tolower(u); });
            if (cmp.compare(words[0]) == 0) {

                try {
                    if constexpr (sizeof(T) == 4) {
                        std::array<T, 3> res {
                            std::stof(words[1]),
                            std::stof(words[2]),
                            std::stof(words[3])
                        };
                        return std::optional(res);
                    } else {
                        std::array<T, 3> res {
                            std::stod(words[1]),
                            std::stod(words[2]),
                            std::stod(words[3])
                        };
                        return std::optional(res);
                    }

                } catch (std::invalid_argument) {
                    return std::nullopt;
                }
            }
            return std::nullopt;
        };

        std::vector<T> vertices;

        std::ifstream f(path, std::ios::in);
        if (f.is_open()) {
            std::size_t teller { 0 };
            for (std::string line; std::getline(f, line);) {
                auto res = processLine(line);
                if (res) {
                    for (auto val : *res)
                        vertices.push_back(val);
                }
            }

        } else {
            m_error = "Could not open file: " + path;
        }
        f.close();

        std::vector<Triangle<T>> mesh;
        if (vertices.size() < 9) {
            m_error = "File do not contains any triangles";
            return mesh;
        } else {
            const auto n_triangles = vertices.size() / 9;
            mesh.clear();
            mesh.reserve(n_triangles);
            for (std::size_t i = 0; i < vertices.size(); i += 9) {
                mesh.emplace_back(&vertices[i]);
            }
        }

        return mesh;
    }

private:
    std::string m_error;
    std::string m_filePath;
};

template <Floating T>
class TriangulatedMesh final : public BaseObject<T> {
public:
    TriangulatedMesh()
        : BaseObject<T>()
    {
    }
    TriangulatedMesh(const std::vector<Triangle<T>>& triangles, const std::size_t max_tree_dept = 8)
        : BaseObject<T>()
    {
        setData(triangles, max_tree_dept);
    }
    TriangulatedMesh(const std::string& path, const std::size_t max_tree_dept = 8)
        : BaseObject<T>()
    {
        STLReader<T> reader;
        auto triangles = reader(path);
        setData(triangles, max_tree_dept);
    }
    const KDTree<Triangle<T>>& kdtree() const { return m_kdtree; }
    std::vector<Triangle<T>> getTriangles() const
    {
        return m_kdtree.items();
    }
    void translate(const std::array<T, 3>& dist) override
    {
        m_kdtree.translate(dist);
        this->m_aabb = m_kdtree.AABB();
    }

    void setData(std::vector<Triangle<T>> triangles, const std::size_t max_tree_dept)
    {
        m_kdtree = KDTree(triangles, max_tree_dept);
        this->m_aabb = m_kdtree.AABB();
    }
    std::array<T, 3> center() const override
    {
        std::array<T, 3> center { 0, 0, 0 };
        const auto triangles = getTriangles();
        std::for_each(triangles.cbegin(), triangles.cend(), [&](const auto& tri) {
            const auto c = tri.center();
            for (std::size_t i = 0; i < 3; ++i) {
                center[i] += c[i];
            }
        });
        for (std::size_t i = 0; i < 3; ++i) {
            center[i] /= triangles.size();
        }
        return center;
    }
    std::optional<T> intersect(Particle<T>& p) const override
    {
        return m_kdtree.intersect(p, this->m_aabb);
    }
    std::optional<T> intersect(Particle<T>& p, const std::array<T, 2>& tbox) const override
    {
        return m_kdtree.intersect(p, this->m_aabb, tbox);
    }
    T transport(Particle<T>& p, RandomState& state) const override
    {
        return 0;
    }

private:
    KDTree<Triangle<T>> m_kdtree;
};
}