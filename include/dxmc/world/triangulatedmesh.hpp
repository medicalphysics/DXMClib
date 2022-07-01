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
class TriangulatedMesh {
public:
    TriangulatedMesh() { }
    TriangulatedMesh(const std::vector<std::array<T, 3>>& vertices, const std::vector<std::array<std::size_t, 3>>& faceIndices)
        : m_vertices(vertices)
        , m_faceIdx(faceIndices)
    {
        reduce();
    }

    const std::array<T, 3>& getVertice(std::size_t index) const { return m_vertices[index]; }
    const std::array<std::size_t, 3>& getFaceIndex(std::size_t index) const { return m_faceIdx[index]; }
    const std::size_t nVertices() const { return m_vertices.size(); }
    const std::size_t nFaces() const { return m_faceIdx.size(); }
    const std::vector<std::array<T, 3>>& getVertices() { return m_vertices; }
    const std::vector<std::array<std::size_t, 3>>& getFaceIndices() { return m_faceIdx; }



    template <int FORWARD>
    std::optional<T> intersect(const Particle<T>& particle) const
    {
        T t = std::numeric_limits<T>::max();

        for (std::size_t i = 0; i < m_nFaces; i++) {
            const auto vres = intersect<FORWARD>(particle, i);
            if (vres)
                t = std::min(t, *vres);
            // intersect = vres ? std::min(*vres, intersect) : intersect;
        }
        return t == std::numeric_limits<T>::max() ? std::nullopt : std::optional<T> { t };
    }
    template <int FORWARD>
    std::optional<T> intersect(const Particle<T>& particle, const std::size_t& facetIndex) const
    {
        constexpr T epsilon = std::numeric_limits<T>::epsilon();
        const auto& v1 = m_vertices[m_faceIdx[facetIndex][0]];
        const auto& v2 = m_vertices[m_faceIdx[facetIndex][1]];
        const auto& v3 = m_vertices[m_faceIdx[facetIndex][2]];

        const std::array<T, 3> v1v2 { v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2] };
        const std::array<T, 3> v1v3 { v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2] };
        const auto& pvec = vectormath::cross(particle.dir, v1v3);
        const auto det = vectormath::dot(v1v2, pvec);

        // Negative det gives intersection on backside of face
        if constexpr (FORWARD == 1) {
            if (det < epsilon)
                return std::nullopt;
        } else if constexpr (FORWARD == -1) {
            if (det > -epsilon)
                return std::nullopt;
        } else {
            if (std::abs(det) < epsilon)
                return std::nullopt;
        }
        const T invDet = T { 1 } / det;

        const std::array<T, 3> tvec { particle.pos[0] - v1[0], particle.pos[1] - v1[1], particle.pos[2] - v1[2] };
        const T u = vectormath::dot(tvec, pvec) * invDet;
        if (u < T { 0 } || u > T { 1 })
            return std::nullopt;

        const auto qvec = vectormath::cross(tvec, v1v2);
        const T v = vectormath::dot(particle.dir, qvec) * invDet;
        if (v < T { 0 } || u + v > T { 1 })
            return std::nullopt;

        return std::optional<T> { vectormath::dot(v1v3, qvec) * invDet };
    }

    bool setData(const std::vector<std::array<T, 3>>& vertices, const std::vector<std::array<std::size_t, 3>>& faceIndices)
    {
        m_vertices = vertices;
        m_faceIdx = faceIndices;
        reduce();
        return true;
    }

    bool setFlatData(const std::vector<T>& vertices, const std::vector<std::size_t>& faceIndices)
    {
        const bool valid_n_points = vertices.size() % 3 == 0 && faceIndices.size() % 3 == 0;
        const bool valid_indices = vertices.size() / 3 > std::reduce(std::execution::par_unseq, faceIndices.cbegin(), faceIndices.cend(), std::size_t { 0 }, [](const auto l, const auto r) { return std::max(l, r); });

        if (valid_n_points && valid_indices) {
            m_nVertices = vertices.size() / 3;
            m_nFaces = faceIndices.size() / 3;
            m_vertices.resize(m_nVertices);
            m_faceIdx.resize(m_nFaces);
            for (std::size_t i = 0; i < m_nVertices; ++i) {
                for (std::size_t j = 0; j < 3; ++j) {
                    const auto flat_ind = i * 3 + j;
                    m_vertices[i][j] = vertices[flat_ind];
                }
            }
            for (std::size_t i = 0; i < m_nFaces; ++i) {
                for (std::size_t j = 0; j < 3; ++j) {
                    const auto flat_ind = i * 3 + j;
                    m_faceIdx[i][j] = flat_ind;
                }
            }

            reduce();
            return true;
        }
        return false;
    }

protected:
    void reduce()
    {
        std::vector<std::array<T, 3>> r_vertices;
        std::vector<std::size_t> r_faceIdxFlat(m_faceIdx.size() * 3);

        for (std::size_t i = 0; i < m_faceIdx.size(); i++) {
            for (std::size_t j = 0; j < 3; j++) {
                r_faceIdxFlat[i * 3 + j] = m_faceIdx[i][j];
            }
        }

        auto equal = [](const std::array<T, 3>& l, const std::array<T, 3>& r) -> bool {
            constexpr T epsilon = std::numeric_limits<T>::epsilon();
            return std::abs(l[0] - r[0]) < epsilon && std::abs(l[1] - r[1]) < epsilon && std::abs(l[2] - r[2]) < epsilon;
        };
        for (std::size_t i = 0; i < m_vertices.size(); i++) {

            const auto& l = m_vertices[i];
            const auto handled = std::any_of(std::execution::par, r_vertices.cbegin(), r_vertices.cend(), [&](const auto& r) { return equal(r, l); });
            if (!handled) {
                r_vertices.push_back(l);
            }
            for (std::size_t j = 0; j < r_vertices.size(); j++) {
                const auto& r = r_vertices[j];
                if (equal(l, r)) {
                    std::replace(std::execution::par_unseq, r_faceIdxFlat.begin(), r_faceIdxFlat.end(), i, j);
                }
            }
        }
        for (std::size_t i = 0; i < m_faceIdx.size(); i++) {
            for (std::size_t j = 0; j < 3; j++) {
                m_faceIdx[i][j] = r_faceIdxFlat[i * 3 + j];
            }
        }
        m_vertices = r_vertices;
        m_nVertices = m_vertices.size();
    }

private:
    std::size_t m_nVertices = 0;
    std::size_t m_nFaces = 0;
    std::vector<std::array<T, 3>> m_vertices;
    std::vector<std::array<std::size_t, 3>> m_faceIdx;
};

template <Floating T>
class STLReader {
public:
    STLReader(const std::string& path)
        : m_filePath(path)
    {
    }
    STLReader() { }

    const std::string& message() const { return m_error; }

    TriangulatedMesh<T> operator()()
    {

        return operator()(m_filePath);
    }
    TriangulatedMesh<T> operator()(const std::string& path)
    {
        m_filePath = path;
        TriangulatedMesh<T> mesh;
        // First we chech for binary or ascii file
        std::ifstream f(path);
        if (f.is_open()) {
            std::string line;
            std::getline(f, line);
            const auto pos = line.find("solid", 0);
            f.close();
            if (pos != std::string::npos) {
                mesh = readSTLfileASCII(path);
            } else {
                mesh = readSTLfileBinary(path);
            }
        } else {
            m_error = "Could not open file: " + path;
            f.close();
        }
        return mesh;
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

    TriangulatedMesh<T> readSTLfileBinary(const std::string& path)
    {
        constexpr auto MIN_STL_SIZE = 80 + 4 + 50;

        TriangulatedMesh<T> mesh;
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
            return mesh;
        }

        const std::uint32_t n_triangles = readFromBytes<std::uint32_t>(&buffer[80]);
        const bool valid_file = MIN_STL_SIZE + 50 * (n_triangles - 1) == buffer.size();
        if (!valid_file) {
            m_error = "Error reading triangles";
            return mesh;
        } else {
            const std::vector<float> data = readTrianglesFromBytes(buffer);

            std::vector<T> vertices;
            std::vector<std::size_t> indices;
            std::vector<T> normals;
            std::size_t teller = 0;
            for (std::size_t i = 0; i < data.size(); i = i + 12) {
                for (std::size_t j = i; j < i + 3; j++) {
                    normals.push_back(data[j]);
                }
                for (std::size_t j = i + 3; j < i + 12; j++) {
                    vertices.push_back(data[j]);
                }
                for (std::size_t j = 0; j < 3; j++) {
                    indices.push_back(teller++);
                }
            }

            auto validmesh = mesh.setFlatData(vertices, indices);
            if (!validmesh) {
                m_error = "Mesh not valid for some reason";
            }
        }

        return mesh;
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

    TriangulatedMesh<T> readSTLfileASCII(const std::string& path)
    {
        TriangulatedMesh<T> mesh;

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

        std::vector<std::array<T, 3>> vertices;
        std::vector<std::size_t> faceIdxFlat;

        std::ifstream f(path, std::ios::in);
        if (f.is_open()) {
            std::size_t teller { 0 };
            for (std::string line; std::getline(f, line);) {
                auto res = processLine(line);
                if (res) {
                    vertices.push_back(*res);
                }
            }
            faceIdxFlat.resize(vertices.size());
            std::iota(faceIdxFlat.begin(), faceIdxFlat.end(), 0);
        } else {
            m_error = "Could not open file: " + path;
        }
        f.close();
        if (vertices.size() < 3) {
            m_error = "File do not contains any triangles";
            return mesh;
        }
        std::vector<std::array<std::size_t, 3>> faceIdx(faceIdxFlat.size() / 3);
        for (std::size_t i = 0; i < faceIdx.size(); i++) {
            for (std::size_t j = 0; j < 3; j++) {
                faceIdx[i][j] = faceIdxFlat[i * 3 + j];
            }
        }
        mesh.setData(vertices, faceIdx);

        return mesh;
    }

private:
    std::string m_error;
    std::string m_filePath;
};
}