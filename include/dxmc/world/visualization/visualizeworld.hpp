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
#include "dxmc/beams/beamtype.hpp"
#include "dxmc/floating.hpp"
#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/kdtreeintersectionresult.hpp"
#include "dxmc/world/visualization/vizualizationprops.hpp"
#include "dxmc/world/world.hpp"
#include "dxmc/world/worlditems/worlditembase.hpp"

#include <array>
#include <chrono>
#include <iterator>
#include <numbers>
#include <thread>
#include <unordered_map>
#include <vector>

#ifdef DXMCLIB_USE_LOADPNG
#include "lodepng/dxmclodepngwrapper.hpp"
#endif

namespace dxmc {

template <typename T>
    requires(std::same_as<T, double> || std::same_as<T, float> || std::same_as<T, std::uint8_t>)
struct VisualizationBuffer {
    std::vector<T> buffer;
    std::size_t width = 512;
    std::size_t height = 512;
    std::chrono::duration<double, std::milli> renderTime;

    VisualizationBuffer(std::size_t size)
    {
        width = size;
        height = size;
        buffer.resize(width * height * 4);
    }
    VisualizationBuffer(std::size_t w, std::size_t h)
    {
        width = w;
        height = h;
        buffer.resize(width * height * 4);
    }
    void resize(std::size_t size)
    {
        width = size;
        height = size;
        buffer.resize(width * height * 4);
    }
    void resize(std::size_t w, std::size_t h)
    {
        width = w;
        height = h;
        buffer.resize(width * height * 4);
    }
    void clear()
    {
        if constexpr (std::same_as<T, std::uint8_t>)
            std::fill(begin(), end(), 255);
        else
            std::fill(begin(), end(), T { 1 });
    }
    auto begin() { return buffer.begin(); }
    auto end() { return buffer.end(); }
};

template <typename U, typename T>
concept WorldType = requires(U world, Particle<T> p, KDTreeIntersectionResult<T, WorldItemBase<T>> res) {
    {
        world.intersect(p)
    } -> std::same_as<KDTreeIntersectionResult<T, WorldItemBase<T>>>;
};

template <Floating T>
class VisualizeWorld {
public:
    template <WorldType<T> W>
    VisualizeWorld(const W& world)
    {
        m_center = world.center();
        m_world_aabb = world.AABB();
        suggestFOV();
        updateColorsFromWorld(world);
    }

    template <WorldType<T> W>
    void updateFromWorld(const W& world)
    {
        m_center = world.center();
        m_world_aabb = world.AABB();
        suggestFOV();
        updateColorsFromWorld(world);
    }

    template <typename U>
        requires(std::same_as<U, T> || std::same_as<U, std::uint8_t>)
    static bool savePNG(const std::string& filename, const std::vector<U>& buffer, std::size_t width, std::size_t height)
    {
#ifdef DXMCLIB_USE_LOADPNG
        dxmclodepng::savePNG(filename, buffer, width, height);
        return true;
#else
        return false;
#endif
    }

    template <typename U>
        requires(std::same_as<U, T> || std::same_as<U, std::uint8_t>)
    static bool savePNG(const std::string& filename, const VisualizationBuffer<U>& buffer)
    {
        return savePNG(filename, buffer.buffer, buffer.width, buffer.height);
    }

    void setDistance(T dist)
    {
        m_camera_pos[0] = std::abs(dist);
    }
    void setPolarAngle(T polar)
    {
        m_camera_pos[1] = polar;
    }
    void setAzimuthalAngle(T azimuthal)
    {
        m_camera_pos[2] = -azimuthal;
    }
    void setPolarAngleDeg(T polar)
    {
        setPolarAngle(polar * DEG_TO_RAD<T>());
    }
    void setAzimuthalAngleDeg(T azimuthal)
    {
        setAzimuthalAngle(azimuthal * DEG_TO_RAD<T>());
    }

    void setCameraPosition(const std::array<T, 3>& pos)
    {
        const auto c = vectormath::subtract(m_center, pos);
        const auto l = vectormath::lenght(c);
        m_camera_pos[0] = l;
        const auto k = std::sqrt(c[0] * c[0] + c[1] * c[1]);
        if (c[2] > k)
            m_camera_pos[1] = std::atan2(k, c[2]);
        else
            m_camera_pos[1] = std::acos(c[2] / l);

        if (std::abs(c[0]) > T { 1e-6 })
            m_camera_pos[2] = std::atan2(c[1], c[0]);
        else if (std::abs(c[1]) > T { 1e-6 })
            m_camera_pos[2] = std::acos(c[0] / k);
        else
            m_camera_pos[2] = std::numbers::pi_v<T> / 2;
    }

    void suggestFOV(T zoom = 1)
    {
        const auto [p1, p2] = vectormath::splice(m_world_aabb);
        const auto plen = vectormath::lenght(vectormath::subtract(p1, p2));
        const auto clen = vectormath::lenght(vectormath::subtract(m_camera_pos, m_center));
        m_fov = std::atan(T { 0.5 } * plen / clen) / zoom;
    }

    void addLineProp(const std::array<T, 3>& start, const std::array<T, 3>& dir, T lenght = -1, T radii = 1)
    {
        m_lines.push_back({ start, dir, lenght, radii });
    }

    void clearLineProps()
    {
        m_lines.clear();
    }

    template <typename B>
        requires requires(B obj, std::array<T, 3> pos, std::array<std::array<T, 3>, 2> cos, std::array<T, 2> ang) {
            {
                obj.position()
            } -> std::convertible_to<std::array<T, 3>>;
            {
                obj.directionCosines()
            } -> std::convertible_to<std::array<std::array<T, 3>, 2>>;

            obj.collimationAngles();
        }

    void addLineProp(const B& obj, T lenght = -1, T radii = 1)
    {
        const auto start = obj.position();
        const auto cos = obj.directionCosines();
        const auto dir = vectormath::cross(cos[0], cos[1]);

        constexpr bool has_advanced_collimation = requires(B obj) {
            {
                obj.collimationAngles()
            } -> std::convertible_to<std::array<T, 4>>;
        };
        constexpr bool has_simple_collimation = requires(B obj) {
            {
                obj.collimationAngles()
            } -> std::convertible_to<std::array<T, 2>>;
        };

        std::array<T, 4> angs;
        if constexpr (has_advanced_collimation) {
            angs = obj.collimationAngles();
        } else if constexpr (has_simple_collimation) {
            auto sang = obj.collimationAngles();
            angs = { -sang[0], -sang[1], sang[0], sang[1] };
        } else {
            return;
        }

        auto f = [](const auto& cosines, const auto& dir, T angx, T angy) -> std::array<T, 3> {
            const auto sinx = std::sin(angx);
            const auto siny = std::sin(angy);
            const auto sinz = std::sqrt(1 - sinx * sinx - siny * siny);
            std::array pdir = {
                cosines[0][0] * sinx + cosines[1][0] * siny + dir[0] * sinz,
                cosines[0][1] * sinx + cosines[1][1] * siny + dir[1] * sinz,
                cosines[0][2] * sinx + cosines[1][2] * siny + dir[2] * sinz
            };
            return pdir;
        };

        addLineProp(start, f(cos, dir, angs[2], angs[3]), lenght, radii);
        addLineProp(start, f(cos, dir, angs[0], angs[3]), lenght, radii);
        addLineProp(start, f(cos, dir, angs[2], angs[1]), lenght, radii);
        addLineProp(start, f(cos, dir, angs[0], angs[1]), lenght, radii);
    }

    template <typename U = std::uint8_t>
        requires(std::same_as<U, T> || std::same_as<U, std::uint8_t>)
    static auto createBuffer(std::size_t width = 512, std::size_t height = 512)
    {
        VisualizationBuffer<U> buffer(width, height);
        return buffer;
    }

    template <WorldType<T> W, typename U>
        requires(std::same_as<U, T> || std::same_as<U, std::uint8_t>)
    void generate(W& world, VisualizationBuffer<U>& buffer) const
    {
        const auto t_start = std::chrono::high_resolution_clock::now();
        generate(world, buffer.buffer, buffer.width, buffer.height);
        const auto t_end = std::chrono::high_resolution_clock::now();
        buffer.renderTime = t_end - t_start;
    }

    template <WorldType<T> W, typename U>
        requires(std::same_as<U, T> || std::same_as<U, std::uint8_t>)
    void generate(W& world, std::vector<U>& buffer, int width = 512, int height = 512) const
    {
        if constexpr (std::is_same<U, std::uint8_t>::value)
            std::fill(buffer.begin(), buffer.end(), 255);
        else
            std::fill(buffer.begin(), buffer.end(), U { 1 });

        const auto xcos = vectormath::rotate({ T { 0 }, T { 1 }, T { 0 } }, { T { 0 }, T { 0 }, T { 1 } }, m_camera_pos[1]);
        const auto ycos = vectormath::rotate({ T { 0 }, T { 0 }, T { 1 } }, xcos, m_camera_pos[2] - std::numbers::pi_v<T> / 2);
        const auto dir = vectormath::cross(ycos, xcos);
        const std::array pos = {
            m_camera_pos[0] * std::cos(m_camera_pos[1]) * std::sin(m_camera_pos[2]) + m_center[0],
            m_camera_pos[0] * std::sin(m_camera_pos[1]) * std::sin(m_camera_pos[2]) + m_center[1],
            m_camera_pos[0] * std::cos(m_camera_pos[2]) + m_center[2]
        };

        const auto len = std::tan(m_fov);
        const auto xstart = vectormath::scale(xcos, -len);
        const auto ystart = vectormath::scale(ycos, -(len * height) / width);
        const auto step = 2 * len / (width - 1);

        auto worker = [&](std::size_t start, std::size_t stop) {
            Particle<T> p;
            p.pos = pos;

            for (std::size_t j = start; j < stop; ++j) {
                const auto y = j / width;
                const auto x = j - y * width;
                const auto ind = j * 4;

                const auto xvec = vectormath::add(xstart, vectormath::scale(xcos, x * step));
                const auto yvec = vectormath::add(ystart, vectormath::scale(ycos, y * step));
                p.dir = vectormath::normalized(vectormath::add(dir, xvec, yvec));

                T line_intersection = std::numeric_limits<T>::max();
                std::array<T, 3> line_normal = { 0, 0, 0 };
                for (const auto& line : m_lines) {
                    const auto opt = line.intersect(p);
                    if (opt) {
                        const auto& [l_cand, l_normal] = opt.value();
                        if (l_cand < line_intersection) {
                            line_intersection = l_cand;
                            line_normal = l_normal;
                        }
                    }
                }

                const auto res = world.intersectVisualization(p);

                if (res.valid()) {
                    if (res.intersection < line_intersection) {
                        const auto scaling = 1 + T { 0.5 } * vectormath::dot(p.dir, res.normal);

                        const auto color = colorOfItem<U>(res.item);
                        for (int i = 0; i < 3; ++i) {
                            if constexpr (std::is_same<U, std::uint8_t>::value) {
                                buffer[ind + i] = static_cast<U>(color[i] * scaling);
                            } else {
                                buffer[ind + i] = color[i] * scaling;
                            }
                        }
                    } else {
                        const auto scaling = std::abs(vectormath::dot(p.dir, line_normal));
                        const auto color = colorOfItem<U>(nullptr);
                        for (int i = 0; i < 3; ++i) {
                            if constexpr (std::is_same<U, std::uint8_t>::value) {
                                buffer[ind + i] = static_cast<U>(color[i] * scaling);
                            } else {
                                buffer[ind + i] = color[i] * scaling;
                            }
                        }
                    }
                } else {
                    if (line_intersection < std::numeric_limits<T>::max()) {
                        const auto scaling = std::abs(vectormath::dot(p.dir, line_normal));
                        const auto color = colorOfItem<U>(nullptr);
                        for (int i = 0; i < 3; ++i) {
                            if constexpr (std::is_same<U, std::uint8_t>::value) {
                                buffer[ind + i] = static_cast<U>(color[i] * scaling);
                            } else {
                                buffer[ind + i] = color[i] * scaling;
                            }
                        }
                    }
                }
            }
        };

        const auto nThreads = std::max(static_cast<std::uint32_t>(std::thread::hardware_concurrency()), std::uint32_t { 1 });
        const auto jobsize = static_cast<std::size_t>(width * height / nThreads);
        std::vector<std::jthread> threads;
        threads.reserve(nThreads - 1);
        std::size_t start = 0;
        auto stop = jobsize;
        for (auto i = nThreads - 1; i > 0; --i) {
            threads.emplace_back(std::jthread(worker, start, stop));
            start = stop;
            stop += jobsize;
        }

        worker(start, static_cast<std::size_t>(height * width));
        for (auto& t : threads) {
            t.join();
        }
    }

protected:
    static std::array<std::uint8_t, 3> HSVtoRGB(const std::uint8_t H, const std::uint8_t S, const std::uint8_t V = 255)
    {
        // H in [0, 255], S in [0, 255], V in [0, 255]
        const std::uint8_t region = H / 43;
        const std::uint8_t remainder = (H - (region * 43)) * 6;

        const std::uint8_t p = (V * (255 - S)) >> 8;
        const std::uint8_t q = (V * (255 - ((S * remainder) >> 8))) >> 8;
        const std::uint8_t t = (V * (255 - ((S * (255 - remainder)) >> 8))) >> 8;

        std::array<std::uint8_t, 3> rgb;
        if (region == 0)
            rgb = { V, t, p };
        else if (region == 1)
            rgb = { q, V, p };
        else if (region == 2)
            rgb = { p, V, t };
        else if (region == 3)
            rgb = { p, q, V };
        else if (region == 4)
            rgb = { t, p, V };
        else
            rgb = { V, p, q };
        return rgb;
    }

    static std::array<T, 3> HSVtoRGB(T h, T s = 1, T v = 1)
    {
        const T c = v * s;
        const T h_prime = (h * 3) / std::numbers::pi_v<T>;
        const T x = c * (1 - std::abs(std::fmod(h_prime, 2) - 1));
        const T m = v - c;

        std::array<T, 3> rgb;
        if (h_prime >= 0 && h_prime < 1) {
            rgb = { c, x, T { 0 } };
        } else if (h_prime >= 1 && h_prime < 2) {
            rgb = { x, c, T { 0 } };
        } else if (h_prime >= 2 && h_prime < 3) {
            rgb = { T { 0 }, c, x };
        } else if (h_prime >= 3 && h_prime < 4) {
            rgb = { T { 0 }, x, c };
        } else if (h_prime >= 4 && h_prime < 5) {
            rgb = { x, T { 0 }, c };
        } else {
            rgb = { c, T { 0 }, x };
        }
        return rgb;
    }

    template <WorldType<T> W>
    void updateColorsFromWorld(const W& world)
    {
        const std::vector<const WorldItemBase<T>*> items = world.getItemPointers();
        const auto N = items.size();
        m_colors_char.resize(N);
        m_colors_float.resize(N);
        for (std::size_t i = 0; i < N; ++i) {
            m_colors_char[i] = HSVtoRGB(static_cast<std::uint8_t>((i * 55) % 255), 255, 200);
            m_colors_float[i] = HSVtoRGB((2 * i * std::numbers::pi_v<T>) / N, T { 1 }, T { 0.8 });
            m_colorIndex[items[i]] = i;
        }
    }

    template <typename U = std::uint8_t>
        requires(std::same_as<U, T> || std::same_as<U, std::uint8_t>)
    std::array<U, 3> colorOfItem(const WorldItemBase<T>* item) const
    {
        if (auto search = m_colorIndex.find(item); search != m_colorIndex.end()) {
            const auto index = search->second;
            if constexpr (std::same_as<U, std::uint8_t>)
                return m_colors_char[index];
            else
                return m_colors_float[index];
        } else {
            if constexpr (std::same_as<U, std::uint8_t>)
                return std::array<uint8_t, 3> { 0, 0, 0 };
            else
                return std::array<T, 3> { 0, 0, 0 };
        }
    }

private:
    std::array<T, 6> m_world_aabb;
    std::array<T, 3> m_center = { 0, 0, 0 };
    std::array<T, 3> m_camera_pos = { 100, 1.4, 1 };
    T m_fov = -1;
    std::vector<std::array<T, 3>> m_colors_float;
    std::vector<std::array<std::uint8_t, 3>> m_colors_char;
    std::unordered_map<const WorldItemBase<T>*, std::size_t> m_colorIndex;
    std::vector<visualizationprops::Line<T>> m_lines;
};

}