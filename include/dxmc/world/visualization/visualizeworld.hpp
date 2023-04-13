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
#include "dxmc/particle.hpp"
#include "dxmc/vectormath.hpp"
#include "dxmc/world/kdtreeintersectionresult.hpp"
#include "dxmc/world/world.hpp"
#include "dxmc/world/worlditems/worlditembase.hpp"

#include <array>
#include <iterator>
#include <numbers>
#include <vector>

namespace dxmc {
template <typename U, typename T>
concept WorldType = requires(U world, Particle<T> p, KDTreeIntersectionResult<T, WorldItemBase<T>> res) {
                        Floating<T>;
                        {
                            world.intersect(p)
                            } -> std::same_as<KDTreeIntersectionResult<T, WorldItemBase<T>>>;
                    };
template <Floating T>
class VisualizeWorld {
public:
    template <WorldType<T> W>
    VisualizeWorld(W& world)
    {
        m_center = world.center();
        if (m_fov < 0)
            suggestFOV(world.AABB());
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
        m_camera_pos[2] = azimuthal;
    }

    // template <WorldType<T> W, typename U>
    //   requires(std::same_as<U, T> || std::same_as<U, std::uint8_t>)
    template <WorldType<T> W>
    void generate(W& world, std::vector<T>& buffer, int width = 512, int height = 512)
    {
        using U = T;
        std::vector<std::pair<WorldItemBase<T>*, std::array<U, 3>>> colors;
        auto items = world.getItemPointers();
        colors.reserve(items.size());
        for (std::size_t i = 0; i < items.size(); ++i) {
            if constexpr (std::is_same<U, std::uint8_t>::value) {
                const std::uint8_t c = static_cast<std::uint8_t>(i + 55);
                colors.emplace_back(std::make_pair(items[i], HSVtoRGB(c, 255, 200)));
            } else {
                const auto c = (i * 2 * std::numbers::pi_v<T>) / items.size();
                colors.emplace_back(std::make_pair(items[i], HSVtoRGB(c, T { 1 }, T { 0.85 })));
            }
        }

        const auto step = std::max(m_fov / width, m_fov / height);
        Particle<T> p;
        p.pos = {
            m_camera_pos[0] * std::cos(m_camera_pos[1]) * std::sin(m_camera_pos[2]) + m_center[0],
            m_camera_pos[0] * std::sin(m_camera_pos[1]) * std::sin(m_camera_pos[2]) + m_center[1],
            m_camera_pos[0] * std::cos(m_camera_pos[2]) + m_center[2]
        };
        auto dir = vectormath::subtract(m_center, p.pos);
        vectormath::normalize(dir);
        const auto xcos = vectormath::rotate({ T { 0 }, T { 1 }, T { 0 } }, { T { 0 }, T { 0 }, T { 1 } }, m_camera_pos[1]);
        const auto ycos = vectormath::rotate({ T { 0 }, T { 0 }, T { 1 } }, xcos, m_camera_pos[2] - std::numbers::pi_v<T> / 2);
        const auto pdir = vectormath::cross(ycos, xcos);

        for (int y = 0; y < height; ++y) {
            const auto y_flat = y * width * 4;
            const auto yang = (-height / 2 + y) * step;
            const auto dy = step * y - height * step / 2;
            for (int x = 0; x < width; ++x) {
                const auto ind = y_flat + x * 4;
                const auto xang = (-width / 2 + x) * step;
                const auto dx = step * x - width * step / 2;
                p.dir = vectormath::rotate(vectormath::rotate(dir, ycos, xang), xcos, yang);

                const auto res = world.intersectVisualization(p);
                if (res.valid()) {
                    const auto scaling = vectormath::dot(p.dir, res.normal);
                    for (const auto& cmap : colors) {
                        if (cmap.first == res.item) {
                            for (int i = 0; i < 3; ++i)
                                buffer[ind + i] = cmap.second[i] * scaling;
                        }
                    }
                } else {
                    for (int i = 0; i < 3; ++i) {
                        if constexpr (std::is_same<U, std::uint8_t>::value) {
                            buffer[ind + i] = 255;
                        } else {
                            buffer[ind + i] = T { 1 };
                        }
                    }
                }
            }
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
        const auto c = v * s;
        const auto h_prime = (h * 3) / std::numbers::pi_v<T>;
        const auto x = c * (1 - std::abs(std::fmod(h_prime, 2) - 1));
        const auto m = v - c;

        std::array<T, 3> rgb;
        if (h_prime >= 0 && h_prime < 1) {
            rgb = { c, x, 0 };
        } else if (h_prime >= 1 && h_prime < 2) {
            rgb = { x, c, 0 };
        } else if (h_prime >= 2 && h_prime < 3) {
            rgb = { 0, c, x };
        } else if (h_prime >= 3 && h_prime < 4) {
            rgb = { 0, x, c };
        } else if (h_prime >= 4 && h_prime < 5) {
            rgb = { x, 0, c };
        } else {
            rgb = { c, 0, x };
        }
        return rgb;
    }

    /* static std::array<T, 3> HSVtoRGB(const T H, const T S = T { 1 }, const T V = .8)
    {
        // H in [0, 2pi], S in [0, 1], V in [0, 1]
        const auto C = V * S;

        constexpr auto sep = 3 / std::numbers::phi_v<T>;
        const auto Hb = std::clamp(H, T { 0 }, 2 * std::numbers::phi_v<T>) * sep;
        const auto r = std::fmod(Hb, T { 2 });
        const auto X = C * (1 - std::abs(r - 1));

        std::array<T, 3> rgb;
        if (Hb < 1) {
            rgb = { C, X, 0 };
        } else if (Hb < 2) {
            rgb = { X, C, 0 };
        } else if (Hb < 3) {
            rgb = { 0, C, X };
        } else if (Hb < 4) {
            rgb = { 0, X, C };
        } else if (Hb < 5) {
            rgb = { X, 0, C };
        } else {
            rgb = { C, 0, X };
        }
        return rgb;
    }*/

    void suggestFOV(const std::array<T, 6>& aabb)
    {
        const auto [p1, p2] = vectormath::splice(aabb);
        const auto len = vectormath::lenght(vectormath::subtract(p1, p2));

        m_fov = std::atan(T { 0.5 } * len / m_camera_pos[0]);
    }

private:
    std::array<T, 3> m_center = { 0, 0, 0 };
    std::array<T, 3> m_camera_pos = { 100, 1.4, 1 };
    std::array<T, 3> m_camera_xcosine = { 1, 0, 0 };
    std::array<T, 3> m_camera_ycosine = { 0, 1, 0 };
    T m_fov = -1;
};

}