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
    template <WorldType W, U T>
        requires std::is_floating_point<U>::value || std::same_as<U, std::uint8_t>::value
    void generate(W& world, std::vector<U>& buffer, int width = 512, int height = 512)
    {
        if (m_fov < 0)
            suggestFOV(world.AABB());
    }

protected:
    void suggestFOV(const std::array<T, 6>& aabb)
    {
        const auto [p1, p2] = vectormath::splice(aabb);
        auto d1 = vectormath::subtract(p1, m_camera_pos);
        vectormath::normalize(d1);
        auto d2 = vectormath::subtract(p2, m_camera_pos);
        vectormath::normalize(d2);

        const auto dir = vectormath::cross(m_camera_xcosine, m_camera_ycosine);

        const auto x = std::max(std::abs(vectormath::dot(dir, d1)),std::abs(vectormath::dot(dir, d2));
        
        m_fov = std::asin(std::max(x, y)) * 2;
    }

private:
    std::array<T, 3> m_camera_pos = { 0, 0, -1000 };
    std::array<T, 3> m_camera_xcosine = { 1, 0, 0 };
    std::array<T, 3> m_camera_ycosine = { 0, 1, 0 };
    T m_fov = -1;
};
namespace visualization {

    template <typename U, typename T>
    concept WorldType = requires(U world, Particle<T> p, KDTreeIntersectionResult<T, WorldItemBase<T>> res) {
                            Floating<T>;
                            {
                                world.intersect(p)
                                } -> std::same_as<KDTreeIntersectionResult<T, WorldItemBase<T>>>;
                        };

    template <Floating T>
    std::array<T, 3> projectToPlane(const std::array<T, 3>& planePos, const std::array<T, 3>& planeNormal, const std::array<T, 3>& point)
    {

        const auto np = vectormath::dot(planeNormal, point);
        const auto no = vectormath::dot(planeNormal, planePos);
        const auto d = vectormath::scale(np - no, planeNormal);
        return vectormath::subtract(planePos, d);
    }

    template <Floating T>
    std::pair<T, T> fieldOfView(const std::array<T, 6>& aabb, const std::array<T, 3>& cameraPos, const std::array<T, 3>& cameraCosinex, const std::array<T, 3>& cameraCosiney)
    {
        const auto dir = vectormath::cross(cameraCosinex, cameraCosiney);
        const std::array<T, 3> c = {
            (aabb[0] + aabb[3]) / 2,
            (aabb[1] + aabb[4]) / 2,
            (aabb[2] + aabb[5]) / 2
        };
        const auto [c0, c2] = vectormath::splice(aabb);
        const auto p0 = std::max(std::abs(vectormath::lenght(vectormath::subtract(c0, c))),
            std::abs(vectormath::lenght(vectormath::subtract(c2, c))));
        const auto p2 = std::max(std::abs(vectormath::lenght(vectormath::subtract(c0, c))),
            std::abs(vectormath::lenght(vectormath::subtract(c2, c))));

        const auto dist = vectormath::lenght(vectormath::subtract(c, cameraPos));

        const auto fx = std::tan(p0 / dist);
        const auto fy = std::tan(p2 / dist);

        return std::make_pair(fx, fy);
    }

    template <Floating T, WorldType<T> W>
    std::vector<T> rayTraceGeometryDistance(W& world, const std::array<T, 3>& cameraPos, const std::array<T, 3>& cameraCosinex, const std::array<T, 3>& cameraCosiney, std::size_t resolution = 128)
    {
        const auto aabb = world.AABB();
        const std::array<T, 3> c = {
            (aabb[0] + aabb[3]) / 2,
            (aabb[1] + aabb[4]) / 2,
            (aabb[2] + aabb[5]) / 2
        };

        const auto [fov_x, fov_y] = fieldOfView(world.AABB(), cameraPos, cameraCosinex, cameraCosiney);
        const auto fov = std::max(fov_x, fov_y);
        const auto pc = vectormath::subtract(c, cameraPos);

        const auto dd = std::tan(fov / resolution / vectormath::lenght(pc));
        const auto dd_start = -dd * resolution / 2;

        const auto dir = vectormath::cross(cameraCosinex, cameraCosiney);

        std::vector<T> data(resolution * resolution, T { 0 });

        for (std::size_t y = 0; y < resolution; ++y) {
            const auto y_flat = y * resolution;
            for (std::size_t x = 0; x < resolution; ++x) {
                const auto flat = y_flat + x;
                Particle<T> p = { .pos = cameraPos, .dir = dir };
                p.dir = vectormath::rotate(p.dir, cameraCosiney, dd_start + dd * x);
                p.dir = vectormath::rotate(p.dir, cameraCosinex, dd_start + dd * y);
                const auto r = world.intersect(p);
                if (r.valid())
                    data[flat] = r.intersection;
            }
        }

        return data;
    }

    std::array<std::uint8_t, 3> HSVtoRGB(const std::uint8_t H, const std::uint8_t S, const std::uint8_t V = 255)
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

    template <Floating T>
    std::array<T, 3> HSVtoRGB(const T H, const T S = T { 1 }, const T V = .8)
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
    }

    template <Floating T, WorldType<T> W>
    std::vector<T> rayTraceGeometry(W& world, const std::array<T, 3>& cameraPos, std::size_t resolution = 128, const T zoom = 1)
    {

        auto dir = vectormath::subtract(world.center(), cameraPos);
        vectormath::normalize(dir);
        const auto dInd = vectormath::argmax3(dir);
        std::array<T, 3> x, y;
        if (dInd == 0) {
            x = { 0, 1, 0 };
            y = { 0, 0, 1 };
        } else if (dInd == 1) {
            x = { 1, 0, 0 };
            y = { 0, 0, 1 };
        } else {
            x = { 1, 0, 0 };
            y = { 0, 1, 0 };
        }

        auto cameraCosinex = vectormath::cross(dir, y);
        vectormath::normalize(cameraCosinex);
        auto cameraCosiney = vectormath::cross(dir, x);
        vectormath::normalize(cameraCosiney);

        const auto d2 = vectormath::cross(cameraCosinex, cameraCosiney);

        const auto [fov_x, fov_y] = fieldOfView(world.AABB(), cameraPos, cameraCosinex, cameraCosiney);
        const auto fov = std::max(fov_x, fov_y);

        const auto dd = fov / resolution / zoom;
        const auto dd_start = -dd * resolution / 2 + dd / 2;

        std::vector<T> data(resolution * resolution * 4, T { 0 });
        for (std::size_t i = 3; i < data.size(); i = i + 4) {
            data[i] = 1;
        }

        auto items = world.getItemPointers();
        std::sort(items.begin(), items.end());
        std::vector<std::array<T, 3>> colors(items.size());
        for (std::size_t i = 0; i < colors.size(); ++i) {
            const auto deg = (std::numbers::pi_v<T> * i) / colors.size();
            colors[i] = HSVtoRGB(deg);
        }

        RandomState state;
        constexpr int n_samples = 4;

        for (std::size_t y = 0; y < resolution; ++y) {
            const auto y_flat = y * resolution;
            for (std::size_t x = 0; x < resolution; ++x) {
                const auto flat = (y_flat + x) * 4;
                for (std::size_t i = 0; i < n_samples; ++i) {
                    Particle<T> p = { .pos = cameraPos, .dir = dir };
                    p.dir = vectormath::rotate(p.dir, cameraCosiney, dd_start + dd * x + state.randomUniform<T>(-dd, dd) / 2);
                    p.dir = vectormath::rotate(p.dir, cameraCosinex, dd_start + dd * y + state.randomUniform<T>(-dd, dd) / 2);
                    const auto r = world.intersectVisualization(p);
                    if (r.valid()) {
                        const auto colorIdx = std::lower_bound(items.cbegin(), items.cend(), r.item);
                        if (colorIdx != items.cend()) {
                            const auto& c = colors[colorIdx - items.cbegin()];
                            p.translate(r.intersection);
                            const auto scaling = vectormath::dot(p.dir, r.normal);
                            for (std::size_t i = 0; i < 3; ++i) {
                                data[flat + i] += std::min(c[i] * scaling, T { 1 }) / n_samples;
                            }
                        }
                    }
                }
            }
        }
        for (std::size_t i = 0; i < data.size(); i = i + 4) {
            T s = 0;
            for (auto j = i; j < i + 3; ++j)
                s += data[j];
            if (s == 0)
                for (auto j = i; j < i + 3; ++j)
                    data[j] = 1;
        }
        return data;
    }
}

}