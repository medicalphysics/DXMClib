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

Copyright 2019 Erlend Andersen
*/

#pragma once

#include "dxmc/floating.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <type_traits>
#include <utility>

// Header library for simple 3D vector math
namespace dxmc {
namespace vectormath {

    template <typename T>
    concept Index = std::is_integral_v<T> && std::is_same<bool, T>::value == false;

    template <Floating T>
    inline constexpr T lenght_sqr(const std::array<T, 3>& vec) noexcept
    {
        constexpr T lsqr = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
        return lsqr;
    }

    template <Floating T>
    inline T lenght(const std::array<T, 3>& vec) noexcept
    {
        const T lsqr = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
        return std::sqrt(lsqr);
    }

    template <Floating T>
    std::pair<std::array<T, 3>, std::array<T, 3>> splice(const std::array<T, 6>& a)
    {
        std::array a1 { a[0], a[1], a[2] };
        std::array a2 { a[3], a[4], a[5] };
        return std::make_pair(a1, a2);
    }
    template <Floating T>
    std::array<T, 6> join(const std::array<T, 3>& a, const std::array<T, 3>& b)
    {
        std::array r { a[0], a[1], a[2], b[0], b[1], b[2] };
        return r;
    }

    template <Floating T>
    inline void normalize(T vec[3]) noexcept
    {
        const T lsqr = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
        constexpr T one { 1 };
        const T norm = one / std::sqrt(lsqr);
        vec[0] *= norm;
        vec[1] *= norm;
        vec[2] *= norm;
    }

    template <Floating T>
    inline void normalize(std::array<T, 3>& vec) noexcept
    {
        const T lsqr = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
        constexpr T one { 1 };
        const T norm = one / std::sqrt(lsqr);
        vec[0] *= norm;
        vec[1] *= norm;
        vec[2] *= norm;
    }

    template <Floating T>
    inline std::array<T, 3> add(const std::array<T, 3>& v1, const std::array<T, 3>& v2)
    {
        std::array<T, 3> r { v1[0] + v2[0], v1[1] + v2[1], v1[2] + v2[2] };
        return r;
    }
    template<Floating T>
    inline auto add( std::array<T, 3> v1, T v2)
    {
        v1[0] += v2;
        v1[1] += v2;
        v1[2] += v2;
        return v1;
    }
    template <Floating T>
    inline auto add(T v2 , const std::array<T, 3>& v1)
    {
        return add(v1, v2); 
    }
    template <Floating T>
    inline std::array<T, 3> subtract(const std::array<T, 3>& v1, const std::array<T, 3>& v2)
    {
        std::array<T, 3> r { v1[0] - v2[0], v1[1] - v2[1], v1[2] - v2[2] };
        return r;
    }

    template <Floating T>
    inline T dot(const T v1[3], const T v2[3]) noexcept
    {
        return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    }

    template <Floating T>
    inline T dot(const std::array<T, 3>& v1, const std::array<T, 3>& v2) noexcept
    {
        return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    }

    template <Floating T>
    inline std::array<T, 3> cross(const std::array<T, 3>& v1, const std::array<T, 3> v2) noexcept
    {
        std::array<T, 3> res {
            v1[1] * v2[2] - v1[2] * v2[1],
            v1[2] * v2[0] - v1[0] * v2[2],
            v1[0] * v2[1] - v1[1] * v2[0]
        };
        return res;
    }

    template <Floating T>
    inline std::array<T, 3> cross(const std::array<T, 6>& v1) noexcept
    {
        std::array<T, 3> res {
            v1[1] * v1[5] - v1[2] * v1[4],
            v1[2] * v1[3] - v1[0] * v1[5],
            v1[0] * v1[4] - v1[1] * v1[3]
        };
        return res;
    }

    template <Floating T>
    inline void rotate(T vec[3], const T axis[3], const T angle) noexcept
    {
        const T sang = std::sin(angle);
        const T cang = std::cos(angle);
        constexpr T one { 1 };
        const T midt = (one - cang) * dot(vec, axis);

        T out[3];
        out[0] = cang * vec[0] + midt * axis[0] + sang * (axis[1] * vec[2] - axis[2] * vec[1]);
        out[1] = cang * vec[1] + midt * axis[1] + sang * (-axis[0] * vec[2] + axis[2] * vec[0]);
        out[2] = cang * vec[2] + midt * axis[2] + sang * (axis[0] * vec[1] - axis[1] * vec[0]);

        vec[0] = out[0];
        vec[1] = out[1];
        vec[2] = out[2];
    }
    template <Floating T>
    inline void rotate(T vec[3], const std::array<T, 3>& axis, const T angle) noexcept
    {
        const T sang = std::sin(angle);
        const T cang = std::cos(angle);
        constexpr T one { 1 };
        const T midt = (one - cang) * dot(vec, &axis[0]);

        std::array<T, 3> res {
            cang * vec[0] + midt * axis[0] + sang * (axis[1] * vec[2] - axis[2] * vec[1]),
            cang * vec[1] + midt * axis[1] + sang * (-axis[0] * vec[2] + axis[2] * vec[0]),
            cang * vec[2] + midt * axis[2] + sang * (axis[0] * vec[1] - axis[1] * vec[0])
        };
        vec[0] = res[0];
        vec[1] = res[1];
        vec[2] = res[2];
    }
    template <Floating T>
    inline std::array<T, 3> rotate(const std::array<T, 3>& vec, const std::array<T, 3>& axis, const T angle) noexcept
    {
        const T sang = std::sin(angle);
        const T cang = std::cos(angle);
        constexpr T one { 1 };
        const T midt = (one - cang) * dot(vec, axis);

        std::array<T, 3> res {
            cang * vec[0] + midt * axis[0] + sang * (axis[1] * vec[2] - axis[2] * vec[1]),
            cang * vec[1] + midt * axis[1] + sang * (-axis[0] * vec[2] + axis[2] * vec[0]),
            cang * vec[2] + midt * axis[2] + sang * (axis[0] * vec[1] - axis[1] * vec[0])
        };
        return res;
    }

    template <Floating T>
    inline T angleBetween(const std::array<T, 3>& vec1, const std::array<T, 3>& vec2) noexcept
    {
        // Herons formula for numeric stable angle computation
        // Do not edit parenthesis and such

        const auto vec3 = subtract(vec1, vec2);
        const T a = lenght(vec1);
        const T b = lenght(vec2);
        const T c = lenght(vec3);

        const T u = b >= c ? c - (a - b) : b - (a - c);

        const T nom = ((a - b) + c) * u;
        const T den = (a + (b + c)) * ((a - c) + b);
        return T { 2 } * std::atan(std::sqrt(nom / den));
    }

    template <Floating T>
    inline T angleBetweenOnPlane(const std::array<T, 3>& vec1, const std::array<T, 3>& vec2, const std::array<T, 3>& planeNormal) noexcept
    {
        const auto cr = cross(vec1, vec2);
        return std::atan2(dot(cr, planeNormal), dot(vec1, vec2));
    }

    template <Index U = std::size_t, Floating T>
    inline U argmin3(const std::array<T, 3>& vec) noexcept
    {
        const T x = std::abs(vec[0]);
        const T y = std::abs(vec[1]);
        const T z = std::abs(vec[2]);
        return x <= y ? x <= z ? 0 : 2 : y <= z ? 1
                                                : 2;
    }

    template <Index U = std::size_t, Floating T>
    inline U argmax3(const std::array<T, 3>& vec) noexcept
    {
        const T x = std::abs(vec[0]);
        const T y = std::abs(vec[1]);
        const T z = std::abs(vec[2]);
        return x >= y ? x >= z ? 0 : 2 : y >= z ? 1
                                                : 2;
    }

    template <Floating T>
    inline std::array<T, 3> changeBasis(const std::array<T, 3>& b1, const std::array<T, 3>& b2, const std::array<T, 3>& b3, const std::array<T, 3>& vector) noexcept
    {
        std::array<T, 3> res {
            b1[0] * vector[0] + b2[0] * vector[1] + b3[0] * vector[2],
            b1[1] * vector[0] + b2[1] * vector[1] + b3[1] * vector[2],
            b1[2] * vector[0] + b2[2] * vector[1] + b3[2] * vector[2]
        };
        return res;
    }

    template <Floating T>
    inline std::array<T, 3> changeBasisInverse(const std::array<T, 3>& b1, const std::array<T, 3>& b2, const std::array<T, 3>& b3, const std::array<T, 3>& vector) noexcept
    {
        std::array<T, 3> res {
            b1[0] * vector[0] + b1[1] * vector[1] + b1[2] * vector[2],
            b2[0] * vector[0] + b2[1] * vector[1] + b2[2] * vector[2],
            b3[0] * vector[0] + b3[1] * vector[1] + b3[2] * vector[2]
        };
        return res;
    }

    template <Floating T>
    inline void peturb(std::array<T, 3>& vec, const T theta, const T phi) noexcept
    {
        // rotates a unit vector theta degrees from its current direction
        // phi degrees about a arbitrary axis orthogonal to the direction vector

        // First we find a vector orthogonal to the vector direction
        // T vec_xy[3], k[3] = { 0, 0, 0 };

        const auto minInd = argmin3<std::uint_fast32_t, T>(vec);
        std::array<T, 3> k { 0, 0, 0 };
        k[minInd] = T { 1 };

        const auto vec_xy_raw = vectormath::cross(vec, k);

        // rotating the arbitrary orthogonal axis about vector direction
        const auto vec_xy = rotate(vec_xy_raw, vec, phi);

        // rotating vector about theta with vector addition (they are orthonormal)
        const T tsin = std::sin(theta);
        const T tcos = std::cos(theta);
        vec[0] = vec[0] * tcos + vec_xy[0] * tsin;
        vec[1] = vec[1] * tcos + vec_xy[1] * tsin;
        vec[2] = vec[2] * tcos + vec_xy[2] * tsin;
    }

}
}