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

#include "dxmc/floating.h"

#include <cmath>
#include <cstdint>
//Header library for simple 3D vector math
namespace dxmc {
namespace vectormath {
    //template <typename T>
    //concept Floating = std::is_floating_point_v<T>;

    template <typename T>
    concept Index = std::is_integral_v<T>&& std::is_same<bool, T>::value == false;

    template <Floating T>
    inline T lenght_sqr(T vec[3])
    {
        const T lsqr = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
        return lsqr;
    }

    template <Floating T>
    inline T lenght(T vec[3])
    {
        const T lsqr = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
        return std::sqrt(lsqr);
    }

    template <Floating T>
    inline void normalize(T vec[3])
    {
        const T lsqr = vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2];
        constexpr T one { 1 };
        const T norm = one / std::sqrt(lsqr);
        vec[0] *= norm;
        vec[1] *= norm;
        vec[2] *= norm;
    }

    template <Floating T>
    inline T dot(const T v1[3], const T v2[3])
    {
        return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    }

    template <Floating T>
    inline void cross(const T v1[3], const T v2[3], T res[3])
    {
        res[0] = v1[1] * v2[2] - v1[2] * v2[1];
        res[1] = v1[2] * v2[0] - v1[0] * v2[2];
        res[2] = v1[0] * v2[1] - v1[1] * v2[0];
    }
    template <Floating T>
    inline void cross(const T v1[6], T res[3])
    {
        res[0] = v1[1] * v1[5] - v1[2] * v1[4];
        res[1] = v1[2] * v1[3] - v1[0] * v1[5];
        res[2] = v1[0] * v1[4] - v1[1] * v1[3];
    }

    template <Floating T>
    inline void rotate(T vec[3], const T axis[3], const T angle)
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
    inline void projectToPlane(T vec[3], const T planeNormal[3])
    {
        const T d = dot(vec, planeNormal);
        for (std::size_t i = 0; i < 3; ++i)
            vec[i] = vec[i] - d * planeNormal[i];
    }

    template <Floating T>
    inline T angleBetween(const T vec1[3], const T vec2[3])
    {
        // Herons formula for numeric stable angle computation
        // Do not edit parenthesis and such
        const T vec3[3] = { vec1[0] - vec2[0], vec1[1] - vec2[1], vec1[2] - vec2[2] };
        const T a = lenght(vec1);
        const T b = lenght(vec2);
        const T c = lenght(vec3);

        const T u = b >= c ? c - (a - b) : b - (a - c);

        const T nom = ((a - b) + c) * u;
        const T den = (a + (b + c)) * ((a - c) + b);
        return T { 2.0 } * std::atan(std::sqrt(nom / den));
    }

    template <Floating T>
    inline double angleBetweenOnPlane(T vec1[3], T vec2[3], T planeNormal[3])
    {
        normalize(vec1);
        normalize(vec2);
        normalize(planeNormal);

        T cr[3];
        cross(vec1, vec2, cr);
        return std::atan2(dot(cr, planeNormal), dot(vec1, vec2));
    }

    template <Index U, Floating T>
    inline U argmin3(const T vec[3])
    {
        const T x = std::abs(vec[0]);
        const T y = std::abs(vec[1]);
        const T z = std::abs(vec[2]);
        return x < y ? x < z ? 0 : 2 : y < z ? 1 : 2;
    }

    template <Index U, Floating T>
    inline U argmax3(const T vec[3])
    {
        const T x = std::abs(vec[0]);
        const T y = std::abs(vec[1]);
        const T z = std::abs(vec[2]);
        return x > y ? x > z ? 0 : 2 : y > z ? 1 : 2;
    }

    template <Floating T>
    inline void changeBasis(const T b1[3], const T b2[3], const T b3[3], const T vector[3], T newVector[3])
    {
        newVector[0] = b1[0] * vector[0] + b2[0] * vector[1] + b3[0] * vector[2];
        newVector[1] = b1[1] * vector[0] + b2[1] * vector[1] + b3[1] * vector[2];
        newVector[2] = b1[2] * vector[0] + b2[2] * vector[1] + b3[2] * vector[2];
    }
    template <Floating T>
    inline void changeBasis(const T b1[3], const T b2[3], const T b3[3], T vector[3])
    {
        T newVector[3];
        newVector[0] = b1[0] * vector[0] + b2[0] * vector[1] + b3[0] * vector[2];
        newVector[1] = b1[1] * vector[0] + b2[1] * vector[1] + b3[1] * vector[2];
        newVector[2] = b1[2] * vector[0] + b2[2] * vector[1] + b3[2] * vector[2];
        vector[0] = newVector[0];
        vector[1] = newVector[1];
        vector[2] = newVector[2];
    }

    template <Floating T>
    inline void changeBasisInverse(const T b1[3], const T b2[3], const T b3[3], const T vector[3], T newVector[3])
    {
        newVector[0] = b1[0] * vector[0] + b1[1] * vector[1] + b1[2] * vector[2];
        newVector[1] = b2[0] * vector[0] + b2[1] * vector[1] + b2[2] * vector[2];
        newVector[2] = b3[0] * vector[0] + b3[1] * vector[1] + b3[2] * vector[2];
    }

    template <Floating T>
    inline void changeBasisInverse(const T b1[3], const T b2[3], const T b3[3], T vector[3])
    {
        T newVector[3];
        newVector[0] = b1[0] * vector[0] + b1[1] * vector[1] + b1[2] * vector[2];
        newVector[1] = b2[0] * vector[0] + b2[1] * vector[1] + b2[2] * vector[2];
        newVector[2] = b3[0] * vector[0] + b3[1] * vector[1] + b3[2] * vector[2];
        vector[0] = newVector[0];
        vector[1] = newVector[1];
        vector[2] = newVector[2];
    }

    template <Floating T>
    inline void peturb(T vec[3], const T theta, const T phi)
    {
        // rotates a unit vector theta degrees from its current direction
        // phi degrees about a arbitrary axis orthogonal to the direction vector

        // First we find a vector orthogonal to the vector direction
        T vec_xy[3], k[3] = { 0, 0, 0 };

        const auto minInd = argmin3<std::uint_fast32_t, T>(vec);

        k[minInd] = T { 1.0 };

        vectormath::cross(vec, k, vec_xy);

        // rotating the arbitrary orthogonal axis about vector direction
        rotate(vec_xy, vec, phi);

        //rotating vector about theta with vector addition (they are orthonormal)
        const T tsin = std::sin(theta);
        const T tcos = std::cos(theta);
        vec[0] = vec[0] * tcos + vec_xy[0] * tsin;
        vec[1] = vec[1] * tcos + vec_xy[1] * tsin;
        vec[2] = vec[2] * tcos + vec_xy[2] * tsin;
    }

}
}