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
#include "dxmc/dxmcrandom.hpp"
#include "dxmc/floating.hpp"
#include "dxmc/particle.hpp"

#include <algorithm>
#include <optional>

namespace dxmc {

template <Floating T>
struct ResultObject {
};

// Forward declaration
template <Floating T>
class WorldItemBase;

template <Floating T>
struct IntersectionResult {
    const WorldItemBase<T>* item = nullptr;
    T intersection = 0;
};

template <Floating T>
class WorldItemBase {
public:
    virtual void translate(const std::array<T, 3>& dist) = 0;
    virtual std::array<T, 3> center() const = 0;
    virtual std::array<T, 6> AABB() const = 0;
    virtual IntersectionResult<T> intersect(const Particle<T>& p) const = 0;
    virtual IntersectionResult<T> intersect(const Particle<T>& p, const std::array<T, 2>& tbox) const = 0;
    virtual T transport(Particle<T>& p, RandomState& state) = 0;

    template <int FORWARD = 1>
    static std::optional<std::array<T, 2>> intersectAABB(const Particle<T>& p, const std::array<T, 6>& aabb)
    {
        auto t = []() -> std::array<T, 2> {
            if constexpr (FORWARD == 1)
                return std::array { T { 0 }, std::numeric_limits<T>::max() };
            else
                return std::array { std::numeric_limits<T>::lowest(), std::numeric_limits<T>::max() };
        }();

        for (std::uint_fast8_t i = 0; i < 3; i++) {
            if (std::abs(p.dir[i]) > std::numeric_limits<T>::epsilon()) {
                const auto d_inv = T { 1 } / p.dir[i];
                const auto [tmin, tmax] = std::minmax((aabb[i] - p.pos[i]) * d_inv, (aabb[i + 3] - p.pos[i]) * d_inv);
                t[0] = std::max(t[0], tmin);
                t[1] = std::min(t[1], tmax);
            }
        }
        return t[0] > t[1] ? std::nullopt : std::make_optional(t);
    }
    protected:
        

    private:
    };
}