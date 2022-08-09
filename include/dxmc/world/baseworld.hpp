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
#include "dxmc/world/kdtree.hpp"
#include "dxmc/material.hpp"
#include "dxmc/attenuationlut.hpp"

#include <vector>
#include <array>


namespace dxmc {

template<Floating T>
struct WorldResult {


};


template <Floating T>
class BaseWorld {
public:
    virtual void setMaxPhotonEnergy(const T maxEnergy) = 0;


protected:
private:
    std::array<T, 6> m_aabb = { -1, -1, -1, 1, 1, 1 };
    AttenuationLut<T> m_attenuationLut;
    std::vector<Material<T>> m_materials;
    
};
}