
#pragma once

#include "dxmc/floating.h"

namespace dxmc {
/**
 * @brief Simple struct to describe a photon 
*/
template <Floating T>
struct Particle {
    /**
     * @brief Position vector in three dimensions.
    */
    T pos[3];
    /**
     * @brief Direction vector in three dimension. This vector is threated as a normal vector.
    */
    T dir[3];
    /**
     * @brief Photon energy in keV.
    */
    T energy;
    /**
     * @brief Photon relative weight.
    */
    T weight;
};
}