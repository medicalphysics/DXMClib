
#include "dxmc/particle.hpp"
#include "dxmc/world/baseworld.hpp"

#include <iostream>
int main(int argc, char* argv[])
{

    dxmc::BaseWorld<float> w;

    dxmc::Particle<float> p;
    p.pos = { 0, 0, 0 };
    p.dir = { 0, 0, 1 };

    const auto t=w.intersect(p);

    if (t) {
        std::cout << "Intersection at: " << *t << std::endl;
    } else {
        std::cout << "No intersection"<< std::endl;
    }


    return EXIT_SUCCESS;
}
