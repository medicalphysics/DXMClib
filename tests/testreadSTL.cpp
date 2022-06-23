

#include "dxmc/world/triangulatedmesh.hpp"

#include <iostream>

int main(int argc, char* argv[])
{
    dxmc::STLReader<float> reader;

    auto mesh = reader("untitled_bin.stl");
    std::cout << "Message: " << reader.message() << std::endl;

    const auto& faceIdx1 = mesh.getFaceIndex(0);

    for (const auto idx : faceIdx1) {
        std::cout << "Vertice: ";
        for (auto v : mesh.getVertice(idx)) {
            std::cout << v << ", ";
        }
        std::cout << std::endl;
    }

    dxmc::Particle<float> p;
    p.pos = { .5, .2, -.3 };
    p.dir = { 0, 0, 1 };
    for (std::size_t i = 0; i < mesh.nFaces(); i++) {
        auto res = mesh.intersect<-1>(p, i);
        if (res) {
            std::cout << i << " " << *res << std::endl;
        } else {
            std::cout << i << " none" << std::endl;
        }
    }

    return EXIT_SUCCESS;
}
