

#include "dxmc/world/triangulatedmesh.hpp"

#include <iostream>

int main(int argc, char* argv[])
{
    dxmc::STLReader<double> reader;

    auto mesh = reader("untitled.stl");
    std::cout << "Message: " << reader.message() << std::endl;

    auto triangles = mesh.getTriangles();

    for (auto triangle : triangles) {
        for (auto vert : triangle.vertices()) {
            for (auto p : vert) {
                std::cout << p << ", ";
            }
            std::cout << std::endl;
        }
        auto cent = triangle.calculateCenter();
        for (auto p : cent) {
            std::cout << p << ", ";
        }
        std::cout << std::endl;
        std::cout << std::endl;
    }

    

    return EXIT_SUCCESS;
}
