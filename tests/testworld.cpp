

#include "dxmc/transport.h"
#include "dxmc/world.h"
#include <iostream>

inline std::size_t index(std::size_t i, std::size_t j, std::size_t k, std::array<std::size_t, 3> dim)
{
    return i + dim[0] * j + dim[0] * dim[1] * k;
}

bool testCTDIgeneration(void)
{
    CTDIPhantom w;

    w.validate();
    bool valid = w.isValid();
    return !valid;
}

int main(int argc, char* argv[])
{
    bool success = true;
    success = success && testCTDIgeneration();

    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
