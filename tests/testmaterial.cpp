
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


#include "dxmc/material/atomicshell.hpp"



using namespace dxmc;

bool testatomicshell()
{
    AtomicShell<float> shell;
    return true;
}


int main(int argc, char* argv[])
{
    auto success = testatomicshell();
    if (success)
        return EXIT_SUCCESS;
    return EXIT_FAILURE;
}
