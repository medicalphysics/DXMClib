﻿# DXMClib
C++ library for x-ray dose scoring in voxel geometry. 

DXMClib aims to be an easy to use C++ dose scoring library for voxelized geometry. It is the primary simulation engine of [OpenDXMC](https://github.com/medicalphysics/OpenDXMC), a GUI application for Monte Carlo dose simulation for diagnostic beam energy levels. 

Documentation can be found at https://dxmclib.readthedocs.io/.

It is possible to simulate dose from conventional x-ray and CT examinations in arbitrary materials. DXMClib also includes a x-ray specter generator based on the work by Gavin Poludniowski and Phil Evans; [Calculation of x‐ray spectra emerging from an x‐ray tube. Part I. Electron penetration characteristics in x‐ray targets](https://aapm.onlinelibrary.wiley.com/doi/abs/10.1118/1.2734725) and [Calculation of x‐ray spectra emerging from an x‐ray tube. Part II. X‐ray production and filtration in x‐ray targets](https://aapm.onlinelibrary.wiley.com/doi/abs/10.1118/1.2734726).

DXMClib is dependent on the [xraylib library by Tom Schoonjans](https://github.com/tschoonj/xraylib) for interaction cross sections and atomic form factors. 

### Installation
To compile DXMClib the xraylib library must be installed, both source and binaries for Windows and Linux are provided by the developer. [See the installation instructions for xraylib.](https://github.com/tschoonj/xraylib/wiki/Installation-instructions)

DXMClib uses CMake as build generator, to include DXMClib in a CMake project it is recommended to use CMakes 'FetchContent' module. Example to include DXMClib in your CMakeLists.txt:

    include(FetchContent)
    ## Adding DXMClib package
    FetchContent_Declare(
        libdxmc
        GIT_REPOSITORY https://github.com/medicalphysics/DXMClib.git
        )
    FetchContent_MakeAvailable(libdxmc)

    # Example target you develop
    add_executable(your_executable your_example.cpp)
    # Adding DXMClib headers
    target_include_directories(your_executable PRIVATE ${libdxmc_SOURCE_DIR}/include)
    # Linking to DXMClib
    target_link_libraries(your_executable PRIVATE libdxmc)

dxmclib takes advantage of concepts and std::atomic_ref introduced in C++20. Currently this library is tested on MSVC >= 16.8 and Clang >= 13.0 compilers.   