# Almost all CMake files should start with this
# You should always specify a range with the newest
# and oldest tested versions of CMake. This will ensure
# you pick up the best policies.
cmake_minimum_required(VERSION 3.10...3.20)

# This is your project statement. You should always list languages;
# Listing the version is nice here since it sets lots of useful variables
project(dxmclib VERSION 0.7.0 LANGUAGES CXX)

# If you set any CMAKE_ variables, that can go here.
# (But usually don't do this, except maybe for C++ standard)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")

#Function for adding extra compiler flags if available
include(CheckCXXCompilerFlag)
function(enable_cxx_compiler_flag_if_supported flag)
    string(FIND "${CMAKE_CXX_FLAGS}" "${flag}" flag_already_set)
    if(flag_already_set EQUAL -1)
        check_cxx_compiler_flag("${flag}" flag_supported)
        if(flag_supported)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}" PARENT_SCOPE)
        endif()
        unset(flag_supported CACHE)
    endif()
endfunction()
# We want to add extra compiler warnings if supported
#enable_cxx_compiler_flag_if_supported("-Wall")
#enable_cxx_compiler_flag_if_supported("-Wextra")
#enable_cxx_compiler_flag_if_supported("-pedantic")

option(OPTIMIZE_FOR_NATIVE "Build with -march=native on gcc/clang or AVX_x on MSVC " ON)
if(OPTIMIZE_FOR_NATIVE)
    if(MSVC)
        set(MSVC_OPTION "AVX2" CACHE STRING "Machine instructions support")
        set_property(CACHE MSVC_OPTION PROPERTY STRINGS AVX AVX2 AVX512)
        MESSAGE("-- Trying to build for ${MSVC_OPTION} instruction set")
        enable_cxx_compiler_flag_if_supported("-arch:${MSVC_OPTION}")  
    endif()
    if(UNIX)
        MESSAGE("-- Trying to build for native (-march=native)")
        enable_cxx_compiler_flag_if_supported("-march=native")  
    endif()
endif()

# Find packages go here.
find_package(Threads REQUIRED) # for threading support

# Finding xraylib library
IF(WIN32)
    set(XRAYLIB_INSTALL_PATH "C:/Program Files/xraylib 64-bit" CACHE PATH "Install path for xraylib (https://github.com/tschoonj/xraylib)")
    set(XRAYLIB_INCLUDES "${XRAYLIB_INSTALL_PATH}/Include")
	set(XRAYLIB_LIBRARY_SHARED_NAME "libxrl-11.dll")
	set(XRAYLIB_LIBRARY_SHARED_PATH "C:/Windows/sysnative/${XRAYLIB_LIBRARY_SHARED_NAME}")
    set(XRAYLIB_LIBRARY_SHARED_PATH_64 "C:/Windows/System32/${XRAYLIB_LIBRARY_SHARED_NAME}")
    find_library(XRAYLIB NAMES libxrl.lib libxrl libxrl-11 xrl HINTS "${XRAYLIB_INSTALL_PATH}/LIB")
ENDIF(WIN32)
IF(UNIX AND NOT APPLE)
    set(XRAYLIB_INSTALL_PATH "/usr/local/lib" CACHE PATH "Install path for xraylib (https://github.com/tschoonj/xraylib)")
    set(XRAYLIB_INCLUDES "/usr/local/include/xraylib")
    find_library(XRAYLIB NAMES libxrl libxrl-11 xrl HINTS "${XRAYLIB_INSTALL_PATH}")
ENDIF()

set(HEADER_LIST 
    "${dxmclib_SOURCE_DIR}/include/dxmc/attenuationinterpolator.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/attenuationlut.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/beamfilters.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/dxmcrandom.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/exposure.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/material.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/progressbar.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/source.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/transport.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/tube.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/vectormath.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/triangulatedmesh.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/baseworld.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/kdtree.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world/intersectsimpleobjects.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/floating.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/constants.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/particle.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/interpolation.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/lowenergycorrectionmodel.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc/betheHeitlerCrossSection.hpp"
    "${dxmclib_SOURCE_DIR}/include/dxmc.hpp"    
    )

add_library(libdxmc 
    src/material.cpp    
    ${HEADER_LIST}
	)
target_include_directories(libdxmc PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_include_directories(libdxmc PRIVATE ${XRAYLIB_INCLUDES})

# Make sure you link your targets with this command. It can also link libraries and
# even flags, so linking a target that does not exist will not give a configure-time error.
target_link_libraries(libdxmc PRIVATE 
    Threads::Threads
    ${XRAYLIB}
    )

## Install
install(TARGETS libdxmc 
    EXPORT libdxmc
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    RUNTIME DESTINATION bin
    INCLUDES DESTINATION include
    )

# Due to windows not being windows on 34/64 bits and xraylib puts it's dll in system32 folder we try to install 
# the library from system32 folder regardless of cmake or NSIS being 64 or 32 bits
if(WIN32)
    install(FILES ${XRAYLIB_LIBRARY_SHARED_PATH} DESTINATION bin OPTIONAL)
    install(FILES ${XRAYLIB_LIBRARY_SHARED_PATH_64} DESTINATION bin OPTIONAL)
endif(WIN32)
install(FILES ${HEADER_LIST} DESTINATION include)

# Testing
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    ENABLE_TESTING()
    include(CTest)
endif()
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
    add_subdirectory(tests)
endif()

# Validation
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    add_subdirectory(validation)
endif()


# Docs by doxygen
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    find_package(Doxygen)
    find_package(Sphinx)
    if(Doxygen_FOUND AND Sphinx_FOUND)
	    add_subdirectory(docs)
    else()
	    message(STATUS "Doxygen or Sphinx not found, not building docs")
    endif()
endif()

# Examples
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    option(BUILD_EXAMPLES "Build examples" ON)
    if(BUILD_EXAMPLES)
        add_subdirectory(examples)
    endif()
endif()