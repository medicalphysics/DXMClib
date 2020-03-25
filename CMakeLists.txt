# Almost all CMake files should start with this
# You should always specify a range with the newest
# and oldest tested versions of CMake. This will ensure
# you pick up the best policies.
cmake_minimum_required(VERSION 3.10...3.20)

# This is your project statement. You should always list languages;
# Listing the version is nice here since it sets lots of useful variables
project(dxmclib VERSION 1.1 LANGUAGES CXX)

# If you set any CMAKE_ variables, that can go here.
# (But usually don't do this, except maybe for C++ standard)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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
    "${dxmclib_SOURCE_DIR}/include/dxmc/attenuationlut.h"
    "${dxmclib_SOURCE_DIR}/include/dxmc/beamfilters.h"
    "${dxmclib_SOURCE_DIR}/include/dxmc/dxmcrandom.h"
    "${dxmclib_SOURCE_DIR}/include/dxmc/exposure.h"
    "${dxmclib_SOURCE_DIR}/include/dxmc/material.h"
    "${dxmclib_SOURCE_DIR}/include/dxmc/progressbar.h"
    "${dxmclib_SOURCE_DIR}/include/dxmc/source.h"
    "${dxmclib_SOURCE_DIR}/include/dxmc/transport.h"
    "${dxmclib_SOURCE_DIR}/include/dxmc/tube.h"
    "${dxmclib_SOURCE_DIR}/include/dxmc/vectormath.h"
    "${dxmclib_SOURCE_DIR}/include/dxmc/world.h"
    )

add_library(libdxmc 
    src/world.cpp
    src/source.cpp
    src/exposure.cpp
    src/tube.cpp
    src/dxmcrandom.cpp
    src/beamfilters.cpp
    src/material.cpp
    src/transport.cpp
	src/attenuationlut.cpp
    ${HEADER_LIST}
	)
target_include_directories(libdxmc PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")
target_include_directories(libdxmc PRIVATE ${XRAYLIB_INCLUDES})

option(DXMC_USE_LOWENERGY_COMPTON "Use low energy correction for compton interactions (increases computation time)" ON)
if(${USE_LOWENERGY_COMPTON})
    target_compile_definitions(libdxmc PUBLIC DXMC_USE_LOWENERGY_COMPTON)
endif()

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
    if(Doxygen_FOUND)
	    add_subdirectory(docs)
    else()
	    message(STATUS "Doxygen not found, not building docs")
    endif()
endif()

# Examples
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    option(BUILD_EXAMPLES "Build examples" ON)
    if(BUILD_EXAMPLES)
        add_subdirectory(examples)
    endif()
endif()