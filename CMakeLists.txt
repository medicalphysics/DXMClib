# Almost all CMake files should start with this
# You should always specify a range with the newest
# and oldest tested versions of CMake. This will ensure
# you pick up the best policies.
cmake_minimum_required(VERSION 3.10...3.20)

# This is your project statement. You should always list languages;
# Listing the version is nice here since it sets lots of useful variables
project(dxmclib VERSION 1.0 LANGUAGES CXX)

# If you set any CMAKE_ variables, that can go here.
# (But usually don't do this, except maybe for C++ standard)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set_property(GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS)

# Find packages go here.
find_package(Threads REQUIRED) # for threading support


# Finding xraylib library
IF(WIN32)
    set(XRAYLIB_INSTALL_PATH "C:/Program Files/xraylib 64-bit" CACHE PATH "Install path for xraylib (https://github.com/tschoonj/xraylib)")
    set(XRAYLIB_INCLUDES "${XRAYLIB_INSTALL_PATH}/Include")
	set(XRAYLIB_LIBRARY_SHARED_NAME "libxrl-7.dll")
	set(XRAYLIB_LIBRARY_SHARED_PATH "C:/Windows/sysnative/${XRAYLIB_LIBRARY_SHARED_NAME}")
    set(XRAYLIB_LIBRARY_SHARED_INSTALL_PATH "C:/Windows/system32/${XRAYLIB_LIBRARY_SHARED_NAME}")
    find_library(XRAYLIB_LIBRARY NAMES libxrl.lib libxrl libxrl-7 xrl HINTS "${XRAYLIB_INSTALL_PATH}/LIB")

    add_library(xraylib SHARED IMPORTED)
    set_property(TARGET xraylib PROPERTY IMPORTED_LOCATION "${XRAYLIB_LIBRARY}" "${XRAYLIB_LIBRARY_SHARED_INSTALL_PATH}")

ENDIF(WIN32)
IF(UNIX AND NOT APPLE)
    set(XRAYLIB_INSTALL_PATH "/usr/local/lib" CACHE PATH "Install path for xraylib (https://github.com/tschoonj/xraylib)")
    set(XRAYLIB_INCLUDES "/usr/local/include/xraylib")
    find_library(XRAYLIB_LIBRARY NAMES libxrl libxrl-7 xrl HINTS "${XRAYLIB_INSTALL_PATH}")
ENDIF()


# You should usually split this into folders, but this is a simple example

# This is a "default" library, and will match the *** variable setting.
# Other common choices are STATIC, SHARED, and MODULE
# Including header files here helps IDEs but is not required.
# Output libname matches target name, with the usual extensions on your system

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

target_include_directories(libdxmc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(libdxmc PRIVATE ${XRAYLIB_INCLUDES})

set_target_properties(libdxmc PROPERTIES PUBLIC_HEADER "${HEADER_LIST}")


# Make sure you link your targets with this command. It can also link libraries and
# even flags, so linking a target that does not exist will not give a configure-time error.
target_link_libraries(libdxmc PRIVATE 
    Threads::Threads
    xraylib
    )
include(GNUInstallDirs)
install(TARGETS libdxmc
    #EXPORT libdxmc
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dxmclib
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dxmclib
    PRIVATE_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dxmclib
    )
