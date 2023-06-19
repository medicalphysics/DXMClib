

## setting up lodepng library files
set(LODEPNG_FILES 
    "${dxmclib_SOURCE_DIR}/thirdparty/lodepng/include/dxmclodepngwrapper.hpp"
    "${dxmclib_SOURCE_DIR}/thirdparty/lodepng/src/dxmclodepngwrapper.cpp"
    "${dxmclib_SOURCE_DIR}/thirdparty/lodepng/src/lodepng.cpp"
    "${dxmclib_SOURCE_DIR}/thirdparty/lodepng/src/lodepng.h"
)

add_library(lodepngwrapper STATIC ${LODEPNG_FILES})
target_include_directories(lodepngwrapper PUBLIC ${dxmclib_SOURCE_DIR}/thirdparty/lodepng/include)
