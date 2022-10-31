

# EPICS data parser
option(DXMCLIB_BUILD_EPICS_DATA_PARSER "Build EPICS data parser" ON)
if(DXMCLIB_BUILD_EPICS_DATA_PARSER)
    SET(DXMCLIB_EPICS_DATA_DIRPATH CACHE PATH "Path to directory holding EPDL2017.ALL and EADL2017.ALL in ENDL format")
    if(DXMCLIB_EPICS_DATA_DIRPATH)
        SET(DXMCLIB_EPICS_DATA_EPDLPATH "${DXMCLIB_EPICS_DATA_DIRPATH}/EPDL2017.ALL" CACHE FILEPATH "Path to EPDL2017.ALL" FORCE)
        SET(DXMCLIB_EPICS_DATA_EADLPATH "${DXMCLIB_EPICS_DATA_DIRPATH}/EADL2017.ALL" CACHE FILEPATH "Path to EADL2017.ALL" FORCE)
        message("Using EPDL file: ${DXMCLIB_EPICS_DATA_EPDLPATH}")
        message("Using EADL file: ${DXMCLIB_EPICS_DATA_EADLPATH}")
    endif()
    add_executable(EPICSparser 
        EPICSParser/atomicshell.cpp 
        EPICSParser/atomicelement.cpp 
        EPICSParser/epicsparser.cpp 
        EPICSParser/main.cpp)
    target_compile_definitions(EPICSparser PUBLIC EADLPATH="${DXMCLIB_EPICS_DATA_EADLPATH}")
    target_compile_definitions(EPICSparser PUBLIC EPDLPATH="${DXMCLIB_EPICS_DATA_EPDLPATH}")
    target_include_directories(EPICSparser PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/EPICSParser/include")    
endif()


