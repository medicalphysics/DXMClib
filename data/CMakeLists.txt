

# EPICS data parser
option(BUILD_EPICS_DATA_PARSER "Build EPICS data parser" ON)
if(BUILD_EPICS_DATA_PARSER)
    SET(EPICS_DATA_DIRPATH CACHE PATH "Path to directory holding EPDL2017.ALL and EADL2017.ALL in ENDL format")
    if(EPICS_DATA_DIRPATH)
        SET(EPICS_DATA_EPDLPATH "${EPICS_DATA_DIRPATH}/EPDL2017.ALL" CACHE FILEPATH "Path to EPDL2017.ALL" FORCE)
        SET(EPICS_DATA_EADLPATH "${EPICS_DATA_DIRPATH}/EADL2017.ALL" CACHE FILEPATH "Path to EADL2017.ALL" FORCE)
        message("Using EPDL file: ${EPICS_DATA_EPDLPATH}")
        message("Using EADL file: ${EPICS_DATA_EADLPATH}")
    endif()
    add_executable(EPICSparser atomicshell.cpp atomicelement.cpp epicsparser.cpp main.cpp)
    target_compile_definitions(EPICSparser PUBLIC EADLPATH="${EPICS_DATA_EADLPATH}")
    target_compile_definitions(EPICSparser PUBLIC EPDLPATH="${EPICS_DATA_EPDLPATH}")
    target_include_directories(EPICSparser PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
    target_link_libraries(EPICSparser PRIVATE libdxmc)
endif()


