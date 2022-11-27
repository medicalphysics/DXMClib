

# EPICS data parser
SET(DXMCLIB_EPICS_PHYSICSLISTS "${CMAKE_CURRENT_BINARY_DIR}/physicslists.bin" CACHE INTERNAL "Path to EPICS physicslists")
option(DXMCLIB_BUILD_EPICS_DATA_PARSER "Build EPICS data parser" ON)
if(DXMCLIB_BUILD_EPICS_DATA_PARSER)
    SET(DXMCLIB_EPICS_DATA_DIRPATH ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH "Path to directory holding EPDL2017.ALL and EADL2017.ALL in ENDL format")
    if(DXMCLIB_EPICS_DATA_DIRPATH)
        SET(DXMCLIB_EPICS_DATA_EPDLPATH "${DXMCLIB_EPICS_DATA_DIRPATH}/EPDL2017.ALL" CACHE FILEPATH "Path to EPDL2017.ALL" FORCE)
        SET(DXMCLIB_EPICS_DATA_EADLPATH "${DXMCLIB_EPICS_DATA_DIRPATH}/EADL2017.ALL" CACHE FILEPATH "Path to EADL2017.ALL" FORCE)

        SET(DXMCLIB_EPICS_DOWNLOAD OFF CACHE BOOL "Download EPICS data from https://www-nds.iaea.org/epics/ if necessary.")

        #CHECK IF EPICS DATAFILES EXISTS
        if(NOT EXISTS "${DXMCLIB_EPICS_DATA_EPDLPATH}" AND "${DXMCLIB_EPICS_DOWNLOAD}" )
            MESSAGE("Downloading EPDL")
            FILE(DOWNLOAD "https://www-nds.iaea.org/epics/ENDL2017/EPDL2017.ALL" "${DXMCLIB_EPICS_DATA_EPDLPATH}" SHOW_PROGRESS)
        endif()
        if(NOT EXISTS "${DXMCLIB_EPICS_DATA_EADLPATH}" AND "${DXMCLIB_EPICS_DOWNLOAD}" )
            MESSAGE("Downloading EADL")
            FILE(DOWNLOAD "https://www-nds.iaea.org/epics/ENDL2017/EADL2017.ALL" "${DXMCLIB_EPICS_DATA_EADLPATH}" SHOW_PROGRESS)
        endif()

        if(NOT EXISTS "${DXMCLIB_EPICS_DATA_EADLPATH}" OR NOT EXISTS "${DXMCLIB_EPICS_DATA_EPDLPATH}")
            message(FATAL_ERROR "Could not find data files EPDL2017.ALL or EPDL2017.ALL in ${DXMCLIB_EPICS_DATA_DIRPATH}")
        endif()

        message(STATUS "Using EPDL file: ${DXMCLIB_EPICS_DATA_EPDLPATH}")
        message(STATUS "Using EADL file: ${DXMCLIB_EPICS_DATA_EADLPATH}")
    endif()

    # Option to specify max energy range
    set(DXMCLIB_MAX_ENERGY_KEV "500" CACHE STRING "Max energy in keV for photon energies.")
    set(DXMCLIB_MIN_ENERGY_KEV "1")

    add_executable(EPICSparser         
        EPICSParser/atomicelementhandler.cpp 
        EPICSParser/epicsparser.cpp 
        EPICSParser/main.cpp        
    )    
    target_compile_definitions(EPICSparser PUBLIC DXMCLIB_MAXENERGY="${DXMCLIB_MAX_ENERGY_KEV}")
    target_compile_definitions(EPICSparser PUBLIC DXMCLIB_EADLPATH="${DXMCLIB_EPICS_DATA_EADLPATH}")
    target_compile_definitions(EPICSparser PUBLIC DXMCLIB_EPDLPATH="${DXMCLIB_EPICS_DATA_EPDLPATH}")
    target_compile_definitions(EPICSparser PUBLIC DXMCLIB_PHYSICSLISTSPATH="${DXMCLIB_EPICS_PHYSICSLISTS}")
    target_include_directories(EPICSparser PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/EPICSParser/include)     
    target_link_libraries(EPICSparser PRIVATE material )

    ## Copying Hartree-Fock profiles for pz=0 since this data is lacking from EPICS data
    add_custom_command(
        TARGET EPICSparser PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/EPICSParser/hartreeFockProfiles_0.txt
                ${CMAKE_CURRENT_BINARY_DIR}/hartreeFockProfiles_0.txt
    )

    ## make parsed EPICS data available to all binaries consuming dxmc on build machine
    target_compile_definitions(libdxmc PUBLIC DXMCLIB_PHYSICSLISTSPATH="${DXMCLIB_EPICS_PHYSICSLISTS}")

    ## make max energy exposed to libdxmc
    target_compile_definitions(libdxmc PUBLIC DXMCLIB_MAXENERGY="${DXMCLIB_MAX_ENERGY_KEV}")
    target_compile_definitions(libdxmc PUBLIC DXMCLIB_MINENERGY="${DXMCLIB_MIN_ENERGY_KEV}")

    add_custom_command(OUTPUT "${DXMCLIB_EPICS_PHYSICSLISTS}"
        COMMAND "${CMAKE_CURRENT_BINARY_DIR}/EPICSparser" 
        DEPENDS EPICSparser)

    add_custom_target(libdxmc_physics ALL
        DEPENDS "${DXMCLIB_EPICS_PHYSICSLISTS}")
    add_dependencies(libdxmc libdxmc_physics)
    # Adding physics data to install folder
    install(FILES ${DXMCLIB_EPICS_PHYSICSLISTS} DESTINATION bin)
endif()


