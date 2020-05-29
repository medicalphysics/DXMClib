
option(BUILD_DOCUMENTATION "Create and install the HTML based API documentation (requires Doxygen)" ${DOXYGEN_FOUND})

if(BUILD_DOCUMENTATION)
    
    set(doxyfile_in ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
    set(doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    configure_file(${doxyfile_in} ${doxyfile} @ONLY)

    add_custom_target(doc ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM)

    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION doc)

    IF(WIN32)
        add_custom_target(readme ALL
            COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/sphinx/make.bat html
            COMMENT "Generating help documentation"
            VERBATIM)
        install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/sphinx/build DESTINATION readme)
    ENDIF(WIN32)
    IF(UNIX AND NOT APPLE)
        add_custom_target(readme ALL
            COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/sphinx/make html
            COMMENT "Generating help documentation"
            VERBATIM)
        install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/sphinx/build/html DESTINATION readme)
    ENDIF()
endif()