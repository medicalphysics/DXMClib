include(CTest)
ENABLE_TESTING()

set(TEST_SOURCES 
    "testworld.cpp"
#    "testdxmclib.cpp"
    "testtriangulatedmesh.cpp"
    "testkdtree.cpp"
    "teststatickdtree.cpp"
    "testworlditembase.cpp"
    "testinteractions.cpp"
#    "testtube.cpp"
#    "testsource.cpp"
#    "testbeamfilters.cpp"
    "testmaterial.cpp"
    "testrandom.cpp"
#    "testattenuationlut.cpp"
#    "testitemcollection.cpp"
#    "testvectormath.cpp"
#    "testexposure.cpp"
#    "testinterpolation.cpp"
)

foreach(sourcefile ${TEST_SOURCES})
    get_filename_component(filename ${sourcefile} NAME_WLE)
    add_executable(${filename} ${sourcefile})
    target_include_directories(${filename} PRIVATE ${XRAYLIB_INCLUDES})
    target_link_libraries(${filename} libdxmc  ${XRAYLIB})
    add_test(NAME ${filename} COMMAND ${filename})    
    
    ## Spesific tools for testmaterial 
    if("${filename}" STREQUAL "testmaterial")
        add_custom_command(
        TARGET ${filename} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/testmaterialAtomCoeffVis.py
                ${CMAKE_CURRENT_BINARY_DIR}/testmaterialAtomCoeffVis.py
        )
        add_custom_target(materialtestvis ALL
        DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/testmaterialAtomCoeffVis.py")
        add_dependencies(${filename} materialtestvis)
    endif()
endforeach()
