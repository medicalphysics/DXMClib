include(CTest)
ENABLE_TESTING()

set(TEST_SOURCES 
    "testworld.cpp"
    "testbeams.cpp"
    "testctdiphantom.cpp"
    "testtriangulatedmesh.cpp"
    "testtetrahedalmesh.cpp"
    "testctaecfilter.cpp"
    "testworlditems.cpp" 
    "testbasicshapes.cpp"
    "testtube.cpp"
    "testvisualization.cpp"
    "testfluencescore.cpp"
    "testbowtiefilter.cpp"   
    "testrandom.cpp"
    "testtransport.cpp"
    "testvectormath.cpp"
    "testworldsphere.cpp"
)


foreach(sourcefile ${TEST_SOURCES})
    get_filename_component(filename ${sourcefile} NAME_WLE)
    add_executable(${filename} ${sourcefile})    
    target_link_libraries(${filename} libdxmc )
    add_test(NAME ${filename} COMMAND ${filename})    
    
    #Adding interaction binary data to example 
    dxmclib_add_physics_list(${filename})

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
