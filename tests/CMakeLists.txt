


set(TEST_SOURCES 
    "testdxmclib.cpp"
    "testtriangulatedmesh.cpp"
#    "testkdtree.cpp"
#    "testtube.cpp"
#    "testsource.cpp"
#    "testbeamfilters.cpp"
#    "testmaterial.cpp"
#    "testrandom.cpp"
#    "testattenuationlut.cpp"
    "testitemcollection.cpp"
#    "testvectormath.cpp"
#    "testexposure.cpp"
#    "testinterpolation.cpp"
)

foreach(sourcefile ${TEST_SOURCES})
    get_filename_component(filename ${sourcefile} NAME_WLE)
    add_executable(${filename} ${sourcefile})
    target_include_directories(${filename} PRIVATE ${XRAYLIB_INCLUDES})
    target_link_libraries(${filename} libdxmc  ${XRAYLIB})
    add_test(${filename} ${filename})
    
endforeach()
