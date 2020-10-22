
set(TEST_SOURCES 
    "testtube.cpp"
    "testsource.cpp"
    "testmaterial.cpp"
    "testrandom.cpp"
    "testattenuationlut.cpp"
    "testworld.cpp"
    "testvectormath.cpp"
)

foreach(sourcefile ${TEST_SOURCES})
    get_filename_component(filename ${sourcefile} NAME_WLE)
    add_executable(${filename} ${sourcefile})
    target_link_libraries(${filename} libdxmc)
    add_test(${filename} ${filename})
    set_tests_properties(${filename}
        PROPERTIES
        PASS_REGULAR_EXPRESSION "Test passed")
    set_tests_properties(${filename}
        PROPERTIES
        FAIL_REGULAR_EXPRESSION "(Exception|Test failed)")
    set_tests_properties(${filename}
        PROPERTIES
        TIMEOUT 120)
endforeach()
