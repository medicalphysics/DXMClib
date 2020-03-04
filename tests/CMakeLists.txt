

set(TEST_SOURCES 
    "testvectormath.cpp"
    "testmaterial.cpp"
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

add_executable(test_material testmaterial.cpp )
target_link_libraries(test_material libdxmc)
add_test(test_material test_material)
set_tests_properties(${test_material}
    PROPERTIES
    PASS_REGULAR_EXPRESSION "Test passed")
set_tests_properties(${test_material}
    PROPERTIES
    FAIL_REGULAR_EXPRESSION "(Exception|Test failed)")
set_tests_properties(${test_material}
    PROPERTIES
    TIMEOUT 120)

add_executable(test_beamfilters testbeamfilters.cpp )
target_link_libraries(test_beamfilters libdxmc)
add_test(test_beamfilters test_beamfilters)
set_tests_properties(${test_beamfilters}
    PROPERTIES
    PASS_REGULAR_EXPRESSION "Test passed")
set_tests_properties(${test_beamfilters}
    PROPERTIES
    FAIL_REGULAR_EXPRESSION "(Exception|Test failed)")
set_tests_properties(${test_beamfilters}
    PROPERTIES
    TIMEOUT 120)

