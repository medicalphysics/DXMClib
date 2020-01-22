


#add_executable(test_loadxraylib loadxraylib.cpp )
#target_include_directories(test_loadxraylib PRIVATE ${XRAYLIB_INCLUDES})
#target_link_libraries(test_loadxraylib ${XRAYLIB_LIBRARY})
#add_test(test_loadxraylib test_loadxraylib)
#set_tests_properties(${test_loadxraylib}
#    PROPERTIES
#    PASS_REGULAR_EXPRESSION "Test passed")
#set_tests_properties(${test_loadxraylib}
#    PROPERTIES
#    FAIL_REGULAR_EXPRESSION "(Exception|Test failed)")
#set_tests_properties(${test_loadxraylib}
#    PROPERTIES
#    TIMEOUT 120)



add_executable(test_material testmaterial.cpp )
#target_include_directories(test_material PRIVATE "${dxmclib_SOURCE_DIR}/include")
#target_include_directories(test_material PRIVATE "${libdxmc_INCLUDE_DIRECTORIES}")
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
