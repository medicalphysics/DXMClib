

add_executable(pencilbeam pencilbeam.cpp)
target_link_libraries(pencilbeam PRIVATE libdxmc)

install(TARGETS pencilbeam
	RUNTIME DESTINATION examples)
