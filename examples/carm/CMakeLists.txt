
# testing to see if dxmclib target is defined. If it's not we will download 
# dxmclib at configure time.
IF( NOT TARGET libdxmc)
	include(FetchContent)
	## Adding DXMClib package
	FetchContent_Declare(
		libdxmc
		GIT_REPOSITORY https://github.com/medicalphysics/DXMClib.git
		GIT_TAG master
		)
	FetchContent_MakeAvailable(libdxmc)
ENDIF()


add_executable(carm carm.cpp phantomreader.cpp phantomreader.hpp)
target_link_libraries(carm PRIVATE libdxmc)

set(CARM_ARM_PATH "${CMAKE_CURRENT_SOURCE_DIR}/data/carmCarm.stl")
set(CARM_TABLE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/data/carmTable.stl")

add_custom_command(TARGET carm POST_BUILD
			COMMENT "Copying carmCarm.stl for C-arm example"
			DEPENDS ${CARM_ARM_PATH}
			COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CARM_ARM_PATH} carm.stl
			WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
			)
add_custom_command(TARGET carm POST_BUILD
			COMMENT "Copying carmTable.stl for C-arm example"
			DEPENDS ${CARM_TABLE_PATH}
			COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CARM_TABLE_PATH} table.stl
			WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
			)


dxmclib_add_physics_list(carm)

install(TARGETS carm
	RUNTIME DESTINATION examples)
