

add_executable(validation validation.cpp)
target_link_libraries(validation libdxmc)




set(DATA_PATH "${CMAKE_CURRENT_SOURCE_DIR}/data")
set(CASE5WORLD_BINARY "case5world.bin")
set(CASE5WORLD_COMPRESSED "case5world.tar.gz")

if(EXISTS "${DATA_PATH}/${CASE5WORLD_BINARY}")
	# If this binary file is in sources we compress it
	execute_process(COMMAND
	  ${CMAKE_COMMAND} -E tar -czvf ${CASE5WORLD_COMPRESSED} ${CASE5WORLD_BINARY}
	  WORKING_DIRECTORY ${DATA_PATH}
	  )
	# Then we delete it
	file(REMOVE "${DATA_PATH}/${CASE5WORLD_BINARY}")
endif()


add_custom_command(TARGET validation POST_BUILD
			COMMENT "Extracting case 5 world ${CASE5WORLDFILE} to ${CMAKE_CURRENT_BINARY_DIR}/data/case5world.bin"
			DEPENDS ${CASE5WORLDFILE}
			COMMAND ${CMAKE_COMMAND} -E tar -xvzf "${DATA_PATH}/${CASE5WORLD_COMPRESSED}" ${CASE5WORLD_BINARY}
			WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
			)

set(PLOTSCRIPT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/validationPlots.py")
set(PLOTSCRIPT_DEST "${CMAKE_CURRENT_BINARY_DIR}/validationPlots.py")
add_custom_command(TARGET validation POST_BUILD
	COMMENT "Copy plotting script ${PLOTSCRIPT_PATH} to  ${PLOTSCRIPT_DEST}"
	DEPENDS ${PLOTSCRIPT_PATH}
    COMMAND ${CMAKE_COMMAND} -E copy ${PLOTSCRIPT_PATH} ${PLOTSCRIPT_DEST}
	)
