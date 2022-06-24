

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


add_executable(PMMA pmma.cpp)
target_link_libraries(PMMA PRIVATE libdxmc)

install(TARGETS PMMA
	RUNTIME DESTINATION examples)
