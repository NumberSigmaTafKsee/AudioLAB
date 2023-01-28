cmake_minimum_required(VERSION 3.6 FATAL_ERROR)

project(SoundEditor LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_definitions(
	-c
)

if ( CMAKE_BUILD_TYPE STREQUAL "" )
	# CMake defaults to leaving CMAKE_BUILD_TYPE empty. This screws up
	# differentiation between debug and release builds.
	set(
		CMAKE_BUILD_TYPE
			"Debug"
		CACHE STRING
			"Choose the type of build, options are: None (CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) \"Debug\" \"Release\" \"RelWithDebInfo\" \"MinSizeRel\"."
		FORCE
	)
endif ( CMAKE_BUILD_TYPE STREQUAL "" )

if ( CMAKE_BUILD_TYPE STREQUAL "Release" )
add_definitions( -O2 )
endif ( CMAKE_BUILD_TYPE STREQUAL "Release" )

# hack for buggy CMake's FindwxWidgets script
if (DEFINED ENV_WX_CONFIG)
    set (ENV{WX_CONFIG} ${ENV_WX_CONFIG})
endif()

set(wxWidgets_USE_STATIC 1)

find_package(wxWidgets REQUIRED)

find_package(SDL2 REQUIRED)

#set(SRCS SoundEditor.cpp "window/Window.cpp" "audio/WAVReader.cpp" "audio/Player.cpp")

file(GLOB SRCS "*.cpp" "./*/*.cpp")

include(${wxWidgets_USE_FILE})

include_directories(${SDL2_INCLUDE_DIRS})

if (APPLE)
    # create bundle on apple compiles
    add_executable(${PROJECT_NAME} MACOSX_BUNDLE ${SRCS} )
    # Set a custom plist file for the app bundle - needed for Mac OS Retina display
    set_target_properties(${PROJECT_NAME} PROPERTIES MACOSX_BUNDLE_INFO_PLIST ${CMAKE_CURRENT_SOURCE_DIR}/Info.plist)
else()
    # the WIN32 is needed for Windows in order for it to look for WinMain
    # instead of the main function. This is ignored on other systems,
    # so it works on all platforms
    add_executable(${PROJECT_NAME} WIN32 ${SRCS})
endif()

target_include_directories(${PROJECT_NAME} PUBLIC "../thirdparty/portaudio")

# Define the link libraries
target_link_libraries(${PROJECT_NAME} PUBLIC ${SDL2_LIBRARIES})
target_link_libraries(${PROJECT_NAME} PUBLIC ${wxWidgets_LIBRARIES})
target_link_libraries(${PROJECT_NAME} PUBLIC portaudio)

