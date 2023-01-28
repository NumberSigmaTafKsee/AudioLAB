cmake_minimum_required(VERSION 3.20)
project(wavestudio)

set(CMAKE_CXX_STANDARD 17)
if (${CMAKE_HOST_SYSTEM_NAME} STREQUAL "Linux")
    set(IMGUI_PATH /home/ben/Documents/repos/imgui)
    #set(IMGUIFILEDIALOG_PATH "/home/ben/Documents/repos/ImGuiFileDialog")
    set(FILEBROWSER_ADDON_PATH "/home/ben/Documents/repos/ImGui-Addons/FileBrowser")
    set(GL3W_PATH ${IMGUI_PATH}/examples/libs/gl3w)
    set(BACKENDS_PATH ${IMGUI_PATH}/backends)
    set(LIBSNDFILE_PATH /home/ben/Downloads/libsndfile-1.0.31)
    include_directories(/usr/include/SDL2/ ${GL3W_PATH}
            ${IMGUI_PATH}
            ${IMGUI_PATH}/backends
            ${FILEBROWSER_ADDON_PATH}
            ${LIBSNDFILE}/include
            )
    link_directories(${LIBSNDFILE_PATH}/src/.libs)
    add_compile_definitions(IMGUI_IMPL_OPENGL_LOADER_GL3W _REENTRANT)
endif ()

add_executable(wavestudio
        main.cpp
        ${GL3W_PATH}/GL/gl3w.c
        ${IMGUI_PATH}/backends/imgui_impl_opengl3.cpp
        ${IMGUI_PATH}/backends/imgui_impl_sdl.cpp
        ${IMGUI_PATH}/imgui.cpp
        ${IMGUI_PATH}/imgui_demo.cpp
        ${IMGUI_PATH}/imgui_draw.cpp
        ${IMGUI_PATH}/imgui_tables.cpp
        ${IMGUI_PATH}/imgui_widgets.cpp
        ${FILEBROWSER_ADDON_PATH}/ImGuiFileBrowser.cpp
        wavestudio_application.cpp wavestudio_application.h)
target_link_libraries(wavestudio SDL2 GL dl)

