cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(frankentorch)

find_package(Torch REQUIRED)
add_compile_options(-fmax-errors=1)
add_executable(frankentorch frankentorch.cpp)
include_directories("/home/james5/Audio/AudioLAB/include;/home/james5/Audio/AudioLAB/StdSamples")
target_link_libraries(frankentorch "${TORCH_LIBRARIES}")
set_property(TARGET frankentorch PROPERTY CXX_STANDARD 17)

