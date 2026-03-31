cmake_minimum_required (VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project (PlateTectonicsExamples)

find_package(PNG REQUIRED)

add_executable(simulation simulation.cpp map_drawing.cpp)

target_include_directories(simulation PRIVATE ../src)
target_link_libraries(simulation PRIVATE PlateTectonics PNG::PNG)
