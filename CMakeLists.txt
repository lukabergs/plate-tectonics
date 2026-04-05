cmake_minimum_required(VERSION 3.21)

project(PlateTectonics LANGUAGES CXX)

include(CTest)

option(PLATE_TECTONICS_BUILD_SIMULATION "Build the simulation tool." ON)
option(PLATE_TECTONICS_WARNINGS_AS_ERRORS "Treat warnings as errors." ON)
option(PLATE_TECTONICS_ENABLE_ASAN "Enable AddressSanitizer." OFF)
option(PLATE_TECTONICS_ENABLE_UBSAN "Enable UndefinedBehaviorSanitizer." OFF)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(MSVC)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
endif()

function(plate_tectonics_apply_common_settings target)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4 /permissive- /EHsc)
        if(PLATE_TECTONICS_WARNINGS_AS_ERRORS)
            target_compile_options(${target} PRIVATE /WX)
        endif()
    else()
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
        if(PLATE_TECTONICS_WARNINGS_AS_ERRORS)
            target_compile_options(${target} PRIVATE -Werror)
        endif()
    endif()

    if(PLATE_TECTONICS_ENABLE_ASAN)
        if(MSVC)
            target_compile_options(${target} PRIVATE /fsanitize=address)
            target_link_options(${target} PRIVATE /fsanitize=address)
        else()
            target_compile_options(${target} PRIVATE -fsanitize=address -fno-omit-frame-pointer)
            target_link_options(${target} PRIVATE -fsanitize=address)
        endif()
    endif()

    if(PLATE_TECTONICS_ENABLE_UBSAN AND NOT MSVC)
        target_compile_options(${target} PRIVATE -fsanitize=undefined -fno-omit-frame-pointer)
        target_link_options(${target} PRIVATE -fsanitize=undefined)
    endif()
endfunction()

function(plate_tectonics_set_output_dirs target)
    set_target_properties(${target} PROPERTIES
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
        PDB_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
        COMPILE_PDB_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
    )

    foreach(config Debug Release RelWithDebInfo MinSizeRel)
        string(TOUPPER "${config}" config_upper)
        set_target_properties(${target} PROPERTIES
            ARCHIVE_OUTPUT_DIRECTORY_${config_upper} "${CMAKE_BINARY_DIR}/lib"
            LIBRARY_OUTPUT_DIRECTORY_${config_upper} "${CMAKE_BINARY_DIR}/bin"
            RUNTIME_OUTPUT_DIRECTORY_${config_upper} "${CMAKE_BINARY_DIR}/bin"
            PDB_OUTPUT_DIRECTORY_${config_upper} "${CMAKE_BINARY_DIR}/bin"
            COMPILE_PDB_OUTPUT_DIRECTORY_${config_upper} "${CMAKE_BINARY_DIR}/bin"
        )
    endforeach()
endfunction()

set(plate_tectonics_library_sources
    src/bounds.cpp
    src/geometry.cpp
    src/heightmap.cpp
    src/lithosphere.cpp
    src/mass.cpp
    src/movement.cpp
    src/noise.cpp
    src/plate.cpp
    src/plate_functions.cpp
    src/platecapi.cpp
    src/rectangle.cpp
    src/segment_creator.cpp
    src/segment_data.cpp
    src/segments.cpp
    src/simplerandom.cpp
    src/simplexnoise.cpp
    src/sqrdmd.cpp
    src/topography_codec.cpp
    src/utils.cpp
    src/world_point.cpp
)

add_library(PlateTectonics STATIC ${plate_tectonics_library_sources})
target_include_directories(PlateTectonics PUBLIC "${PROJECT_SOURCE_DIR}/src")
target_compile_features(PlateTectonics PUBLIC cxx_std_20)
plate_tectonics_apply_common_settings(PlateTectonics)
plate_tectonics_set_output_dirs(PlateTectonics)

if(BUILD_TESTING OR PLATE_TECTONICS_BUILD_SIMULATION)
    find_package(PNG REQUIRED)
    add_subdirectory(tools)
endif()

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
