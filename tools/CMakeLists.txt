add_library(PlateTectonicsHeightmapIO STATIC heightmap_io.cpp)
target_include_directories(PlateTectonicsHeightmapIO
    PUBLIC
        "${PROJECT_SOURCE_DIR}/src"
        "${CMAKE_CURRENT_SOURCE_DIR}"
)
target_link_libraries(PlateTectonicsHeightmapIO PUBLIC PlateTectonics PNG::PNG TIFF::TIFF)
plate_tectonics_apply_common_settings(PlateTectonicsHeightmapIO)
plate_tectonics_set_output_dirs(PlateTectonicsHeightmapIO)

if(PLATE_TECTONICS_BUILD_SIMULATION)
    add_executable(simulation simulation.cpp)
    target_compile_definitions(simulation PRIVATE PLATE_TECTONICS_ROOT="${PROJECT_SOURCE_DIR}")
    target_link_libraries(simulation PRIVATE PlateTectonicsHeightmapIO)
    plate_tectonics_apply_common_settings(simulation)
    plate_tectonics_set_output_dirs(simulation)
endif()
