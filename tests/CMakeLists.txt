include(FetchContent)

FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG v1.15.2
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

set(plate_tectonics_test_sources
    test_acceptance.cpp
    test_bounds.cpp
    test_heightmap.cpp
    test_mass.cpp
    test_movement.cpp
    test_plate.cpp
    test_portability.cpp
    test_randomness.cpp
    test_rectangle.cpp
    test_regression.cpp
    test_sqrdmd.cpp
    test_topography_codec.cpp
    test_topography_io.cpp
)

add_executable(PlateTectonicsTests ${plate_tectonics_test_sources})
target_link_libraries(PlateTectonicsTests
    PRIVATE
        PlateTectonics
        PlateTectonicsHeightmapIO
        gtest
        gtest_main
)
plate_tectonics_apply_common_settings(PlateTectonicsTests)
plate_tectonics_set_output_dirs(PlateTectonicsTests)

add_test(NAME PlateTectonicsTests COMMAND $<TARGET_FILE:PlateTectonicsTests>)
set_tests_properties(PlateTectonicsTests PROPERTIES WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
