# Executable variables
SET(EXE_NAME "tests")

SET(EXE_FILES 
    test-main.cpp 
    test-simulation.cpp
    )

SET(EXE_LINKS 
    Catch2::Catch2
    spdlog::spdlog_header_only
    ${Boost_LIBRARIES}
    ${MKL_LIBRARIES}
    data_io
    simulation_system
    forces
    integrators
    test_simulation_system
    )

# Copy data input files for unit tests
CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/../test/input/collinear_swimmer_wall/initial_frame_dt1e-6_Z-height6.gsd"
    "${CMAKE_CURRENT_BINARY_DIR}/input/collinear_swimmer_wall/initial_frame_dt1e-6_Z-height6.gsd"
    COPYONLY
    )

CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/../test/input/collinear_swimmer_wall/initial_frame_dt1e-1_Z-height6.gsd"
    "${CMAKE_CURRENT_BINARY_DIR}/input/collinear_swimmer_wall/initial_frame_dt1e-1_Z-height6.gsd"
    COPYONLY
    )


# Include source directories in header search paths (-I flag)
INCLUDE_DIRECTORIES(../src/data_io)
INCLUDE_DIRECTORIES(../src/simulation_system)
INCLUDE_DIRECTORIES(../src/forces)
INCLUDE_DIRECTORIES(../src/integrators)
INCLUDE_DIRECTORIES(../src/helpers/eigen)

IF(DEFINED CMAKE_CUDA_COMPILER)
    INCLUDE_DIRECTORIES(../src/helpers/cuda)
ENDIF()

# Include test directories in header search paths (-I flag)
INCLUDE_DIRECTORIES(simulation_system)


# Add subdirectories to the build (processes CMakeLists.txt in these dirs)
ADD_SUBDIRECTORY(simulation_system)


# Make all of source code a library
ADD_EXECUTABLE(
    ${EXE_NAME}
    ${EXE_FILES}
    )

# Link other libraries 
TARGET_LINK_LIBRARIES(
    ${EXE_NAME}
    PUBLIC
    ${EXE_LINKS}
    )


# Catch2 CMake integration: https://github.com/catchorg/Catch2/blob/devel/docs/cmake-integration.md
INCLUDE(CTest)
INCLUDE(Catch)
CATCH_DISCOVER_TESTS(tests)

IF(ENABLE_COVERAGE STREQUAL "True")

    # lcov: https://stackoverflow.com/questions/52255589/code-coverage-analysis-using-codecoverage-cmake-with-multiple-targets
    SETUP_TARGET_FOR_COVERAGE_LCOV(
            NAME coverage                 # New target name
            EXECUTABLE ctest -j          # Executable in PROJECT_BINARY_DIR
            DEPENDENCIES ${EXE_NAME}     # Dependencies to build first
        )

    # gcovr: https://stackoverflow.com/a/62730462/13215572
    SETUP_TARGET_FOR_COVERAGE_GCOVR_HTML(
                NAME ${EXE_NAME}_coverage_html
                EXECUTABLE ${EXE_NAME}
                DEPENDENCIES ${EXE_LINKS}
        )

ENDIF()