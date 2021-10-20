# Include directories in header search paths (-I flag)
INCLUDE_DIRECTORIES(../src/data_io)
INCLUDE_DIRECTORIES(../src/simulation_system)
INCLUDE_DIRECTORIES(../src/forces)
INCLUDE_DIRECTORIES(../src/integrators)
INCLUDE_DIRECTORIES(../src/helpers/eigen)

IF(DEFINED CMAKE_CUDA_COMPILER)
    INCLUDE_DIRECTORIES(../src/helpers/cuda)
ENDIF()

# Copy data input files for unit tests
CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/../test/input/collinear_swimmer_wall/initial_frame_dt1e-6_Z-height6.gsd"
    "${CMAKE_CURRENT_BINARY_DIR}/input/collinear_swimmer_wall/initial_frame_dt1e-6_Z-height6.gsd"
    COPYONLY
    )

# Catch2 CMake integration: https://github.com/catchorg/Catch2/blob/devel/docs/cmake-integration.md

# Add test executable
ADD_EXECUTABLE(tests 
               test-main.cpp 
               test-simulation.cpp)

TARGET_LINK_LIBRARIES(
    tests
    PUBLIC
    Catch2::Catch2
    spdlog::spdlog_header_only
    ${Boost_LIBRARIES}
    ${MKL_LIBRARIES}
    data_io
    simulation_system
    forces
    integrators
    )

INCLUDE(CTest)
INCLUDE(Catch)

catch_discover_tests(tests)
