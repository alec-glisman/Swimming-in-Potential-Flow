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
CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/input/data.gsd"
    "${CMAKE_CURRENT_BINARY_DIR}/input/data.gsd"
    COPYONLY
    )

# Add catch as an interface library
ADD_LIBRARY(Catch INTERFACE)

# Add test executable
ADD_EXECUTABLE(tests 
               testsMain.cpp 
               testsSimulationBuild.cpp)

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

ADD_TEST(NAME tests COMMAND tests)

ENABLE_TESTING()
