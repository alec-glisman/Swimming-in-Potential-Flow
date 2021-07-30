# Include directories in header search paths (-I flag)
INCLUDE_DIRECTORIES(../src/data_io)
INCLUDE_DIRECTORIES(../src/simulation_system)

IF(DEFINED CMAKE_CUDA_COMPILER)
    INCLUDE_DIRECTORIES(../src/cuda_helpers)
ENDIF()

# Add catch as an interface library
ADD_LIBRARY(Catch INTERFACE)

# Add test executable
ADD_EXECUTABLE(tests testsMain.cpp testsSimulationBuild.cpp)

IF(CMAKE_BUILD_TYPE STREQUAL "Debug")

    TARGET_LINK_LIBRARIES(
        tests
        PRIVATE
        Eigen3::Eigen
        ${Boost_LIBRARIES}
        Catch2::Catch2
        spdlog::spdlog
        data_io
        simulation_system
        )

ELSEIF(CMAKE_BUILD_TYPE STREQUAL "Release")

    TARGET_LINK_LIBRARIES(
        tests
        PRIVATE
        Eigen3::Eigen
        ${Boost_LIBRARIES}
        Catch2::Catch2
        spdlog::spdlog
        data_io
        simulation_system
        )

ENDIF()

ADD_TEST(NAME tests COMMAND tests)

ENABLE_TESTING()
