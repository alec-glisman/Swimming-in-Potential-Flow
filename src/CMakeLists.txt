# Include directories in header search paths (-I flag)
INCLUDE_DIRECTORIES(data_io)
INCLUDE_DIRECTORIES(simulation_system)
INCLUDE_DIRECTORIES(forces)
INCLUDE_DIRECTORIES(integrators)

IF(DEFINED CMAKE_CUDA_COMPILER)
    INCLUDE_DIRECTORIES(cuda_helpers)
ENDIF()

# Add subdirectories to the build (processes CMakeLists.txt in these dirs)
ADD_SUBDIRECTORY(data_io)
ADD_SUBDIRECTORY(simulation_system)
ADD_SUBDIRECTORY(forces)
ADD_SUBDIRECTORY(integrators)

# Main project executable
ADD_EXECUTABLE(bodies_in_potential_flow main.cpp)

IF(CMAKE_BUILD_TYPE STREQUAL "Debug")

    TARGET_LINK_LIBRARIES(
        bodies_in_potential_flow
        PRIVATE
        ${Boost_LIBRARIES}
        Catch2::Catch2
        spdlog::spdlog
        data_io
        simulation_system
        forces
        integrators
        )

ELSEIF(CMAKE_BUILD_TYPE STREQUAL "Release")

    TARGET_LINK_LIBRARIES(
        bodies_in_potential_flow
        PRIVATE
        ${Boost_LIBRARIES}
        Catch2::Catch2
        spdlog::spdlog
        data_io
        simulation_system
        forces
        integrators
        )

ENDIF()
