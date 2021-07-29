# Include directories in header search paths (-I flag)
# INCLUDE_DIRECTORIES(forces)
# INCLUDE_DIRECTORIES(input)
# INCLUDE_DIRECTORIES(initial_configurations)
# INCLUDE_DIRECTORIES(integrators)
# INCLUDE_DIRECTORIES(simulation_system)
# INCLUDE_DIRECTORIES(utils)

IF(DEFINED CMAKE_CUDA_COMPILER)
    INCLUDE_DIRECTORIES(cuda_helpers)
ENDIF()

# Add subdirectories to the build (processes CMakeLists.txt in these dirs)
# ADD_SUBDIRECTORY(forces)
# ADD_SUBDIRECTORY(input)
# ADD_SUBDIRECTORY(initial_configurations)
# ADD_SUBDIRECTORY(integrators)
# ADD_SUBDIRECTORY(simulation_system)

# Main project executable
ADD_EXECUTABLE(bodies_in_potential_flow main.cpp)

IF(CMAKE_BUILD_TYPE STREQUAL "Debug")

    TARGET_LINK_LIBRARIES(
        bodies_in_potential_flow
        PRIVATE ${Boost_LIBRARIES}
        PRIVATE Catch2::Catch2
        PRIVATE Eigen3::Eigen
        )

ELSEIF(CMAKE_BUILD_TYPE STREQUAL "Release")

    TARGET_LINK_LIBRARIES(
        bodies_in_potential_flow
        PRIVATE ${Boost_LIBRARIES}
        PRIVATE Catch2::Catch2
        PRIVATE Eigen3::Eigen
        )

ENDIF()
