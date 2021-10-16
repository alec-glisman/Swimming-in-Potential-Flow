# Include directories in header search paths (-I flag)
INCLUDE_DIRECTORIES(data_io)
INCLUDE_DIRECTORIES(simulation_system)
INCLUDE_DIRECTORIES(forces)
INCLUDE_DIRECTORIES(integrators)
INCLUDE_DIRECTORIES(helpers/eigen)

IF(DEFINED CMAKE_CUDA_COMPILER)
    INCLUDE_DIRECTORIES(helpers/cuda)
ENDIF()

# Add subdirectories to the build (processes CMakeLists.txt in these dirs)
ADD_SUBDIRECTORY(data_io)
ADD_SUBDIRECTORY(simulation_system)
ADD_SUBDIRECTORY(forces)
ADD_SUBDIRECTORY(integrators)

# Main project executable
ADD_EXECUTABLE(bodies_in_potential_flow main.cpp)

TARGET_LINK_LIBRARIES(
    bodies_in_potential_flow
    PUBLIC
    spdlog::spdlog_header_only
    ${Boost_LIBRARIES}
    ${MKL_LIBRARIES}
    Eigen3::Eigen
    data_io
    simulation_system
    forces
    integrators
    )
