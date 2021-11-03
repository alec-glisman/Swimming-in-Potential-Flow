# Executable variables
SET(EXE_NAME "bodies-in-potential-flow")

SET(EXE_FILES 
    main.cpp 
    )

SET(EXE_LINKS 
    spdlog::spdlog_header_only
    ${MKL_LIBRARIES}
    Eigen3::Eigen
    data_io
    simulation_system
    forces
    integrators
    )


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