# Library variables
SET(LIB_NAME "simulation_system")

SET(LIB_FILES 
    systemData.cpp systemData.hpp 
    engine.cpp engine.hpp
    progressBar.hpp)

SET(LIB_LINKS 
    ${MKL_LIBRARIES} 
    Eigen3::Eigen
    data_io)

# COMPUTE ARCHITECTURES: GeForce RTX 3080: Compute Capability 8.6 GeForce
# GTX 1080 TI: Compute Capability 6.1
IF(DEFINED CMAKE_CUDA_COMPILER)

    SET_TARGET_PROPERTIES(${LIB_NAME} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    SET_PROPERTY(TARGET ${LIB_NAME} PROPERTY CUDA_ARCHITECTURES 86 61)

ENDIF()

# Make all of source code a library
ADD_LIBRARY(
    ${LIB_NAME}
    ${LIB_FILES}
    )

# Link other libraries 
TARGET_LINK_LIBRARIES(
    ${LIB_NAME}
    PUBLIC
    ${LIB_LINKS}
    )
