# Library variables
SET(LIB_NAME "data_io")

SET(LIB_FILES 
    gsd.c gsd.h
    GSDUtil.cpp GSDUtil.hpp)

SET(LIB_LINKS 
    spdlog::spdlog_header_only 
    ${MKL_LIBRARIES} 
    Eigen3::Eigen)

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

# COMPUTE ARCHITECTURES: GeForce RTX 3080: Compute Capability 8.6 GeForce
# GTX 1080 TI: Compute Capability 6.1
IF(DEFINED CMAKE_CUDA_COMPILER)

    SET_TARGET_PROPERTIES(${LIB_NAME} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    SET_PROPERTY(TARGET ${LIB_NAME} PROPERTY CUDA_ARCHITECTURES 86 61)

ENDIF()
