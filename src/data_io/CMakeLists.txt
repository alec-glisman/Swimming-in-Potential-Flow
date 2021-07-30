# Make all of source code a library
IF(DEFINED CMAKE_CUDA_COMPILER)

    ADD_LIBRARY(
        data_io
        gsd.c 
        gsd.h)

    # COMPUTE ARCHITECTURES: GeForce RTX 3080: Compute Capability 8.6 GeForce
    # GTX 1080 TI: Compute Capability 6.1
    SET_TARGET_PROPERTIES(forces PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    SET_PROPERTY(TARGET forces PROPERTY CUDA_ARCHITECTURES 86 61)

ELSE()

    ADD_LIBRARY(
        data_io
        gsd.c 
        gsd.h)

ENDIF()
