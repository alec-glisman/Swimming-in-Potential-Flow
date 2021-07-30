# Make all of source code a library
IF(DEFINED CMAKE_CUDA_COMPILER)

    ADD_LIBRARY(
        data_io
        gsd.c 
        gsd.h
        GSDReader.cpp
        GSDReader.hpp
        )

    TARGET_LINK_LIBRARIES(
        data_io
        PRIVATE
        spdlog::spdlog
    )

    # COMPUTE ARCHITECTURES: GeForce RTX 3080: Compute Capability 8.6 GeForce
    # GTX 1080 TI: Compute Capability 6.1
    SET_TARGET_PROPERTIES(data_io PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    SET_PROPERTY(TARGET data_io PROPERTY CUDA_ARCHITECTURES 86 61)

ELSE()

    ADD_LIBRARY(
        data_io
        gsd.c 
        GSDReader.cpp
        GSDReader.hpp
        )

    TARGET_LINK_LIBRARIES(
        data_io
        PRIVATE
        spdlog::spdlog
    )

ENDIF()
