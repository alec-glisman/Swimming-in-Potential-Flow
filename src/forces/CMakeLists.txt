# Make all of source code a library
IF(DEFINED CMAKE_CUDA_COMPILER)

    ADD_LIBRARY(
        forces
        potentialHydrodynamics.cpp 
        potentialHydrodynamics.hpp
        )

    TARGET_LINK_LIBRARIES(
        forces
        PRIVATE
        spdlog::spdlog
        Eigen3::Eigen
    )

    # COMPUTE ARCHITECTURES: GeForce RTX 3080: Compute Capability 8.6 GeForce
    # GTX 1080 TI: Compute Capability 6.1
    SET_TARGET_PROPERTIES(forces PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    SET_PROPERTY(TARGET forces PROPERTY CUDA_ARCHITECTURES 86 61)

ELSE()

    ADD_LIBRARY(
        forces
        potentialHydrodynamics.cpp 
        potentialHydrodynamics.hpp
        )

    TARGET_LINK_LIBRARIES(
        forces
        PRIVATE
        spdlog::spdlog
        Eigen3::Eigen
    )

ENDIF()
