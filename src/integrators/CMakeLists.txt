# Make all of source code a library
IF(DEFINED CMAKE_CUDA_COMPILER)

    # COMPUTE ARCHITECTURES: GeForce RTX 3080: Compute Capability 8.6 GeForce
    # GTX 1080 TI: Compute Capability 6.1
    SET_TARGET_PROPERTIES(integrators PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    SET_PROPERTY(TARGET integrators PROPERTY CUDA_ARCHITECTURES 86 61)

ENDIF()

ADD_LIBRARY(
    integrators
    rungeKutta4.cpp 
    rungeKutta4.hpp
    )

TARGET_LINK_LIBRARIES(
    integrators
    PUBLIC
    spdlog::spdlog_header_only
    ${MKL_LIBRARIES}
    Eigen3::Eigen
)
