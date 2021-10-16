# Make all of source code a library
IF(DEFINED CMAKE_CUDA_COMPILER)

    # COMPUTE ARCHITECTURES: GeForce RTX 3080: Compute Capability 8.6 GeForce
    # GTX 1080 TI: Compute Capability 6.1
    SET_TARGET_PROPERTIES(helper_eigen PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    SET_PROPERTY(TARGET helper_eigen PROPERTY CUDA_ARCHITECTURES 86 61)

ENDIF()

ADD_LIBRARY(
    helper_eigen
    helper_eigenTensorConversion.hpp
    )

TARGET_LINK_LIBRARIES(
    helper_eigen
    PUBLIC
    ${MKL_LIBRARIES}
    Eigen3::Eigen
)
