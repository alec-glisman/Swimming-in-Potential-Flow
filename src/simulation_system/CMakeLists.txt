# Make all of source code a library
IF(DEFINED CMAKE_CUDA_COMPILER)

    ADD_LIBRARY(
        simulation_system
        systemData.cpp 
        systemData.hpp 
        engine.cpp
        engine.hpp
    )
    
    TARGET_LINK_LIBRARIES(
        simulation_system
        PRIVATE
        Eigen3::Eigen
        spdlog::spdlog
        data_io
    )
    
    # COMPUTE ARCHITECTURES: GeForce RTX 3080: Compute Capability 8.6 GeForce
    # GTX 1080 TI: Compute Capability 6.1
    SET_TARGET_PROPERTIES(simulation_system PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    SET_PROPERTY(TARGET simulation_system PROPERTY CUDA_ARCHITECTURES 86 61)

ELSE()

    ADD_LIBRARY(
        simulation_system
        systemData.cpp 
        systemData.hpp 
        engine.cpp 
        engine.hpp
    )
        
    TARGET_LINK_LIBRARIES(
        simulation_system
        PRIVATE
        Eigen3::Eigen
        spdlog::spdlog
        data_io
    )

ENDIF()
