# SECTION: User information
#   FLAGS ACCEPTED:
#       CMAKE_BUILD_TYPE ["Release", "Debug", "Profile"]
#       ENABLE_TESTING ["True", "False"]
#       CMAKE_EXPORT_COMPILE_COMMANDS ["True"] (for coc.nvim)
#   GENERATORS TESTED:
#       "Unix Makefiles"   
#!SECTION



# SECTION: Startup Messages
MESSAGE("WELCOME TO BODIES-IN-POTENTIAL-FLOW")
MESSAGE("-- USER environment variable is set to: " $ENV{USER})
MESSAGE("-- HOSTNAME environment variable is set to: " $ENV{HOSTNAME})
# !SECTION



# SECTION: Compiler and language configuration
IF(APPLE)  # Build on laptop PC

    # Special considerations for MacOS Builds
    MESSAGE(STATUS "** Building on MacOS, do not use frameworks **")
    # Try to avoid frameworks
    SET(CMAKE_MACOSX_RPATH OFF)
    SET(CMAKE_FIND_FRAMEWORK NEVER)
    SET(CMAKE_OSX_SYSROOT "")
    # Compiler options
    OPTION(USE_CLANG "build application with clang" OFF) # OFF is the default

    # Use GNU-10 Compilers on MacOS
    MESSAGE(STATUS "** Compiling with gcc-11 on MacOS **")
    CMAKE_POLICY(SET CMP0048 NEW)
    SET(CMAKE_CXX_STANDARD 17)
    SET(CMAKE_CXX_COMPILER "/usr/local/bin/g++-11" CACHE PATH "" FORCE)
    SET(CMAKE_C_COMPILER "/usr/local/bin/gcc-11" CACHE PATH "" FORCE)

    ADD_DEFINITIONS(-DMACOS_COMPILE)


ELSEIF($ENV{USER} STREQUAL "aglisman")  # Build on desktop PC

    MESSAGE(STATUS "** Compiling with gcc-11 on PC-Ubuntu **")

    SET(CMAKE_CXX_STANDARD 17)
    SET(CMAKE_CUDA_STANDARD 14)
    
    ADD_DEFINITIONS(-DMACOS_COMPILE)

    SET(CMAKE_CXX_COMPILER
        "/usr/bin/g++-11"
        CACHE PATH "" FORCE)
    SET(CMAKE_C_COMPILER
        "/usr/bin/gcc-11"
        CACHE PATH "" FORCE)
        
    SET(CMAKE_CUDA_COMPILER "/opt/nvidia/hpc_sdk/Linux_x86_64/21.5/cuda/bin/nvcc")

    
ELSEIF($ENV{USER} STREQUAL "alec")  # Build on Brady group servers

    MESSAGE(STATUS "** Compiling on Brady group servers **")      
    
    SET(CMAKE_CXX_STANDARD 14)
    CMAKE_POLICY(SET CMP0048 OLD)

    INCLUDE_DIRECTORIES(
        "/home/alec/local-new/intel/compilers_and_libraries/linux/mkl/include")
    INCLUDE_DIRECTORIES(
        "/home/alec/local-new/intel/compilers_and_libraries/linux")
    INCLUDE_DIRECTORIES("/usr/lib/x86_64-linux-gnu")


ELSE()

    MESSAGE(STATUS "** Compiling on unknown machine **")
    
    SET(CMAKE_CXX_STANDARD 17)
    SET(CMAKE_CUDA_STANDARD 14)
    
    SET(CMAKE_CXX_COMPILER
        "/usr/bin/g++-11"
        CACHE PATH "" FORCE)
    SET(CMAKE_C_COMPILER
        "/usr/bin/gcc-11"
        CACHE PATH "" FORCE)

ENDIF()



# Mandate given C++/CUDA requirements
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CUDA_STANDARD_REQUIRED ON)
# COMPUTE ARCHITECTURES: 
#     GeForce RTX 3080: Compute Capability 8.6 GeForce GTX,
#     1080 Ti: Compute Capability 6.1
SET(CMAKE_CUDA_ARCHITECTURES 86 61)



# CMake requirements
CMAKE_MINIMUM_REQUIRED(VERSION 3.5.1 FATAL_ERROR)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
                      ${CMAKE_CURRENT_SOURCE_DIR}/include/cmake_modules)

# Do not use <PackageName>_ROOT variables set in rc files
CMAKE_POLICY(SET CMP0074 NEW)
#!SECTION



# SECTION: Build configuration
IF(CMAKE_BUILD_TYPE STREQUAL "Debug")

    MESSAGE(STATUS "** Debug mode **")
    # SET(CMAKE_VERBOSE_MAKEFILE ON)

ELSEIF(CMAKE_BUILD_TYPE STREQUAL "Release")

    MESSAGE(STATUS "** Release mode **")

ELSEIF(CMAKE_BUILD_TYPE STREQUAL "Profile")

    MESSAGE(STATUS "** Profile mode **")

ENDIF()

MESSAGE("")
# !SECTION



# SECTION: Project configuration
# Add project command
IF(APPLE)

    PROJECT(
        bodies_in_potential_flow
        VERSION "0.1.0"
        DESCRIPTION
            "Study dynamics of inertial active matter in a potential fluid (irrotational, incompressible)"
        HOMEPAGE_URL
            "https://github.com/alec-glisman/bodies-in-potential-flow"
        LANGUAGES CXX C)

ELSE()  # Use CUDA when not on MacOS

    PROJECT(
        bodies_in_potential_flow
        VERSION "0.1.0"
        DESCRIPTION
            "Study dynamics of inertial active matter in a potential fluid (irrotational, incompressible)"
        HOMEPAGE_URL
            "https://github.com/alec-glisman/bodies-in-potential-flow"
        LANGUAGES CXX C CUDA)

ENDIF()
# !SECTION



# SECTION: Compilation flags
SET(CMAKE_CXX_FLAGS_DEBUG
    "-g -O0 -fno-omit-frame-pointer  -Wno-unknown-pragmas -Wall -Wextra -Wlogical-op -Wredundant-decls -Wshadow -Woverloaded-virtual -Wuninitialized -Warray-bounds -Wcast-align -Wcast-qual -Wdisabled-optimization -Wdiv-by-zero -Wendif-labels -Wformat-extra-args -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wimport -Winit-self -Winline -Winvalid-pch -Wlogical-op -Wmissing-declarations -Wno-missing-format-attribute -Wmissing-include-dirs -Wmultichar -Wpacked -Wpointer-arith -Wreturn-type -Wsequence-point -Wsign-compare -Wstrict-aliasing -Wstrict-aliasing=2 -Wswitch  -Wno-unused -Wvariadic-macros -Wwrite-strings"
)
SET(CMAKE_C_FLAGS_DEBUG
    "-g -O0 -fno-omit-frame-pointer -Wno-unknown-pragmas -Wall -Wextra -Wlogical-op -Wredundant-decls -Wshadow -Wuninitialized -Warray-bounds -Wcast-align -Wcast-qual -Wdisabled-optimization -Wdiv-by-zero -Wendif-labels -Wformat-extra-args -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wimport -Winit-self -Winline -Winvalid-pch -Wlogical-op -Wmissing-declarations -Wno-missing-format-attribute -Wmissing-include-dirs -Wmultichar -Wpacked -Wpointer-arith -Wreturn-type -Wsequence-point -Wsign-compare -Wstrict-aliasing -Wstrict-aliasing=2 -Wswitch -Wno-unused -Wvariadic-macros -Wwrite-strings"
)


SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -O3 -mfpmath=sse -funroll-loops -lboost_system -Wno-deprecated-declarations")
SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=native -O3 -mfpmath=sse -funroll-loops -lboost_system -Wno-deprecated-declarations")


IF(CMAKE_BUILD_TYPE STREQUAL "Profile")

    SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-omit-frame-pointer -march=native -O3 -mfpmath=sse -funroll-loops -lboost_system -Wno-deprecated-declarations")
    SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -fno-omit-frame-pointer -march=native -O3 -mfpmath=sse -funroll-loops -lboost_system -Wno-deprecated-declarations")

ENDIF()
# !SECTION



# SECTION: External libraries
# BOOST libraries
SET(Boost_USE_STATIC_LIBS OFF)
SET(Boost_USE_DEBUG_LIBS OFF)
SET(Boost_USE_RELEASE_LIBS ON)
SET(Boost_USE_MULTITHREADED ON)
SET(Boost_USE_STATIC_RUNTIME OFF)

FIND_PACKAGE(Boost REQUIRED)
IF(Boost_FOUND)

    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
    LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
    MESSAGE( "** Boost **" )
    MESSAGE(STATUS "Found: ${Boost_FOUND}")
    MESSAGE(STATUS "Version: ${Boost_VERSION}")
    MESSAGE(STATUS "Include Dir: ${Boost_INCLUDE_DIR}")

ELSE()

    MESSAGE(FATAL_ERROR "Boost was not located")

ENDIF()
MESSAGE("")

# Optimization libraries
MESSAGE( "** MKL **" )
FIND_PACKAGE(MKL REQUIRED)
MESSAGE( STATUS "Found: ${MKL_FOUND}" )
MESSAGE( STATUS "Include Dir: ${MKL_INCLUDE_DIR}" )
MESSAGE("")

# Linear algebra libraries
MESSAGE( "** EIGEN3 **" )
FIND_PACKAGE(Eigen3 REQUIRED)
MESSAGE( STATUS "Found: ${EIGEN3_FOUND}" )
MESSAGE( STATUS "Include Dir: ${EIGEN3_INCLUDE_DIR}" )
MESSAGE("")

# Logging library
MESSAGE( "** SPDLOG **" )
FIND_PACKAGE(Spdlog REQUIRED)
MESSAGE( STATUS "Found: ${SPDLOG_FOUND}" )
MESSAGE( STATUS "Include Dir: ${SPDLOG_INCLUDE_DIRS}" )
MESSAGE( STATUS "Libraries: ${SPDLOG_LIBRARIES}" )
MESSAGE("")

# Unit testing libraries
MESSAGE( "** CATCH2 **" )
FIND_PACKAGE(Catch2 REQUIRED)
MESSAGE( STATUS "Found: ${Catch2_FOUND}" )
MESSAGE( STATUS "Include Dir: ${Catch2_INCLUDE_DIRS}" )
INCLUDE(CTest)
INCLUDE(Catch)
MESSAGE("")
# !SECTION



# SECTION: Internal libraries
# Add other CMake files
ADD_SUBDIRECTORY(src)

IF(ENABLE_TESTING STREQUAL "True" OR CMAKE_BUILD_TYPE STREQUAL "Debug")

    ADD_SUBDIRECTORY(tests)

ENDIF()
# !SECTION
