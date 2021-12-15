# SECTION: User information
#   FLAGS ACCEPTED:
#       CMAKE_BUILD_TYPE ["Release", "Debug", "Profile"]
#       ENABLE_TESTING ["True", "False"]
#       ENABLE_COVERAGE ["True", "False"]
#   GENERATORS TESTED:
#       "Unix Makefiles"   
#       "Ninja"
#   COMPILERS TESTED:
#       gcc-11
#!SECTION

# SECTION: Color output variables
if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(ColourBold  "${Esc}[1m")
  set(Red         "${Esc}[31m")
  set(Green       "${Esc}[32m")
  set(Yellow      "${Esc}[33m")
  set(Blue        "${Esc}[34m")
  set(Magenta     "${Esc}[35m")
  set(Cyan        "${Esc}[36m")
  set(White       "${Esc}[37m")
  set(BoldRed     "${Esc}[1;31m")
  set(BoldGreen   "${Esc}[1;32m")
  set(BoldYellow  "${Esc}[1;33m")
  set(BoldBlue    "${Esc}[1;34m")
  set(BoldMagenta "${Esc}[1;35m")
  set(BoldCyan    "${Esc}[1;36m")
  set(BoldWhite   "${Esc}[1;37m")
endif()
# !SECTION


# SECTION: Startup Messages
MESSAGE("${BoldYellow}" "WELCOME TO BODIES-IN-POTENTIAL-FLOW" "${ColourReset}")
MESSAGE("${BoldGreen}" "USER environment variable is set to: " "$ENV{USER}" "${ColourReset}")
MESSAGE("${BoldGreen}" "HOSTNAME environment variable is set to: " "$ENV{HOSTNAME}" "${ColourReset}" "\n")
# !SECTION



# SECTION: Compiler and language configuration
IF(APPLE)  # Build on laptop PC    

    MESSAGE(STATUS "${BoldGreen}" "Compiling on MacOS" "${ColourReset}")

    SET(CMAKE_MACOSX_RPATH 1)
    SET(CMAKE_BUILD_RPATH  
        ${CMAKE_BUILD_RPATH}
        "/opt/intel/oneapi/mkl/latest/lib"
        "/opt/intel/oneapi/tbb/latest/lib")
    
    SET(CMAKE_C_COMPILER /usr/local/bin/gcc-11 CACHE PATH "")
    SET(CMAKE_CXX_COMPILER /usr/local/bin/g++-11 CACHE PATH "")
    
    SET(GCOV_PATH "/usr/local/bin/gcov-11")

ELSEIF("$ENV{USER}" STREQUAL "aglisman")  # Build on desktop PC

    MESSAGE(STATUS "${BoldGreen}" "Compiling on PC-Ubuntu" "${ColourReset}")

    SET(CMAKE_CUDA_COMPILER "/opt/nvidia/hpc_sdk/Linux_x86_64/21.9/cuda/bin/nvcc")

    SET(CMAKE_C_COMPILER /usr/bin/gcc-11 CACHE PATH "")
    SET(CMAKE_CXX_COMPILER /usr/bin/g++-11 CACHE PATH "")
    
    SET(GCOV_PATH "/usr/bin/gcov-11")

    
ELSEIF("$ENV{USER}" STREQUAL "alec")  # Build on Brady group servers

    MESSAGE(STATUS "${BoldGreen}" "Compiling on Brady group servers" "${ColourReset}")      
    
    SET(CMAKE_C_COMPILER /usr/bin/gcc-11 CACHE PATH "")
    SET(CMAKE_CXX_COMPILER /usr/bin/g++-11 CACHE PATH "")


ELSE()

    MESSAGE(STATUS "${BoldGreen}" "Compiling on unknown machine" "${ColourReset}")
    
ENDIF()



# Export compile commands
SET(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# C++/CUDA requirements
SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_CUDA_STANDARD 14)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)
SET(CMAKE_CUDA_STANDARD_REQUIRED ON)

# CUDA COMPUTE ARCHITECTURES: 
#     GeForce RTX 3080:    Compute Capability 8.6
#     GeForce GTX 1080 Ti: Compute Capability 6.1
SET(CMAKE_CUDA_ARCHITECTURES 86 61)

# CMake requirements
CMAKE_MINIMUM_REQUIRED(VERSION 3.5.1 FATAL_ERROR)
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/include/cmake_modules")

# vcpkg integrartion
if(DEFINED ENV{VCPKG_ROOT} AND NOT DEFINED CMAKE_TOOLCHAIN_FILE)
  set(CMAKE_TOOLCHAIN_FILE "$ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake"
      CACHE STRING "")
endif()
#!SECTION



# SECTION: Build configuration
IF(CMAKE_BUILD_TYPE STREQUAL "Debug")
    MESSAGE(STATUS "${BoldWhite}" "Debug mode" "${ColourReset}")
    # SET(CMAKE_VERBOSE_MAKEFILE ON)

ELSEIF(CMAKE_BUILD_TYPE STREQUAL "Release")
    MESSAGE(STATUS "${BoldWhite}" "Release mode" "${ColourReset}")

ELSEIF(CMAKE_BUILD_TYPE STREQUAL "Profile")
    MESSAGE(STATUS "${BoldWhite}" "Profile mode" "${ColourReset}")

ENDIF()

MESSAGE("")
# !SECTION



# SECTION: Project configuration
# Add project command
IF(APPLE)
    SET(LANGS CXX C)
    
ELSE()  # Use CUDA when not on MacOS
    SET(LANGS CXX C)
    # SET(LANGS CXX CUDA)  # TODO: (future) Add cuda back later (removed for GitHub actions)

ENDIF()


PROJECT(
        bodies-in-potential-flow
        VERSION "0.1.0"
        DESCRIPTION
            "Study dynamics of inertial active matter in a potential fluid (irrotational, incompressible)"
        LANGUAGES ${LANGS})


# GCC requirements
IF("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    IF(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11)
        MESSAGE(WARNING "GCC version must be at least 11! Using: ${CMAKE_CXX_COMPILER_VERSION}")
    ENDIF()

ELSE()
    MESSAGE(WARNING "You are using an unsupported compiler! Compilation has only been tested with GCC.")
    MESSAGE(WARNING "Compiler ID: ${CMAKE_CXX_COMPILER_ID}")

ENDIF()

# !SECTION



# SECTION: Compilation flags

# Debug build flags
STRING(CONCAT CXX_DEBUG "-g -O0 -fno-omit-frame-pointer -fopenmp"
    " -Wno-unknown-pragmas -Wall -Wextra -Wlogical-op -Wredundant-decls -Wshadow -Woverloaded-virtual -Wuninitialized"
    " -Warray-bounds -Wcast-align -Wcast-qual -Wdisabled-optimization -Wdiv-by-zero -Wendif-labels -Wformat-extra-args"
    " -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wimport -Winit-self -Winline -Winvalid-pch -Wlogical-op"
    " -Wmissing-declarations -Wno-missing-format-attribute -Wmissing-include-dirs -Wmultichar -Wpacked -Wpointer-arith"
    " -Wreturn-type -Wsequence-point -Wsign-compare -Wstrict-aliasing -Wstrict-aliasing=2 -Wswitch -Wno-unused"
    " -Wno-sign-compare -Wno-cast-qual -Wno-redundant-decls -Wvariadic-macros -Wwrite-strings")
STRING(CONCAT C_DEBUG "-g -O0 -fno-omit-frame-pointer -fopenmp"
    " -Wno-unknown-pragmas -Wall -Wextra -Wlogical-op -Wredundant-decls -Wshadow -Wuninitialized"
    " -Warray-bounds -Wcast-align -Wcast-qual -Wdisabled-optimization -Wdiv-by-zero -Wendif-labels -Wformat-extra-args"
    " -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wimport -Winit-self -Winline -Winvalid-pch -Wlogical-op" 
    " -Wmissing-declarations -Wno-missing-format-attribute -Wmissing-include-dirs -Wmultichar -Wpacked"
    " -Wpointer-arith -Wreturn-type -Wsequence-point -Wsign-compare -Wstrict-aliasing -Wstrict-aliasing=2 -Wswitch -Wno-unused"
    " -Wno-sign-compare -Wno-cast-qual -Wno-redundant-decls -Wvariadic-macros -Wwrite-strings")

SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${CXX_DEBUG}")
SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${C_DEBUG}")

# Release build flags
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -O3 -mfpmath=sse -funroll-loops -fopenmp  -Wno-deprecated-declarations")
SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=native -O3 -mfpmath=sse -funroll-loops -fopenmp -Wno-deprecated-declarations")

# Profile build flags
IF(CMAKE_BUILD_TYPE STREQUAL "Profile")

    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O0 -fno-omit-frame-pointer")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -fno-omit-frame-pointer")

    SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -g -fno-omit-frame-pointer -march=native -mfpmath=sse -funroll-loops -Wno-deprecated-declarations")
    SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -fno-omit-frame-pointer -march=native -mfpmath=sse -funroll-loops -Wno-deprecated-declarations")

ENDIF()
# !SECTION



# SECTION: External libraries
MESSAGE("")

# Logging library
MESSAGE( "${BoldMagenta}" "External Dependency: spdlog" "${ColourReset}" )
FIND_PACKAGE(Spdlog 1.9.1 REQUIRED)
MESSAGE( STATUS "${Cyan}" "Found: ${SPDLOG_FOUND}" )
MESSAGE( STATUS "${Cyan}" "Include Dir: ${SPDLOG_INCLUDE_DIRS}" "${ColourReset}")
MESSAGE( STATUS "${Cyan}" "Libraries: ${SPDLOG_LIBRARIES}" "${ColourReset}")
MESSAGE("")

# Optimization libraries
MESSAGE( "${BoldMagenta}" "External Dependency: Intel TBB" "${ColourReset}" )
FIND_PACKAGE(TBB)    
MESSAGE( STATUS "${Cyan}" "Found: ${TBB_FOUND}" "${ColourReset}")
MESSAGE( STATUS "${Cyan}" "Include Dir: ${TBB_INCLUDE_DIRS}" "${ColourReset}")
MESSAGE( STATUS "${Cyan}" "Library Dir: ${TBB_LIBRARIES}" "${ColourReset}")
MESSAGE("")

MESSAGE( "${BoldMagenta}" "External Dependency: Intel MKL" "${ColourReset}" )
SET(MKL_INTERFACE_FULL intel_lp64)
IF(${TBB_FOUND})
    SET(MKL_THREADING tbb_thread)
ELSE()
    SET(MKL_THREADING gnu_thread)
ENDIF()

FIND_PACKAGE(MKL)    
MESSAGE( STATUS "${Cyan}" "Found: ${MKL_FOUND}" "${ColourReset}")
MESSAGE( STATUS "${Cyan}" "Include Dir: ${MKL_INCLUDE_DIR}" "${ColourReset}")
MESSAGE( STATUS "${Cyan}" "Library Dir: ${MKL_LIBRARIES}" "${ColourReset}")
MESSAGE("")

# MKL and TBB specific flags
IF(${MKL_FOUND} AND ${TBB_FOUND})
    
    ADD_COMPILE_OPTIONS(-ltbb -lpthread -m64 -lmkl_intel_lp64 -lmkl_tbb_thread -lmkl_core -lm -ldl)

ELSEIF(${MKL_FOUND})
    
    ADD_COMPILE_OPTIONS(-m64 -lmkl_intel_lp64 -lmkl_core -lm -ldl)

ENDIF()


# Linear algebra libraries
MESSAGE( "${BoldMagenta}" "External Dependency: Eigen3" "${ColourReset}" )
FIND_PACKAGE(Eigen3 3.3.9 REQUIRED)
MESSAGE( STATUS "${Cyan}" "Found: ${EIGEN3_FOUND}" "${ColourReset}")
MESSAGE( STATUS "${Cyan}" "Include Dir: ${EIGEN3_INCLUDE_DIR}" "${ColourReset}")
MESSAGE("")

IF(ENABLE_COVERAGE STREQUAL "True")

    MESSAGE( "${BoldMagenta}" "External Dependency: Code coverage" "${ColourReset}" )
    INCLUDE(CodeCoverage)
    APPEND_COVERAGE_COMPILER_FLAGS()
    MESSAGE("")

    # optional excludes - None needed here
    SET(COVERAGE_EXCLUDES
        "/home/aglisman/.linuxbrew/**"
        "/usr/include/**"
        "include/**"
        "/usr/lib/**"
        "/opt/intel/**"
        "/Applications/**"
        "${PROJECT_SOURCE_DIR}/src/main.cpp"
        "test/**"
        "vcpkg/**"
        "vcpkg_installed/**"
        "build/**"
        )

ENDIF()


IF(ENABLE_TESTING STREQUAL "True")

    # Unit testing libraries
    MESSAGE( "${BoldMagenta}" "External Dependency: Catch2" "${ColourReset}" )
    FIND_PACKAGE(Catch2 2.13.6 REQUIRED)
    MESSAGE( STATUS "${Cyan}" "Found: ${Catch2_FOUND}" "${ColourReset}")
    MESSAGE( STATUS "${Cyan}" "Include Dir: ${Catch2_INCLUDE_DIRS}" "${ColourReset}")
    MESSAGE("")

    ENABLE_TESTING()

ENDIF()

# !SECTION


# SECTION: Internal libraries
# ADD_SUBDIRECTORY(src)

# IF(ENABLE_TESTING STREQUAL "True")

#     ADD_SUBDIRECTORY(test)

# ENDIF()
# !SECTION
