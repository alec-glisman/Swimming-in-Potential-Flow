# SECTION: User information
#   FLAGS ACCEPTED:
#       CMAKE_BUILD_TYPE ["Release", "Debug", "Profile"]
#       ENABLE_TESTING ["True", "False"]
#       ENABLE_COVERAGE ["True", "False"]
#   GENERATORS TESTED:
#       "Unix Makefiles"   
#   COMPILERS TESTED:
#       gcc
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
    set(CMAKE_SKIP_BUILD_RPATH FALSE)
    set(CMAKE_BUILD_RPATH "/opt/intel/oneapi/mkl/latest/lib")

ELSEIF("$ENV{USER}" STREQUAL "aglisman")  # Build on desktop PC

    MESSAGE(STATUS "${BoldGreen}" "Compiling on PC-Ubuntu" "${ColourReset}")

    SET(CMAKE_CUDA_COMPILER "/opt/nvidia/hpc_sdk/Linux_x86_64/21.9/cuda/bin/nvcc")

    set(CMAKE_C_COMPILER /usr/bin/gcc-11 CACHE PATH "")
    set(CMAKE_CXX_COMPILER /usr/bin/g++-11 CACHE PATH "")

    
ELSEIF("$ENV{USER}" STREQUAL "alec")  # Build on Brady group servers

    MESSAGE(STATUS "${BoldGreen}" "Compiling on Brady group servers" "${ColourReset}")      
    
    INCLUDE_DIRECTORIES("/home/alec/local-new/intel/compilers_and_libraries/linux/mkl/include")
    INCLUDE_DIRECTORIES("/home/alec/local-new/intel/compilers_and_libraries/linux")
    INCLUDE_DIRECTORIES("/usr/lib/x86_64-linux-gnu")


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
SET(CMAKE_MODULE_PATH 
    ${CMAKE_MODULE_PATH}
    ${CMAKE_CURRENT_SOURCE_DIR}/include/cmake_modules)
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



# GCC requirements
IF("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")

    IF(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 11)

        MESSAGE(FATAL_ERROR "GCC version must be at least 11! Using: ${CMAKE_CXX_COMPILER_VERSION}")

    ENDIF()

ELSE()

    MESSAGE(WARNING "You are using an unsupported compiler! Compilation has only been tested with GCC.")
    MESSAGE(WARNING "Compiler ID: ${CMAKE_CXX_COMPILER_ID}")

ENDIF()
# !SECTION



# SECTION: Compilation flags
SET(CMAKE_CXX_FLAGS_DEBUG
    "-g -O0 -fno-omit-frame-pointer -fopenmp -lboost_system -ltbb -lpthread -m64 -lmkl_intel_lp64 -lmkl_tbb_thread -lmkl_core -lm -ldl -Wno-unknown-pragmas -Wall -Wextra -Wlogical-op -Wredundant-decls -Wshadow -Woverloaded-virtual -Wuninitialized -Warray-bounds -Wcast-align -Wcast-qual -Wdisabled-optimization -Wdiv-by-zero -Wendif-labels -Wformat-extra-args -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wimport -Winit-self -Winline -Winvalid-pch -Wlogical-op -Wmissing-declarations -Wno-missing-format-attribute -Wmissing-include-dirs -Wmultichar -Wpacked -Wpointer-arith -Wreturn-type -Wsequence-point -Wsign-compare -Wstrict-aliasing -Wstrict-aliasing=2 -Wswitch -Wno-unused -Wno-sign-compare -Wno-cast-qual -Wno-redundant-decls -Wvariadic-macros -Wwrite-strings"
)
SET(CMAKE_C_FLAGS_DEBUG
    "-g -O0 -fno-omit-frame-pointer -fopenmp -lboost_system -ltbb -lpthread -m64 -lmkl_intel_lp64 -lmkl_tbb_thread -lmkl_core -lm -ldl -Wno-unknown-pragmas -Wall -Wextra -Wlogical-op -Wredundant-decls -Wshadow -Wuninitialized -Warray-bounds -Wcast-align -Wcast-qual -Wdisabled-optimization -Wdiv-by-zero -Wendif-labels -Wformat-extra-args -Wformat-nonliteral -Wformat-security -Wformat-y2k -Wimport -Winit-self -Winline -Winvalid-pch -Wlogical-op -Wmissing-declarations -Wno-missing-format-attribute -Wmissing-include-dirs -Wmultichar -Wpacked -Wpointer-arith -Wreturn-type -Wsequence-point -Wsign-compare -Wstrict-aliasing -Wstrict-aliasing=2 -Wswitch -Wno-unused -Wno-sign-compare -Wno-cast-qual -Wno-redundant-decls -Wvariadic-macros -Wwrite-strings"
)


SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -O3 -mfpmath=sse -funroll-loops -fopenmp -lboost_system -ltbb -lpthread -m64 -lmkl_intel_lp64 -lmkl_tbb_thread -lmkl_core -lm -ldl -Wno-deprecated-declarations")
SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=native -O3 -mfpmath=sse -funroll-loops -fopenmp -lboost_system -ltbb -lpthread -m64 -lmkl_intel_lp64 -lmkl_tbb_thread -lmkl_core -lm -ldl -Wno-deprecated-declarations")


IF(CMAKE_BUILD_TYPE STREQUAL "Profile")

    SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fno-omit-frame-pointer -fopenmp -march=native -O3 -mfpmath=sse -funroll-loops -lboost_system -ltbb -lpthread -m64 -lmkl_intel_lp64 -lmkl_tbb_thread -lmkl_core -lm -ldl -Wno-deprecated-declarations")
    SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -fno-omit-frame-pointer -fopenmp -march=native -O3 -mfpmath=sse -funroll-loops -lboost_system -ltbb -lpthread -m64 -lmkl_intel_lp64 -lmkl_tbb_thread -lmkl_core -lm -ldl -Wno-deprecated-declarations")

ENDIF()

IF(ENABLE_COVERAGE STREQUAL "True" OR CMAKE_BUILD_TYPE STREQUAL "Debug")

    ADD_COMPILE_OPTIONS("--coverage")
    ADD_LINK_OPTIONS("--coverage")
   
    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_RELEASE} -fprofile-arcs -ftest-coverage")
    SET(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_RELEASE} -fprofile-arcs -ftest-coverage")
    
    SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fprofile-arcs -ftest-coverage")
    SET(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -fprofile-arcs -ftest-coverage")

ENDIF()

# !SECTION



# SECTION: External libraries
MESSAGE("")

# BOOST libraries
MESSAGE( "${BoldMagenta}" "External Dependency: Boost" "${ColourReset}" )
SET(Boost_USE_STATIC_LIBS OFF)
SET(Boost_USE_DEBUG_LIBS OFF)
SET(Boost_USE_RELEASE_LIBS ON)
SET(Boost_USE_MULTITHREADED ON)
SET(Boost_USE_STATIC_RUNTIME OFF)

FIND_PACKAGE(Boost 1.71.0 REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
MESSAGE(STATUS "${Cyan}" "Found: ${Boost_FOUND}" "${ColourReset}")
MESSAGE(STATUS "${Cyan}" "Version: ${Boost_VERSION}" "${ColourReset}")
MESSAGE(STATUS "${Cyan}" "Include Dir: ${Boost_INCLUDE_DIR}" "${ColourReset}")
MESSAGE("")

# Optimization libraries
MESSAGE( "${BoldMagenta}" "External Dependency: Intel MKL" "${ColourReset}" )
SET(MKL_THREADING tbb_thread)
SET(MKL_INTERFACE_FULL intel_lp64)
FIND_PACKAGE(MKL REQUIRED)    
MESSAGE( STATUS "${Cyan}" "Found: ${MKL_FOUND}" "${ColourReset}")
MESSAGE( STATUS "${Cyan}" "Include Dir: ${MKL_INCLUDE_DIR}" "${ColourReset}")
MESSAGE( STATUS "${Cyan}" "Library Dir: ${MKL_LIBRARIES}" "${ColourReset}")
MESSAGE("")

MESSAGE( "${BoldMagenta}" "External Dependency: Intel TBB" "${ColourReset}" )
FIND_PACKAGE(TBB REQUIRED)    
MESSAGE( STATUS "${Cyan}" "Found: ${TBB_FOUND}" "${ColourReset}")
MESSAGE( STATUS "${Cyan}" "Include Dir: ${TBB_INCLUDE_DIRS}" "${ColourReset}")
MESSAGE( STATUS "${Cyan}" "Library Dir: ${TBB_LIBRARIES}" "${ColourReset}")
MESSAGE("")

# Linear algebra libraries
MESSAGE( "${BoldMagenta}" "External Dependency: Eigen3" "${ColourReset}" )
FIND_PACKAGE(Eigen3 3.4.0 REQUIRED)
MESSAGE( STATUS "${Cyan}" "Found: ${EIGEN3_FOUND}" "${ColourReset}")
MESSAGE( STATUS "${Cyan}" "Include Dir: ${EIGEN3_INCLUDE_DIR}" "${ColourReset}")
MESSAGE("")

# Logging library
MESSAGE( "${BoldMagenta}" "External Dependency: spdlog" "${ColourReset}" )
FIND_PACKAGE(Spdlog 1.9.1 REQUIRED)
MESSAGE( STATUS "${Cyan}" "Found: ${SPDLOG_FOUND}" )
MESSAGE( STATUS "${Cyan}" "Include Dir: ${SPDLOG_INCLUDE_DIRS}" "${ColourReset}")
MESSAGE( STATUS "${Cyan}" "Libraries: ${SPDLOG_LIBRARIES}" "${ColourReset}")
MESSAGE("")

# !SECTION



# SECTION: Internal libraries
ADD_SUBDIRECTORY(src)

IF(ENABLE_TESTING STREQUAL "True" OR CMAKE_BUILD_TYPE STREQUAL "Debug")

    ENABLE_TESTING()

    # Unit testing libraries
    MESSAGE( "${BoldMagenta}" "External Dependency: Catch2" "${ColourReset}" )
    FIND_PACKAGE(Catch2 2.13.6 REQUIRED)
    MESSAGE( STATUS "${Cyan}" "Found: ${Catch2_FOUND}" "${ColourReset}")
    MESSAGE( STATUS "${Cyan}" "Include Dir: ${Catch2_INCLUDE_DIRS}" "${ColourReset}")
    INCLUDE(CTest)
    INCLUDE(Catch)
    MESSAGE("")

    ADD_SUBDIRECTORY(tests)

ENDIF()
# !SECTION
