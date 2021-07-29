# Find the spdlog include directory
# The following variables are set if spdlog is found.
#  CATCH_FOUND        - True when the spdlog include directory is found.
#  CATCH2_INCLUDE_DIR  - The path to where the spdlog include files are.
# If spdlog is not found, spdlog_FOUND is set to false.

find_package(PkgConfig)

if(NOT EXISTS "${CATCH2_INCLUDE_DIR}")

    SET(Catch2_CHECK_PATHS 
        "/usr/local/include"
        "/usr/local/homebrew/include" # Mac OS X
        "/opt/local/var/macports/software" # Mac OS X.
        "/opt/local/include"
        "/usr/include"
        "$ENV{CPLUS_INCLUDE_PATH}"
        "$ENV{CPATH}"
        "${CMAKE_SOURCE_DIR}/thirdparty/catch"
        )

    find_path(CATCH2_INCLUDE_DIR
        NAMES 
        DOC "spdlog library header files"
        HINTS "${CMAKE_SOURCE_DIR}/include/Catch2/include" 
        )

endif()

if(EXISTS "${CATCH2_INCLUDE_DIR}")

    include(FindPackageHandleStandardArgs)
    mark_as_advanced(CATCH2_INCLUDE_DIR)
  
endif()

if(EXISTS "${CATCH2_INCLUDE_DIR}")

  set(CATCH2_FOUND 1)

else()

  set(CATCH2_FOUND 0)

endif()