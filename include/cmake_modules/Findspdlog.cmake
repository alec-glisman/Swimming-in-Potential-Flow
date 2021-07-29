# Find the spdlog include directory
# The following variables are set if spdlog is found.
#  spdlog_FOUND        - True when the spdlog include directory is found.
#  spdlog_INCLUDE_DIR  - The path to where the spdlog include files are.
# If spdlog is not found, spdlog_FOUND is set to false.

find_package(PkgConfig)

if(NOT EXISTS "${spdlog_INCLUDE_DIR}")

    find_path(spdlog_INCLUDE_DIR
        NAMES spdlog/spdlog.h
        DOC "spdlog library header files"
        HINTS "/usr/local/include"
              "/usr/local/Cellar/spdlog/1.8.5"  # MacOS
              "${CMAKE_SOURCE_DIR}/include/spdlog/include" 
    )

endif()

if(EXISTS "${spdlog_INCLUDE_DIR}")

    include(FindPackageHandleStandardArgs)
    mark_as_advanced(spdlog_INCLUDE_DIR)
  
endif()

if(EXISTS "${spdlog_INCLUDE_DIR}")

  set(spdlog_FOUND 1)

else()

  set(spdlog_FOUND 0)

endif()