//
// Created by Alec Glisman on 07/30/21
//

#ifndef BODIES_IN_POTENTIAL_FLOW_SYSTEM_DATA_H
#define BODIES_IN_POTENTIAL_FLOW_SYSTEM_DATA_H

#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif

/* Include all internal project dependencies */
#include "gsd.h" // GSD File

/* Include all external project dependencies */
// Logging
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/spdlog.h>

// STL
#include <memory> // for std::unique_ptr and std::shared_ptr
#include <string> // std::string

class systemData
{
  public:
    systemData(std::string inputGSDFile, std::string outputDir);

    ~systemData();

  private:
    std::string m_inputGSDFile;
    std::string m_outputDir;
    std::string m_logFile;

    std::shared_ptr<spdlog::logger> m_logger;

    std::shared_ptr<gsd_handle> m_handle{new gsd_handle};
    int                         m_return_val{-1};
};

#endif
