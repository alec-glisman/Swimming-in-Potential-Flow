//
// Created by Alec Glisman on 07/30/21
// Building off of GSDReader class in HOOMD-Blue
//

#ifndef BODIES_IN_POTENTIAL_FLOW_GSD_READER_H
#define BODIES_IN_POTENTIAL_FLOW_GSD_READER_H

#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif

/* Include all internal project dependencies */
#include <gsd.h>          // GSD File
#include <systemData.hpp> // simulation data

/* Include all external project dependencies */
// Logging
#include <spdlog/sinks/basic_file_sink.h>
// STL
#include <memory> // for std::unique_ptr and std::shared_ptr
#include <string> // std::string

/* Forward declarations */
class systemData;

class GSDReader
{
  public:
    //! Loads in the file and parses the data
    GSDReader(std::shared_ptr<systemData> sys);

    //! Destructor
    ~GSDReader();

  private:
    // helper functions to read sections of the file
    void
    readHeader();

    void
    readParticles();

    void
    readTopology();

    // classes
    std::shared_ptr<systemData> system;

    // logging
    std::string                     m_logFile;
    std::string                     m_logName{"GSDReader"};
    std::shared_ptr<spdlog::logger> m_logger;
};

#endif // BODIES_IN_POTENTIAL_FLOW_GSD_READER_H
