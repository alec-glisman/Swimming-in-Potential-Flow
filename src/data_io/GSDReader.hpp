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
#include <gsd.h> // GSD File

/* Include all external project dependencies */
// Logging
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/spdlog.h>
// STL
#include <memory> // for std::unique_ptr and std::shared_ptr
#include <string> // std::string

class GSDReader
{
  public:
    //! Loads in the file and parses the data
    GSDReader(std::string inputGSDFile, std::string outputDir);

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

    // GSD
    std::string m_inputGSDFile;

    // logging
    std::string m_outputDir;
    std::string m_logFile;
    std::string m_logName{"GSDReader"};
};

#endif // BODIES_IN_POTENTIAL_FLOW_GSD_READER_H
