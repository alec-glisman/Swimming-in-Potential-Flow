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
#include <systemData.hpp>

/* Include all external project dependencies */
// Logging
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/spdlog.h>
// STL
#include <memory> // for std::unique_ptr and std::shared_ptr
#include <string> // std::string

/* Forward declarations */
class systemData;

class GSDReader
{
  public:
    // Loads in the file and parses the data
    GSDReader(std::shared_ptr<systemData> sys);

    GSDReader(std::shared_ptr<systemData> sys, uint64_t frame);

    //! Destructor
    ~GSDReader();

  private:
    bool
    readChunk(void* data, uint64_t frame, const char* name, size_t expected_size,
              unsigned int cur_n);

    void
    readHeader();

    void
    readParticles();

    // classes
    std::shared_ptr<systemData> system;

    // GSD
    std::string m_inputGSDFile;
    uint64_t    m_frame; //!< Cached frame

    // logging
    std::string m_outputDir;
    std::string m_logFile;
    std::string m_logName{"GSDReader"};
};

#endif // BODIES_IN_POTENTIAL_FLOW_GSD_READER_H
