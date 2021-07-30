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
#include <memory> // for std::unique_ptr and std::shared_ptr
#include <string> // std::string

class GSDReader
{
  public:
    //! Loads in the file and parses the data
    GSDReader(const std::string& name, const uint64_t frame);

    //! Destructor
    ~GSDReader();

  private:
    std::string m_name;   //!< Cached file name
    uint64_t    m_frame;  //!< Cached frame
    gsd_handle  m_handle; //!< Handle to the file

    // helper functions to read sections of the file
    void
    readHeader();

    void
    readParticles();

    void
    readTopology();
};

#endif // BODIES_IN_POTENTIAL_FLOW_GSD_READER_H
