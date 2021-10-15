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
#include <memory>    // for std::unique_ptr and std::shared_ptr
#include <stdexcept> // std::errors
#include <string>    // std::string

/* Forward declarations */
class systemData;

/**
 * @class GSDUtil
 *
 * @brief Wrapper class for `gsd.h` to load data from input GSD file into `systemData` class
 *
 */
class GSDUtil
{
  public:
    /**
     * @brief Construct a new GSDUtil object
     *
     * @param sys `systemData` class to collect data from. Its attributes will be overwritten using setter functions.
     */
    explicit GSDUtil(std::shared_ptr<systemData> sys);

    /**
     * @brief Construct a new GSDUtil object
     *
     * @param sys `systemData` class to collect data from. Its attributes will be overwritten using setter functions.
     * @param frame specific GSD frame number to parse
     */
    explicit GSDUtil(std::shared_ptr<systemData> sys, uint64_t frame);

    /**
     * @brief Destroy the GSDUtil object
     *
     */
    ~GSDUtil();

    /**
     * @brief Removes all frames and data from GSD file. This returns GSD to a "new file" state.
     *
     * @warning Make sure you are positive that you want to call this function. It is dangerous, as it deletes all data
     * inside GSD file.
     */
    void
    truncateGSD();

    /**
     * @brief Appends frame to GSD file using data from `systemData` class
     *
     * @details Ends GSD frame after data is written
     *
     */
    void
    writeFrame();

  private:
    /**
     * @brief Helper function that checks data parsing from GSD is successful.
     *
     */
    void
    checkGSDReturn();

    /**
     * @brief Finds and then loads data from GSD
     *
     * @details Expected size is in units of bytes
     *     u_int8: 1
     *     float (np.single, np.float32): 4
     *     double (np.double, np.float64): 8
     *
     * @param data pointer to variable function will load data into
     * @param frame GSD frame number to load data from
     * @param name Path of data chunk to load data from in GSD schema
     * @param expected_size n number of bytes that data is.
     * @param cur_n number of rows in chunk
     * @return true Data chunk found and loaded
     * @return false Data chunk not found
     */
    bool
    readChunk(void* data, uint64_t frame, const char* name, size_t expected_size, unsigned int cur_n = 0);

    /**
     * @brief Read header information from GSD frame `m_frame`
     *
     */
    void
    readHeader();

    /**
     * @brief Read parameter information from GSD frame `m_frame`
     *
     */
    void
    readParameters();

    /**
     * @brief Read particle information from GSD frame `m_frame`
     *
     */
    void
    readParticles();

    /**
     * @brief Reads system specific variables from GSD frame `m_frame`
     *
     * @review_swimmer change variables loaded based on variables needed in `systemData` for specific systems.
     */
    void
    readSystemSpecifics();

    /**
     * @brief Writes header information to GSD
     *
     */
    void
    writeHeader();

    /**
     * @brief Writes parameter information to GSD
     *
     */
    void
    writeParameters();

    /**
     * @brief Writes particle information to GSD
     *
     */
    void
    writeParticles();

    // classes
    /// shared pointer reference to `systemData` class
    std::shared_ptr<systemData> m_system;

    // GSD
    /// GSD frame number
    uint64_t m_frame;

    // logging
    /// path of logfile for spdlog to write to
    std::string m_logFile;
    /// filename of logfile for spdlog to write to
    const std::string m_logName{"GSDUtil"};

    /* SECTION: getters/setters */
  public:
    uint64_t
    frame() const
    {
        return m_frame;
    }
    /* !SECTION */
};

#endif // BODIES_IN_POTENTIAL_FLOW_GSD_READER_H
