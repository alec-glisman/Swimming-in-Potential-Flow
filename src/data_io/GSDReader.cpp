//
// Created by Alec Glisman on 07/30/21
// Building off of GSDReader class in HOOMD-Blue
//

#include <GSDReader.hpp>

GSDReader::GSDReader(std::shared_ptr<systemData> sys)
{
    // Save classes
    system = sys;

    // Initialize logger
    m_logFile = system->outputDir() + "/logs/" + m_logName + "-log.txt";
    m_logger  = spdlog::basic_logger_mt(m_logName, m_logFile);
    spdlog::get(m_logName)->info("Initializing GSD reader");

    readHeader();
    readParticles();
    readTopology();
}

GSDReader::~GSDReader() = default;
