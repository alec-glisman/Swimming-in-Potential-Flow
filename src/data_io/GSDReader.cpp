//
// Created by Alec Glisman on 07/30/21
// Building off of GSDReader class in HOOMD-Blue
//

#include <GSDReader.hpp>

GSDReader::GSDReader(std::string inputGSDFile, std::string outputDir)
{
    // Copy strings
    m_inputGSDFile = inputGSDFile;
    m_outputDir    = outputDir;

    // Initialize logger
    m_logFile   = m_outputDir + "/logs/" + m_logName + "-log.txt";
    auto logger = spdlog::basic_logger_mt(m_logName, m_logFile);
    spdlog::get(m_logName)->info("Initializing GSD reader");

    // readHeader();
    // readParticles();
    // readTopology();
}

GSDReader::~GSDReader() = default;
