//
// Created by Alec Glisman on 07/30/21
//

#include <systemData.hpp>

systemData::systemData(std::string inputGSDFile, std::string outputDir)
    : m_inputGSDFile(std::move(inputGSDFile)), m_outputDir(std::move(outputDir))
{
    // Initialize logger
    m_logFile   = m_outputDir + "/logs/" + m_logName + "-log.txt";
    auto logger = spdlog::basic_logger_mt(m_logName, m_logFile);
    spdlog::get(m_logName)->info("Initializing system data");

    // Load GSD frame
    spdlog::get(m_logName)->info("Loading input GSD file");
    m_return_val = gsd_open(m_handle.get(), m_inputGSDFile.c_str(), GSD_OPEN_READWRITE);
    check_gsd_return();

    // Begin to parse GSD file
    gsdReader = std::make_shared<GSDReader>(m_inputGSDFile, m_outputDir);

    // Parameters: degrees of freedom
    spdlog::get(m_logName)->info("Parameters: degrees of freedom");
}

systemData::~systemData()
{
    gsd_close(m_handle.get());
}

void
systemData::check_gsd_return()
{
    if (m_return_val != 0)
    {
        spdlog::get(m_logName)->critical("m_return_val = {0}", m_return_val);
    }
}
