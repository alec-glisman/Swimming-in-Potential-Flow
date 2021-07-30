//
// Created by Alec Glisman on 07/30/21
//

#include <systemData.hpp>

systemData::systemData(std::string inputGSDFile, std::string outputDir)
    : m_inputGSDFile(std::move(inputGSDFile)), m_outputDir(std::move(outputDir))
{
    // Initialize logger
    m_logFile = m_outputDir + "/logs/systemData-log.txt";
    m_logger  = spdlog::basic_logger_mt("systemData", m_logFile);
    spdlog::get("systemData")->info("Initializing system data");

    // Load GSD data
    m_return_val = gsd_open(m_handle.get(), m_inputGSDFile.c_str(), GSD_OPEN_READWRITE);
}

systemData::~systemData() = default;