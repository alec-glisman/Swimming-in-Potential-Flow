//
// Created by Alec Glisman on 07/30/21
//

#include <systemData.hpp>

systemData::systemData(std::string inputGSDFile, std::string outputDir)
    : m_inputGSDFile(std::move(inputGSDFile)), m_outputDir(std::move(outputDir))
{
    // Initialize logger
    // auto r = spdlog::basic_logger_mt("file_logger", m_outputDir + "logs/basic-log.txt");

    // Load GSD data
    m_return_val = gsd_open(m_handle.get(), m_inputGSDFile.c_str(), GSD_OPEN_READWRITE);
}

systemData::~systemData() = default;