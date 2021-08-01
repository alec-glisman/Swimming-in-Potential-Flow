//
// Created by Alec Glisman on 07/31/21
//

#include <potentialHydrodynamics.hpp>

potentialHydrodynamics::potentialHydrodynamics(systemData& sys)
{
    // save classes
    system = &sys;

    // Initialize logger
    m_logFile   = system->outputDir() + "/logs/" + m_logName + "-log.txt";
    auto logger = spdlog::basic_logger_mt(m_logName, m_logFile);
    spdlog::get(m_logName)->info("Initializing potential hydrodynamics");
}

potentialHydrodynamics::~potentialHydrodynamics()
{
    spdlog::get(m_logName)->info("Destructing potential hydrodynamics");
}