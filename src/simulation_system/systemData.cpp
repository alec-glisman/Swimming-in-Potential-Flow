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
    spdlog::get(m_logName)->info("Output path: {0}", m_outputDir);

    // Parse GSD file
    std::shared_ptr<systemData> shd_ptr = std::make_shared<systemData>(*this);
    gsdReader                           = std::make_shared<GSDReader>(shd_ptr);

    // TODO: Input checking
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
        spdlog::get(m_logName)->error("m_return_val = {0}", m_return_val);
        // throw std::runtime_error("Error parsing GSD file");
    }
    if (m_return_bool == false)
    {
        spdlog::get(m_logName)->error("m_return_bool = {0}", m_return_bool);
        // throw std::runtime_error("Error parsing GSD file");
    }
}

void
systemData::resizeTensors()
{
    int len = m_num_dim * m_num_particles;

    m_positions     = std::make_shared<Eigen::VectorXd>(len);
    m_velocities    = std::make_shared<Eigen::VectorXd>(len);
    m_accelerations = std::make_shared<Eigen::VectorXd>(len);
}