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
    m_gsdUtil = std::make_shared<GSDUtil>(*this);

    // Input checking
    assert(m_fluid_density >= 0.0 && "Fluid density must be non-negative.");
    assert(m_particle_density >= 0.0 && "Particle density must be non-negative");
    assert(m_num_particles > 0 && "Must have at least one particle to simulate.");
    assert(m_tf > 0. && "Must have positive integration time.");
    assert(m_dt > 0.0 && "Integration time step must be positive");
    assert(m_wca_epsilon >= 0.0 && "WCA_epsilon must be non-negative");
    assert(m_wca_sigma >= 0.0 && "WCA_sigma must be non-negative");

    assert((*m_positions).size() == m_num_dim * m_num_particles &&
           "Position vector has incorrect length, not 3N");
    assert((*m_velocities).size() == m_num_dim * m_num_particles &&
           "Velocity vector has incorrect length, not 3N");
    assert((*m_accelerations).size() == m_num_dim * m_num_particles &&
           "Acceleration vector has incorrect length, not 3N");
}

systemData::~systemData()
{
    spdlog::get(m_logName)->info("systemData desstructor called");
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