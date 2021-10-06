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

    spdlog::get(m_logName)->info("Constructor complete");
    spdlog::get(m_logName)->flush();
}

systemData::~systemData()
{
    spdlog::get(m_logName)->info("systemData destructor called");
    gsd_close(m_handle.get());
    spdlog::get(m_logName)->flush();
}

void
systemData::parseGSD()
{
    spdlog::get(m_logName)->info("Starting to parse GSD file");

    // parse GSD file and load data into *this
    m_gsdUtil    = std::make_shared<GSDUtil>(shared_from_this());
    m_GSD_parsed = true;

    // verify data is not a-physical
    checkInput();
}

void
systemData::check_gsd_return()
{
    if (m_return_val != 0)
    {
        spdlog::get(m_logName)->error("m_return_val = {0}", m_return_val);
        spdlog::get(m_logName)->flush();
        throw std::runtime_error("Error parsing GSD file");
    }
    if (m_return_bool == false)
    {
        spdlog::get(m_logName)->error("m_return_bool = {0}", m_return_bool);
        spdlog::get(m_logName)->flush();
        throw std::runtime_error("Error parsing GSD file");
    }
}

void
systemData::checkInput()
{
    spdlog::get(m_logName)->info("Input checking assertions");

    assert(m_positions_particles.size() == 4 * m_num_particles &&
           "Particle orientation (unit quaternions) vector has incorrect length, not 4N.");
    assert(m_positions_particles.size() == 3 * m_num_particles &&
           "Particle position vector has incorrect length, not 3N.");
    assert(m_velocities_particles.size() == 3 * m_num_particles &&
           "Particle velocity vector has incorrect length, not 3N.");
    assert(m_accelerations_particles.size() == 3 * m_num_particles &&
           "Particle acceleration vector has incorrect length, not 3N.");

    assert(m_positions_bodies.size() == 7 * m_num_particles &&
           "Body position vector has incorrect length, not 7N.");
    assert(m_velocities_bodies.size() == 7 * m_num_particles &&
           "Body velocity vector has incorrect length, not 7N.");
    assert(m_accelerations_bodies.size() == 7 * m_num_particles &&
           "Body acceleration vector has incorrect length, not 7N.");

    assert(m_particle_type_id.size() == m_num_particles &&
           "Particle type ID vector has incorrect length, not N.");

    assert(m_num_particles > 0 && "Must have at least one particle to simulate.");
    assert(m_num_bodies > 0 && "Must have at least one body to simulate.");
    assert(m_num_spatial_dim == 3 &&
           "Simulation framwork currently requires 3 spatial dimensions.");

    assert(m_tf > 0. && "Must have positive integration time.");
    assert(m_dt > 0.0 && "Integration time step must be positive.");
    assert(m_tf > m_t && "Final time must be greater than current time.");
    assert(m_tau > 0.0 && "Characteristic time must be greater than zero.");

    assert(m_fluid_density >= 0.0 && "Fluid density must be non-negative.");
    assert(m_particle_density >= 0.0 && "Particle density must be non-negative");

    assert(m_wca_epsilon >= 0.0 && "WCA_epsilon must be non-negative.");
    assert(m_wca_sigma >= 0.0 && "WCA_sigma must be non-negative.");
}
