//
// Created by Alec Glisman on 07/31/21
//

#include <engine.hpp>

engine::engine(std::shared_ptr<systemData> sys)
{
    // save classes
    m_system = sys;

    // Initialize logger
    m_logFile   = m_system->outputDir() + "/logs/" + m_logName + "-log.txt";
    auto logger = spdlog::basic_logger_mt(m_logName, m_logFile);
    spdlog::get(m_logName)->info("Initializing engine");

    // validate system loaded GSD data
    spdlog::get(m_logName)->info("Checking input systemData class loaded GSD data: {0}",
                                 m_system->gSDParsed());
    assert(m_system->gSDParsed() == true &&
           "GSD data not loaded before calling engine constructor");

    // Initialize forces
    spdlog::get(m_logName)->info("Initializing potential hydrodynamics");
    m_potHydro = std::make_shared<potentialHydrodynamics>(m_system);

    // Initialize integrator
    spdlog::get(m_logName)->info("Initializing integrator");
    m_rk4Integrator = std::make_shared<rungeKutta4>(m_system, m_potHydro);

    // Initialize progressBar
    spdlog::get(m_logName)->info("Initializing progressBar");
    int num_step = (int)ceil(m_system->tf() / m_system->dt());
    spdlog::get(m_logName)->info("Numer of integration steps: {0}", num_step);
    unsigned int barWidth = 70;
    m_progressBar = std::make_shared<ProgressBar>(static_cast<unsigned int>(num_step), barWidth);

    spdlog::get(m_logName)->info("Constructor complete");
    spdlog::get(m_logName)->flush();
}

engine::~engine()
{
    spdlog::get(m_logName)->info("Destructing engine");
    spdlog::get(m_logName)->flush();
}

void
engine::run()
{
    spdlog::get(m_logName)->critical("Starting engine run");
    m_progressBar->display(); // display the bar

    int tot_step     = (int)ceil(m_system->tf() / m_system->dt());
    int write_step   = (int)ceil(m_system->tf() / m_system->dt() / m_system->numStepsOutput());
    int display_step = (int)ceil(m_system->tf() / m_system->dt() * m_outputPercentile);

    spdlog::get(m_logName)->info("Write step: {0}", write_step);
    spdlog::get(m_logName)->info("Display step: {0}", display_step);

    while (m_system->timestep() < tot_step)
    {
        // Integrate system forward in time
        integrate();

        // Update time for next step
        m_system->setT(m_system->t() + m_system->dt());
        m_system->setTimestep(m_system->timestep() + 1);
        ++(*m_progressBar);

        // Output data
        if ((m_system->timestep() % write_step == 0) || (m_system->t() >= m_system->tf()))
        {
            spdlog::get(m_logName)->info("Writing frame at t = {0}", m_system->t());
            m_system->gsdUtil()->writeFrame();
        }
        if (m_system->timestep() % display_step == 0)
        {
            m_progressBar->display(); // display the bar
        }
        spdlog::get(m_logName)->flush();
    }

    // Final data writing and shut down
    spdlog::get(m_logName)->info("Ending engine run");
    spdlog::get(m_logName)->info("Writing frame at t = {0}", m_system->t());
    m_system->gsdUtil()->writeFrame();
    m_progressBar->display();
    spdlog::get(m_logName)->flush();
}

void
engine::integrate()
{
    m_rk4Integrator->integrate();
}