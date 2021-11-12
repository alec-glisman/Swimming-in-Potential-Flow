//
// Created by Alec Glisman on 07/31/21
//

#include <Engine.hpp>

Engine::Engine(std::shared_ptr<SystemData> sys)
{
    // save classes
    m_system = sys;

    // Initialize logger
    m_logFile   = m_system->outputDir() + "/logs/" + m_logName + "-log.txt";
    auto logger = spdlog::basic_logger_mt(m_logName, m_logFile);
    spdlog::get(m_logName)->info("Initializing Engine");

    // validate system loaded GSD data
    spdlog::get(m_logName)->info("Checking input SystemData class loaded GSD data: {0}", m_system->gSDParsed());
    if (m_system->gSDParsed() == false)
    {
        throw std::runtime_error("GSD data not loaded into SystemData class before calling Engine constructor.");
    }

    // Initialize forces
    spdlog::get(m_logName)->info("Initializing potential hydrodynamics");
    m_potHydro = std::make_shared<PotentialHydrodynamics>(m_system);

    // Initialize integrator
    spdlog::get(m_logName)->info("Initializing integrator");
    m_rk4Integrator = std::make_shared<RungeKutta4>(m_system, m_potHydro);

    // Initialize ProgressBar
    spdlog::get(m_logName)->info("Initializing ProgressBar");
    int num_step = (int)ceil(m_system->tf() / m_system->dt());
    spdlog::get(m_logName)->info("Numer of integration steps: {0}", num_step);
    unsigned int barWidth = 70;
    m_ProgressBar         = std::make_shared<ProgressBar>(static_cast<unsigned int>(num_step), barWidth);

    spdlog::get(m_logName)->info("Constructor complete");
    spdlog::get(m_logName)->flush();
}

Engine::~Engine()
{
    spdlog::get(m_logName)->info("Destructing Engine");
    spdlog::get(m_logName)->flush();
    spdlog::drop(m_logName);
}

void
Engine::run()
{
    spdlog::get(m_logName)->critical("Starting Engine run");
    m_ProgressBar->display(); // display the progress bar

    // calculate number of steps in integration
    double t_remaining{m_system->tf() - m_system->t()};

    int tot_step     = (int)ceil(t_remaining / m_system->dt());
    int write_step   = (int)ceil(t_remaining / m_system->dt() / m_system->numStepsOutput());
    int display_step = (int)ceil(t_remaining / m_system->dt() * m_outputPercentile);

    spdlog::get(m_logName)->info("Write step: {0}", write_step);
    spdlog::get(m_logName)->info("Display step: {0}", display_step);

    // create eigen3 thread-pool and device
    Eigen::ThreadPool       thread_pool = Eigen::ThreadPool(m_simulation_cores);
    Eigen::ThreadPoolDevice all_cores_device(&thread_pool, m_simulation_cores);

    while (m_system->timestep() < tot_step)
    {
        // Integrate system forward in time
        integrate(all_cores_device);

        // Update time for next step
        m_system->setT(m_system->t() + m_system->dt());
        m_system->setTimestep(m_system->timestep() + 1);
        ++(*m_ProgressBar);

        // Output data
        if ((m_system->timestep() % write_step == 0) || (m_system->t() >= m_system->tf()))
        {
            spdlog::get(m_logName)->info("Normalizing quaternions at t = {0}", m_system->t());
            m_system->normalizeQuaternions();

            spdlog::get(m_logName)->info("Writing frame at t = {0}", m_system->t());
            m_system->gsdUtil()->writeFrame();
            m_system->logData();
            spdlog::get(m_logName)->flush();
        }
        if (m_system->timestep() % display_step == 0)
        {
            m_ProgressBar->display(); // display the progress bar
        }
    }

    // Final data writing and shut down
    spdlog::get(m_logName)->info("Ending Engine run");
    spdlog::get(m_logName)->info("Writing frame at t = {0}", m_system->t());
    m_system->gsdUtil()->writeFrame();
    m_ProgressBar->display();
    spdlog::get(m_logName)->flush();
}

void
Engine::integrate(Eigen::ThreadPoolDevice& device)
{
    m_rk4Integrator->integrate(device);
}
