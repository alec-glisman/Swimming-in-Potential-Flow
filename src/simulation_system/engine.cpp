//
// Created by Alec Glisman on 07/31/21
//

#include <engine.hpp>

engine::engine(systemData& sys)
{
    // save classes
    system = &sys;

    // Initialize logger
    m_logFile   = system->outputDir() + "/logs/" + m_logName + "-log.txt";
    auto logger = spdlog::basic_logger_mt(m_logName, m_logFile);
    spdlog::get(m_logName)->info("Initializing engine");

    // Create integrator
    spdlog::get(m_logName)->info("Creating integrator");
    rk4Integrator = std::make_shared<rungeKutta4>(sys);
}

engine::~engine()
{
    spdlog::get(m_logName)->info("Destructing engine");
}

void
engine::run()
{
    spdlog::get(m_logName)->info("Starting engine run");
    int write_step = (int)ceil(system->tf() / system->dt() / system->numStepsOutput());

    while (system->t() <= system->tf())
    {
        // Integrate system forward in time
        integrate();

        // Update time for next step
        system->setT(system->t() + system->dt());
        system->setTimestep(system->timestep() + 1);

        // Output data
        if ((system->timestep() % write_step == 0) || (system->t() >= system->tf()))
        {
            spdlog::get(m_logName)->info("Writing frame at t = {0}", system->t());
            system->gsdUtil()->writeFrame();
        }
    }
}

void
engine::acceleration_update()
{
}

void
engine::integrate()
{
}