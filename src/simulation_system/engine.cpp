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

    // Initialize forces
    spdlog::get(m_logName)->info("Initializing potential hydrodynamics");
    potHydro = std::make_shared<potentialHydrodynamics>(sys);

    // Initialize integrator
    spdlog::get(m_logName)->info("Initializing integrator");
    rk4Integrator = std::make_shared<rungeKutta4>(sys, *potHydro);

    // Initialize progressBar
    spdlog::get(m_logName)->info("Initializing progressBar");
    unsigned int num_step = ceil(system->t() / system->dt());
    unsigned int barWidth = 70;
    progressBar           = std::make_shared<ProgressBar>(num_step, barWidth);
}

engine::~engine()
{
    spdlog::get(m_logName)->info("Destructing engine");
}

void
engine::run()
{
    spdlog::get(m_logName)->info("Starting engine run");
    progressBar->display(); // display the bar

    int write_step   = (int)ceil(system->tf() / system->dt() / system->numStepsOutput());
    int display_step = (int)ceil(system->tf() / system->dt() * outputPercentile);

    while (system->t() <= system->tf())
    {
        // Integrate system forward in time
        integrate();

        // Update time for next step
        system->setT(system->t() + system->dt());
        system->setTimestep(system->timestep() + 1);
        ++(*progressBar);

        // Output data
        if ((system->timestep() % write_step == 0) || (system->t() >= system->tf()))
        {
            spdlog::get(m_logName)->info("Writing frame at t = {0}", system->t());
            system->gsdUtil()->writeFrame();
        }
        if (system->timestep() % display_step == 0)
        {
            progressBar->display(); // display the bar
        }
    }
}

void
engine::integrate()
{
    rk4Integrator->integrate();
}