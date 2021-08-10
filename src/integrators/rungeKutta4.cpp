//
// Created by Alec Glisman on 07/31/21
//

#include <rungeKutta4.hpp>

rungeKutta4::rungeKutta4(systemData& sys, potentialHydrodynamics& hydro)
{
    // save classes
    system   = &sys;
    potHydro = &hydro;

    // Initialize logger
    m_logFile   = system->outputDir() + "/logs/" + m_logName + "-log.txt";
    auto logger = spdlog::basic_logger_mt(m_logName, m_logFile);
    spdlog::get(m_logName)->info("Initializing Runge-Kutta 4th order integrator");

    // Create integrator
    spdlog::get(m_logName)->info("Creating integrator");
}

rungeKutta4::~rungeKutta4()
{
    spdlog::get(m_logName)->info("Destructing Runge-Kutta 4th order");
}

void
rungeKutta4::integrate()
{
    /* Time unit conversion */
    double t  = system->timestep() * system->tau();
    double dt = system->dt() * system->tau();

    /* Step 1: k1 = f(t, h) */
    Eigen::VectorXd x1 = system->positions;
    Eigen::VectorXd v1 = system->velocities;
    Eigen::VectorXd a1 = system->accelerations;

    /* Step 2: k2 = f(t + h/2, y + h*k1/2) */
    Eigen::VectorXd x2          = x1 + (0.5 * dt) * v1;
    system->positions.noalias() = x2;
    potHydro->update();
    Eigen::VectorXd v2 = v1 + (0.5 * dt) * a1;
    Eigen::VectorXd a2;
    acceleration_update();

    /* Step 3: k3 = f(t + h/2, y + h/2*k2) */
    Eigen::VectorXd x3          = x1 + (0.5 * dt) * v2;
    system->positions.noalias() = x3;
    potHydro->update();
    Eigen::VectorXd v3 = v1 + (0.5 * dt) * a2;
    Eigen::VectorXd a3;
    acceleration_update();

    /* Step 4: k4 = f(t + h, y + h*k3) */
    Eigen::VectorXd x4          = x1 + dt * v3;
    system->positions.noalias() = x4;
    potHydro->update();
    Eigen::VectorXd v4 = v1 + dt * a3;
    Eigen::VectorXd a4;
    acceleration_update();

    /* Output calculated values */
    system->positions.noalias()     = x1 + (dt / 6.0) * (v1 + 2.0 * v2 + 2.0 * v3 + v4);
    system->velocities.noalias()    = v1 + (dt / 6.0) * (a1 + 2.0 * a2 + 2.0 * a3 + a4);
    system->accelerations.noalias() = a4;
}

void
rungeKutta4::acceleration_update()
{
    // TODO
}