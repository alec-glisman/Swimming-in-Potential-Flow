//
// Created by Alec Glisman on 07/31/21
//

#include <rungeKutta4.hpp>

rungeKutta4::rungeKutta4(systemData& sys, potentialHydrodynamics& hydro)
{
    // save classes
    m_system   = &sys;
    m_potHydro = &hydro;

    // Initialize logger
    m_logFile   = m_system->outputDir() + "/logs/" + m_logName + "-log.txt";
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
    double t  = m_system->timestep() * m_system->tau();
    double dt = m_system->dt() * m_system->tau();

    /* Step 1: k1 = f(t, h) */
    Eigen::VectorXd x1 = m_system->positions;
    Eigen::VectorXd v1 = m_system->velocities;
    Eigen::VectorXd a1 = m_system->accelerations;

    /* Step 2: k2 = f(t + h/2, y + h*k1/2) */
    Eigen::VectorXd x2          = x1 + (0.5 * dt) * v1;
    m_system->positions.noalias() = x2;
    m_potHydro->update();
    Eigen::VectorXd v2 = v1 + (0.5 * dt) * a1;
    Eigen::VectorXd a2 = Eigen::VectorXd::Zero(m_system->numDim() * m_system->numParticles());
    acceleration_update(a2);

    /* Step 3: k3 = f(t + h/2, y + h/2*k2) */
    Eigen::VectorXd x3          = x1 + (0.5 * dt) * v2;
    m_system->positions.noalias() = x3;
    m_potHydro->update();
    Eigen::VectorXd v3 = v1 + (0.5 * dt) * a2;
    Eigen::VectorXd a3 = Eigen::VectorXd::Zero(m_system->numDim() * m_system->numParticles());
    acceleration_update(a3);

    /* Step 4: k4 = f(t + h, y + h*k3) */
    Eigen::VectorXd x4          = x1 + dt * v3;
    m_system->positions.noalias() = x4;
    m_potHydro->update();
    Eigen::VectorXd v4 = v1 + dt * a3;
    Eigen::VectorXd a4 = Eigen::VectorXd::Zero(m_system->numDim() * m_system->numParticles());
    acceleration_update(a4);

    /* Output calculated values */
    m_system->positions.noalias()     = x1 + (dt / 6.0) * (v1 + 2.0 * v2 + 2.0 * v3 + v4);
    m_system->velocities.noalias()    = v1 + (dt / 6.0) * (a1 + 2.0 * a2 + 2.0 * a3 + a4);
    m_system->accelerations.noalias() = a4;
}

void
rungeKutta4::acceleration_update(Eigen::VectorXd& acc)
{
    // TODO
}