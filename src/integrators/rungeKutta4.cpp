//
// Created by Alec Glisman on 07/31/21
//

#include <rungeKutta4.hpp>

rungeKutta4::rungeKutta4(systemData& sys)
{
    // save classes
    system = &sys;

    // Initialize logger
    m_logFile   = system->outputDir() + "/logs/" + m_logName + "-log.txt";
    auto logger = spdlog::basic_logger_mt(m_logName, m_logFile);
    spdlog::get(m_logName)->info("Initializing Runge-Kutta 4th order integrator");

    // Create integrator
    spdlog::get(m_logName)->info("Creating integrator");
}

rungeKutta4::rungeKutta4(systemData& sys, potentialHydrodynamics& hydro) : rungeKutta4(sys)
{
    // save classes
    potHydro = &hydro;
}

rungeKutta4::~rungeKutta4()
{
    spdlog::get(m_logName)->info("Destructing Runge-Kutta 4th order");
}

void
rungeKutta4::integrate()
{
    /* Time unit conversion */
    double t  = system->timestep() * system->Tau();
    double dt = system->dt() * system->Tau();

    /* Step 1: k1 = f(t, h) */
    Eigen::VectorXd x1 = *system->positions();
    Eigen::VectorXd v1 = *system->velocities();
    Eigen::VectorXd a1 = *system->accelerations();

    /* Step 2: k2 = f(t + h/2, y + h*k1/2) */
    Eigen::VectorXd x2       = x1 + (0.5 * dt) * v1;
    sys->positions.noalias() = x2;
    updateHydrodynamics();
    Eigen::VectorXd v2 = v1 + (0.5 * dt) * a1;
    Eigen::VectorXd a2;
    accVec(x2, v2, a2, t + (0.5 * dt));

    /* Step 3: k3 = f(t + h/2, y + h/2*k2) */
    Eigen::VectorXd x3       = x1 + (0.5 * dt) * v2;
    sys->positions.noalias() = x3;
    updateHydrodynamics();
    Eigen::VectorXd v3 = v1 + (0.5 * dt) * a2;
    Eigen::VectorXd a3;
    accVec(x3, v3, a3, tND + (0.5 * dt));

    /* Step 4: k4 = f(t + h, y + h*k3) */
    Eigen::VectorXd x4       = x1 + (dt)*v3;
    sys->positions.noalias() = x4;
    updateHydrodynamics();
    Eigen::VectorXd v4 = v1 + (dt)*a3;
    Eigen::VectorXd a4;
    accVec(x4, v4, a4, tND + (dt));

    /* Output calculated values */
    sys->positions.noalias()     = x1 + (dt / 6.0) * (v1 + 2.0 * v2 + 2.0 * v3 + v4);
    sys->velocities.noalias()    = v1 + (dt / 6.0) * (a1 + 2.0 * a2 + 2.0 * a3 + a4);
    sys->accelerations.noalias() = a4;
}

void
rungeKutta4::acceleration_update()
{
}