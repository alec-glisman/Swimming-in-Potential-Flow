//
// Created by Alec Glisman on 07/31/21
//

#include <rungeKutta4.hpp>

rungeKutta4::rungeKutta4(std::shared_ptr<systemData>             sys,
                         std::shared_ptr<potentialHydrodynamics> hydro)
{
    // save classes
    m_system   = sys;
    m_potHydro = hydro;

    // initialize logger
    m_logFile   = m_system->outputDir() + "/logs/" + m_logName + "-log.txt";
    auto logger = spdlog::basic_logger_mt(m_logName, m_logFile);
    spdlog::get(m_logName)->info("Initializing Runge-Kutta 4th order integrator");

    // Create integrator
    spdlog::get(m_logName)->info("Creating integrator");

    // initialize member variables
    m_dt = m_system->dt() * m_system->tau();
    spdlog::get(m_logName)->info("dt (dimensional): {0}", m_dt);
    m_c1_2_dt = m_c1_2 * m_dt;
    spdlog::get(m_logName)->info("1/2 * dt (dimensional): {0}", m_c1_2_dt);
    m_c1_6_dt = m_c1_6 * m_dt;
    spdlog::get(m_logName)->info("1/6 * dt (dimensional): {0}", m_c1_6_dt);
}

rungeKutta4::~rungeKutta4()
{
    spdlog::get(m_logName)->info("Destructing Runge-Kutta 4th order");
}

void
rungeKutta4::integrate()
{
    /* Step 1: k1 = f(t, h) */
    Eigen::VectorXd x1 = m_system->positions;
    Eigen::VectorXd v1 = m_system->velocities;
    Eigen::VectorXd a1 = m_system->accelerations;

    /* Step 2: k2 = f(t + h/2, y + h*k1/2) */
    Eigen::VectorXd x2 = v1;
    x2 *= m_c1_2_dt;
    x2.noalias() += x1;
    m_system->positions.noalias() = x2;
    m_potHydro->update();

    Eigen::VectorXd v2 = a1;
    v2 *= m_c1_2_dt;
    v2.noalias() += v1;

    Eigen::VectorXd a2 = Eigen::VectorXd::Zero(m_system->numDim() * m_system->numParticles());
    accelerationUpdate(a2);

    /* Step 3: k3 = f(t + h/2, y + h/2*k2) */
    Eigen::VectorXd x3 = v2;
    x3 *= m_c1_2_dt;
    x3.noalias() += x1;
    m_system->positions.noalias() = x3;
    m_potHydro->update();

    Eigen::VectorXd v3 = a2;
    v3 *= m_c1_2_dt;
    v3.noalias() += v1;

    Eigen::VectorXd a3 = Eigen::VectorXd::Zero(m_system->numDim() * m_system->numParticles());
    accelerationUpdate(a3);

    /* Step 4: k4 = f(t + h, y + h*k3) */
    Eigen::VectorXd x4 = v3;
    x4 *= m_c1_2_dt;
    x4.noalias() += x1;
    m_system->positions.noalias() = x4;
    m_potHydro->update();

    Eigen::VectorXd v4 = a3;
    v4 *= m_c1_2_dt;
    v4.noalias() += v1;

    Eigen::VectorXd a4 = Eigen::VectorXd::Zero(m_system->numDim() * m_system->numParticles());
    accelerationUpdate(a4);

    /* Output calculated values */
    m_system->positions.noalias() = v1;
    m_system->positions.noalias() += m_c2 * v2;
    m_system->positions.noalias() += m_c2 * v3;
    m_system->positions.noalias() += v4;
    m_system->positions *= m_c1_6_dt;
    m_system->positions.noalias() += x1;

    m_system->velocities.noalias() = a1;
    m_system->velocities.noalias() += m_c2 * a2;
    m_system->velocities.noalias() += m_c2 * a3;
    m_system->velocities.noalias() += a4;
    m_system->velocities *= m_c1_6_dt;
    m_system->velocities.noalias() += v1;

    m_system->accelerations.noalias() = a4;
}

void
rungeKutta4::accelerationUpdate(Eigen::VectorXd& acc)
{
    /*! Solve linear equation of the form: Ax = b, where x is the unknown acceleration vector
     * A is the total mass matrix, and b are the known parts of the forces */
    Eigen::VectorXd f = Eigen::VectorXd::Zero(m_system->numDim() * m_system->numParticles());

    if (m_system->particleDensity() >= 0) // hydrodynamic force
    {
        f.noalias() += m_potHydro->fHydroNoInertia();
    }

    acc.noalias() = m_potHydro->mTotal().llt().solve(f);
}

const Eigen::VectorXd&
rungeKutta4::articulationVel()
{
    double dimensional_time = m_system->t() * m_system->tau();
    m_velArtic              = Eigen::VectorXd::Zero(3 * m_system->numParticles());

    double part0_x_vel = -bondU0 * cos(bondOmega * dimensional_time);
    double part2_x_vel = bondU0 * cos(bondOmega * dimensional_time + bondPhaseAngle);

    m_velArtic(0)     = part0_x_vel;
    m_velArtic(3 * 2) = part2_x_vel;

    return m_velArtic;
}

const Eigen::VectorXd&
rungeKutta4::articulationAcc()
{
    double dimensional_time = m_system->t() * m_system->tau();
    m_accArtic              = Eigen::VectorXd::Zero(3 * m_system->numParticles());

    double part0_x_acc = bondU0 * bondOmega * sin(bondOmega * dimensional_time);
    double part2_x_acc = -bondU0 * bondOmega * sin(bondOmega * dimensional_time + bondPhaseAngle);

    m_accArtic(0)     = part0_x_acc;
    m_accArtic(3 * 2) = part2_x_acc;

    return m_accArtic;
}

const Eigen::Vector3d&
rungeKutta4::rLoc()
{
    m_RLoc = m_system->positions.segment<3>(3 * 1);
    return m_RLoc;
}