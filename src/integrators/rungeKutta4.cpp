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

    spdlog::get(m_logName)->info("Constructor complete");
    spdlog::get(m_logName)->flush();
}

rungeKutta4::~rungeKutta4()
{
    spdlog::get(m_logName)->info("Destructing Runge-Kutta 4th order");
    spdlog::get(m_logName)->flush();
}

void
rungeKutta4::integrate()
{
    integrateSecondOrder(); // Udwadia-Kalaba method only gives acceleration components
}

void
rungeKutta4::integrateSecondOrder()
{
    /* ANCHOR: Solve system of form: y'(t) = f( y(t),  t )
     * @REFERENCE:
     * https://www.physicsforums.com/threads/using-runge-kutta-method-for-position-calc.553663/post-3634957
     */

    /* Step 1: k1 = f( y(t_0),  t_0 ),
     * initial conditions at current step */
    const Eigen::VectorXd v1 = m_system->velocitiesParticles();
    const Eigen::VectorXd x1 = m_system->positionsParticles();
    Eigen::VectorXd       a1 = Eigen::VectorXd::Zero(3 * m_system->numParticles());
    accelerationUpdate(a1, m_system->t());

    /* Step 2: k2 = f( y(t_0) + k1 * dt/2,  t_0 + dt/2 )
     * time rate-of-change k1 evaluated halfway through time step (midpoint) */
    Eigen::VectorXd v2 = v1;
    v2.noalias() += m_c1_2_dt * a1;
    Eigen::VectorXd x2 = x1;
    x2.noalias() += m_c1_2_dt * v2;

    m_system->setVelocitiesParticles(v2);
    m_system->setPositionsParticles(x2);

    Eigen::VectorXd a2 = Eigen::VectorXd::Zero(3 * m_system->numParticles());
    accelerationUpdate(a2, m_system->t() + 0.50 * m_system->dt());

    /* Step 3: k3 = f( y(t_0) + k2 * dt/2,  t_0 + dt/2 )
     * time rate-of-change k2 evaluated halfway through time step (midpoint) */
    Eigen::VectorXd v3 = v1;
    v3.noalias() += m_c1_2_dt * a2;
    Eigen::VectorXd x3 = x1;
    x3.noalias() += m_c1_2_dt * v3;

    m_system->setVelocitiesParticles(v3);
    m_system->setPositionsParticles(x3);

    Eigen::VectorXd a3 = Eigen::VectorXd::Zero(3 * m_system->numParticles());
    accelerationUpdate(a3, m_system->t() + 0.50 * m_system->dt());

    /* Step 4: k4 = f( y(t_0) + k3 * dt,  t_0 + dt )
     * time rate-of-change k3 evaluated at end of step (endpoint) */
    Eigen::VectorXd v4 = v1;
    v4.noalias() += m_dt * a3;
    Eigen::VectorXd x4 = x1;
    x4.noalias() += m_dt * v4;

    m_system->setVelocitiesParticles(v4);
    m_system->setPositionsParticles(x4);

    Eigen::VectorXd a4 = Eigen::VectorXd::Zero(3 * m_system->numParticles());
    accelerationUpdate(a4, m_system->t() + m_system->dt());

    /* ANCHOR: Calculate kinematics at end of time step */
    Eigen::VectorXd v_out = a1;
    v_out.noalias() += 2.0 * a2;
    v_out.noalias() += 2.0 * a3;
    v_out.noalias() += a4;
    v_out *= m_c1_6_dt;
    v_out.noalias() += v1;

    Eigen::VectorXd x_out = v1;
    x_out.noalias() += 2.0 * v2;
    x_out.noalias() += 2.0 * v3;
    x_out.noalias() += v4;
    x_out *= m_c1_6_dt;
    x_out.noalias() += x1;

    m_system->setVelocitiesParticles(v_out);
    m_system->setPositionsParticles(x_out);

    Eigen::VectorXd a_out = Eigen::VectorXd::Zero(3 * m_system->numParticles());
    accelerationUpdate(a_out, m_system->t() + m_system->dt());
    m_system->setAccelerationsParticles(a_out);
}

void
rungeKutta4::integrateFirstOrder()
{
    const double    time{m_system->t() * m_system->tau()};
    Eigen::VectorXd unused = Eigen::VectorXd::Zero(3 * m_system->numParticles());

    /* Step 1: k1 = f( y(t_0),  t_0 ),
     * initial conditions at current step */
    const Eigen::VectorXd x1 = m_system->positionsParticles();
    accelerationUpdate(unused, time);
    const Eigen::VectorXd v1 = m_system->velocitiesParticles();

    /* Step 2: k2 = f( y(t_0) + k1 * dt/2,  t_0 + dt/2 )
     * time rate-of-change k1 evaluated halfway through time step (midpoint) */
    Eigen::VectorXd x2 = x1;
    x2.noalias() += m_c1_2_dt * v1;

    m_system->setPositionsParticles(x2);

    accelerationUpdate(unused, time);
    const Eigen::VectorXd v2 = m_system->velocitiesParticles();

    /* Step 3: k3 = f( y(t_0) + k2 * dt/2,  t_0 + dt/2 )
     * time rate-of-change k2 evaluated halfway through time step (midpoint) */
    Eigen::VectorXd x3 = x1;
    x3.noalias() += m_c1_2_dt * v2;

    m_system->setPositionsParticles(x3);

    accelerationUpdate(unused, time);
    const Eigen::VectorXd v3 = m_system->velocitiesParticles();

    /* Step 4: k4 = f( y(t_0) + k3 * dt,  t_0 + dt )
     * time rate-of-change k3 evaluated at end of step (endpoint) */
    Eigen::VectorXd x4 = x1;
    x4.noalias() += m_dt * v3;

    m_system->setPositionsParticles(x4);

    accelerationUpdate(unused, time);
    const Eigen::VectorXd v4 = m_system->velocitiesParticles();

    /* ANCHOR: Calculate kinematics at end of time step */
    Eigen::VectorXd x_out = v1;
    x_out.noalias() += 2.0 * v2;
    x_out.noalias() += 2.0 * v3;
    x_out.noalias() += v4;
    x_out *= m_c1_6_dt;
    x_out.noalias() += x1;

    m_system->setPositionsParticles(x_out);

    accelerationUpdate(unused, time);
}

void
rungeKutta4::accelerationUpdate(Eigen::VectorXd& acc, double time)
{
    m_potHydro->update();              // update hydrodynamic tensors and forces
    m_system->update(time); // update Udwadia linear system

    udwadiaKalaba(acc); // solve for acceleration components
}

void
rungeKutta4::udwadiaKalaba(Eigen::VectorXd& acc)
{
    /* NOTE: Following the formalism developed in Udwadia & Kalaba (1992) Proc. R. Soc. Lond. A
     * Solve system of the form M_total * acc = Q + Q_con
     * Q is the forces present in unconstrained system
     * Q_con is the generalized constraint forces */

    // calculate Q
    Eigen::VectorXd Q = Eigen::VectorXd::Zero(3 * m_system->numParticles());
    if (m_system->fluidDensity() > 0) // hydrodynamic force
    {
        Q.noalias() += m_potHydro->fHydroNoInertia();
    }

    /* calculate Q_con = K (b - A * M_total^{-1} * Q)
     * Linear constraint system: A * acc = b
     * Linear proportionality: K = M_total^{1/2} * (A * M_total^{-1/2})^{+};
     * + is Moore-Penrose inverse */

    // calculate M^{1/2} & M^{-1/2}
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(m_potHydro->mTotal());
    if (eigensolver.info() != Eigen::Success)
    {
        spdlog::get(m_logName)->error("Computing eigendecomposition of M_total failed at t={0}",
                                      m_system->t());
        throw std::runtime_error("Computing eigendecomposition of M_total failed");
    }
    const Eigen::MatrixXd M_total_halfPower         = eigensolver.operatorSqrt();
    const Eigen::MatrixXd M_total_negativeHalfPower = eigensolver.operatorInverseSqrt();

    // calculate K
    const Eigen::MatrixXd AM_nHalf = m_system->udwadiaA() * M_total_negativeHalfPower;
    const Eigen::MatrixXd AM_nHalf_pInv =
        AM_nHalf.completeOrthogonalDecomposition().pseudoInverse();
    const Eigen::MatrixXd K = M_total_halfPower * AM_nHalf_pInv;

    // calculate Q_con
    const Eigen::MatrixXd M_total_inv   = m_potHydro->mTotal().inverse();
    const Eigen::MatrixXd M_total_invQ  = M_total_inv * Q;
    const Eigen::VectorXd AM_total_invQ = m_system->udwadiaA() * M_total_invQ;
    Eigen::VectorXd       b_tilde       = m_system->udwadiaB();
    b_tilde.noalias() -= AM_total_invQ;
    const Eigen::VectorXd Q_con = K * b_tilde;

    // calculate accelerations
    Eigen::VectorXd Q_total = Q;
    Q_total.noalias() += Q_con;
    acc.noalias() = m_potHydro->mTotal().llt().solve(Q_total);
}
