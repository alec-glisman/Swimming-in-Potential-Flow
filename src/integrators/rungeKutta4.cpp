//
// Created by Alec Glisman on 07/31/21
//

#include <rungeKutta4.hpp>

rungeKutta4::rungeKutta4(std::shared_ptr<systemData> sys, std::shared_ptr<potentialHydrodynamics> hydro)
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

    m_3M = 3 * m_system->numBodies();

    spdlog::get(m_logName)->info("Constructor complete");
    spdlog::get(m_logName)->flush();
}

rungeKutta4::~rungeKutta4()
{
    spdlog::get(m_logName)->info("Destructing Runge-Kutta 4th order");
    spdlog::get(m_logName)->flush();
}

void
rungeKutta4::integrate(Eigen::ThreadPoolDevice& device)
{
    integrateSecondOrder(device); // Udwadia-Kalaba method only gives acceleration components
}

void
rungeKutta4::integrateSecondOrder(Eigen::ThreadPoolDevice& device)
{
    /* Step 1: k1 = f( y(t_0),  t_0 ),
     * initial conditions at current step */
    const double          t1{m_system->t()};
    const Eigen::VectorXd v1 = m_system->velocitiesParticles();
    const Eigen::VectorXd x1 = m_system->positionsParticles();

    Eigen::VectorXd a1 = Eigen::VectorXd::Zero(m_3M);
    accelerationUpdate(t1, a1, device);

    /* Step 2: k2 = f( y(t_0) + k1 * dt/2,  t_0 + dt/2 )
     * time rate-of-change k1 evaluated halfway through time step (midpoint) */
    const double    t2{t1 + 0.50 * m_system->dt()};
    Eigen::VectorXd v2 = v1;
    v2.noalias() += m_c1_2_dt * a1;
    Eigen::VectorXd x2 = x1;
    x2.noalias() += m_c1_2_dt * v2;

    Eigen::VectorXd a2 = Eigen::VectorXd::Zero(m_3M);
    accelerationUpdate(t2, a2, device);

    /* Step 3: k3 = f( y(t_0) + k2 * dt/2,  t_0 + dt/2 )
     * time rate-of-change k2 evaluated halfway through time step (midpoint) */
    const double    t3{t2};
    Eigen::VectorXd v3 = v1;
    v3.noalias() += m_c1_2_dt * a2;
    Eigen::VectorXd x3 = x1;
    x3.noalias() += m_c1_2_dt * v3;

    Eigen::VectorXd a3 = Eigen::VectorXd::Zero(m_3M);
    accelerationUpdate(t3, a3, device);

    /* Step 4: k4 = f( y(t_0) + k3 * dt,  t_0 + dt )
     * time rate-of-change k3 evaluated at end of step (endpoint) */
    const double    t4{t1 + m_system->dt()};
    Eigen::VectorXd v4 = v1;
    v4.noalias() += m_dt * a3;
    Eigen::VectorXd x4 = x1;
    x4.noalias() += m_dt * v4;

    Eigen::VectorXd a4 = Eigen::VectorXd::Zero(m_3M);
    accelerationUpdate(t4, a4, device);

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

    Eigen::VectorXd a_out = Eigen::VectorXd::Zero(m_3M);
    accelerationUpdate(t4, a_out, device);

    // reset system time to t1 of this integration step as `engine` class manages updating system time at end of each
    // step
    m_system->setT(t1);
}

void
rungeKutta4::accelerationUpdate(const double t, Eigen::VectorXd& acc, Eigen::ThreadPoolDevice& device)
{
    // NOTE: Order of all 4 function calls must remain the same
    m_system->setT(t);
    m_system->update(device);
    m_potHydro->update(device);

    udwadiaKalaba(acc);
}

void
rungeKutta4::udwadiaKalaba(Eigen::VectorXd& acc)
{
    /* NOTE: Following the formalism developed in Udwadia & Kalaba (1992) Proc. R. Soc. Lond. A
     * Solve system of the form M_total * acc = Q + Q_con
     * Q is the forces present in unconstrained system
     * Q_con is the generalized constraint forces */

    // calculate Q
    Eigen::VectorXd Q = Eigen::VectorXd::Zero(7 * m_system->numBodies());

    // hydrodynamic force
    if (m_system->fluidDensity() > 0)
    {
        Q.noalias() += m_potHydro->fHydroNoInertia();
    }

    /* calculate Q_con = K (b - A * M_total^{-1} * Q)
     * Linear constraint system: A * acc = b
     * Linear proportionality: K = M_total^{1/2} * (A * M_total^{-1/2})^{+};
     * + is Moore-Penrose inverse */

    // calculate M^{1/2} & M^{-1/2}
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(m_potHydro->mTilde()); // = M
    if (eigensolver.info() != Eigen::Success)
    {
        spdlog::get(m_logName)->error("Computing eigendecomposition of effective total mass matrix failed at t={0}",
                                      m_system->t());
        throw std::runtime_error("Computing eigendecomposition of effective total mass matrix failed");
    }
    const Eigen::MatrixXd M_total_halfPower         = eigensolver.operatorSqrt();
    const Eigen::MatrixXd M_total_negativeHalfPower = eigensolver.operatorInverseSqrt();

    // calculate K
    const Eigen::MatrixXd AM_nHalf      = m_system->udwadiaA() * M_total_negativeHalfPower;
    const Eigen::MatrixXd AM_nHalf_pInv = AM_nHalf.completeOrthogonalDecomposition().pseudoInverse();
    const Eigen::MatrixXd K             = M_total_halfPower * AM_nHalf_pInv;

    // calculate Q_con
    const Eigen::MatrixXd M_total_inv   = m_potHydro->mTilde().inverse();
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
