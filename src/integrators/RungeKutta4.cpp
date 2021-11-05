//
// Created by Alec Glisman on 07/31/21
//

#include <RungeKutta4.hpp>

RungeKutta4::RungeKutta4(std::shared_ptr<SystemData> sys, std::shared_ptr<PotentialHydrodynamics> hydro)
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

    m_7M = 7 * m_system->numBodies();

    spdlog::get(m_logName)->info("Constructor complete");
    spdlog::get(m_logName)->flush();
}

RungeKutta4::~RungeKutta4()
{
    spdlog::get(m_logName)->info("Destructing Runge-Kutta 4th order");
    spdlog::get(m_logName)->flush();
    spdlog::drop(m_logName);
}

void
RungeKutta4::integrate(Eigen::ThreadPoolDevice& device)
{
    integrateSecondOrder(device); // Udwadia-Kalaba method only gives acceleration components
}

void
RungeKutta4::integrateSecondOrder(Eigen::ThreadPoolDevice& device)
{
    /* Step 1: k1 = f( y(t_0),  t_0 ),
     * initial conditions at current step */
    const double    t1{m_system->t()};
    Eigen::VectorXd v1 = m_system->velocitiesBodies();
    Eigen::VectorXd x1 = m_system->positionsBodies();

    Eigen::VectorXd a1 = Eigen::VectorXd::Zero(m_7M);
    accelerationUpdate(t1, x1, v1, a1, device);

    /* Step 2: k2 = f( y(t_0) + k1 * dt/2,  t_0 + dt/2 )
     * time rate-of-change k1 evaluated halfway through time step (midpoint) */
    const double    t2{t1 + 0.50 * m_system->dt()};
    Eigen::VectorXd v2 = v1;
    v2.noalias() += m_c1_2_dt * a1;
    Eigen::VectorXd x2 = x1;
    x2.noalias() += m_c1_2_dt * v2;

    Eigen::VectorXd a2 = Eigen::VectorXd::Zero(m_7M);
    accelerationUpdate(t2, x2, v2, a2, device);

    /* Step 3: k3 = f( y(t_0) + k2 * dt/2,  t_0 + dt/2 )
     * time rate-of-change k2 evaluated halfway through time step (midpoint) */
    const double    t3{t2};
    Eigen::VectorXd v3 = v1;
    v3.noalias() += m_c1_2_dt * a2;
    Eigen::VectorXd x3 = x1;
    x3.noalias() += m_c1_2_dt * v3;

    Eigen::VectorXd a3 = Eigen::VectorXd::Zero(m_7M);
    accelerationUpdate(t3, x3, v3, a3, device);

    /* Step 4: k4 = f( y(t_0) + k3 * dt,  t_0 + dt )
     * time rate-of-change k3 evaluated at end of step (endpoint) */
    const double    t4{t1 + m_system->dt()};
    Eigen::VectorXd v4 = v1;
    v4.noalias() += m_dt * a3;
    Eigen::VectorXd x4 = x1;
    x4.noalias() += m_dt * v4;

    Eigen::VectorXd a4 = Eigen::VectorXd::Zero(m_7M);
    accelerationUpdate(t4, x4, v4, a4, device);

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

    Eigen::VectorXd a_out = Eigen::VectorXd::Zero(m_7M);
    accelerationUpdate(t4, x_out, v_out, a_out, device);

    // reset system time to t1 as `Engine` class manages updating system time at end of each step
    m_system->setT(t1);
}

void
RungeKutta4::accelerationUpdate(const double t, Eigen::VectorXd& pos, Eigen::VectorXd& vel, Eigen::VectorXd& acc,
                                Eigen::ThreadPoolDevice& device)
{
    // NOTE: Order of function calls must remain the same
    m_system->setT(t);

    // TODO: Related to Issue #3
    if (m_system->imageSystem())
    {
        imageBodyPosVel(pos, vel);
    }
    m_system->setPositionsBodies(pos);
    m_system->setVelocitiesBodies(vel);

    m_system->update(device);
    m_potHydro->update(device);

    // TODO: Related to Issue #3
    if (m_system->imageSystem())
    {
        // only update Udwadia system for "real" system not image dipoles
        const int num_img_bodies{m_system->numBodies() / 2};
        const int num_img_bodies_7{7 * num_img_bodies};

        Eigen::VectorXd acc_real_body = acc.segment(0, num_img_bodies_7);

        udwadiaKalaba(acc_real_body);

        // update acceleration components using constraints
        acc.segment(0, num_img_bodies_7).noalias() = acc_real_body;
        imageBodyAcc(acc);
    }
    else
    {
        // update Udwadia system for all bodies
        udwadiaKalaba(acc);
    }

    m_system->setAccelerationsBodies(acc);
}

void
RungeKutta4::imageBodyPosVel(Eigen::VectorXd& pos, Eigen::VectorXd& vel)
{
    const int num_img_bodies = m_system->numBodies() / 2;

    for (int img_body_id = num_img_bodies; img_body_id < m_system->numBodies(); img_body_id++)
    {
        // indices
        const int body_id_7{7 * (img_body_id - num_img_bodies)};
        const int img_body_id_7{7 * img_body_id};

        /* ANCHOR: Position image system: mirror image about xy-plane (transform z --> -z) */
        pos.segment<7>(img_body_id_7) = pos.segment<7>(body_id_7);

        // (linear components) flip z component, leave x-y unchanged
        pos.segment<1>(img_body_id_7 + 2) *= -1;

        // (quaternion components) flip x-y components of quaternion vector part (see note in header file: look for "C1
        // and C2 are mirrored along the Z-axis")
        pos.segment<2>(img_body_id_7 + 4) *= -1;

        /* ANCHOR: Velocity image system: invert x-y components to have image dipoles  */
        vel.segment<7>(img_body_id_7) = vel.segment<7>(body_id_7);

        // (linear components) flip x-y components, leave z unchanged
        vel.segment<2>(img_body_id_7) *= -1;

        // (quaternion components) same transform as positional transform
        vel.segment<2>(img_body_id_7 + 4) *= -1;
    }
}

void
RungeKutta4::imageBodyAcc(Eigen::VectorXd& acc)
{
    const int num_img_bodies = m_system->numBodies() / 2;

    for (int img_body_id = num_img_bodies; img_body_id < m_system->numBodies(); img_body_id++)
    {
        // indices
        const int body_id_7{7 * (img_body_id - num_img_bodies)};
        const int img_body_id_7{7 * img_body_id};

        /* ANCHOR: Acceleration image system: invert x-y components to have image dipoles  */
        acc.segment<7>(img_body_id_7) = acc.segment<7>(body_id_7);

        // (linear components) flip x-y components, leave z unchanged
        acc.segment<2>(img_body_id_7) *= -1;

        // (quaternion components) same transform as positional transform
        acc.segment<2>(img_body_id_7 + 4) *= -1;
    }
}

void
RungeKutta4::udwadiaKalaba(Eigen::VectorXd& acc)
{
    /* NOTE: Following the formalism developed in Udwadia & Kalaba (1992) Proc. R. Soc. Lond. A
     * Solve system of the form M_eff * acc = Q + Q_con
     * Q is the forces present in unconstrained system
     * Q_con is the generalized constraint forces */

    int body_dof{m_system->numBodies()}; // = m

    if (m_system->imageSystem())
    {
        body_dof /= 2;
    }

    const int body_dof_7{body_dof * 7};

    // calculate Q
    Eigen::VectorXd Q = Eigen::VectorXd::Zero(body_dof_7); // (7m, 1)

    // hydrodynamic force
    if (m_system->fluidDensity() > 0)
    {
        Q.noalias() += m_potHydro->fHydroNoInertia().segment(0, body_dof_7);
    }

    /* calculate Q_con = K (b - A * M_eff^{-1} * Q)
     * Linear constraint system: A * acc = b
     * Linear proportionality: K = M_eff^{1/2} * (A * M_eff^{-1/2})^{+};
     * + is Moore-Penrose inverse */

    const Eigen::MatrixXd M_eff = m_potHydro->mTotalBodyCoords().block(0, 0, body_dof_7, body_dof_7); // (7m, 7m)

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(M_eff);

    // compute eigen-decomposition of M_eff
    if (eigensolver.info() != Eigen::Success)
    {
        spdlog::get(m_logName)->error("Computing eigendecomposition of effective total mass matrix failed at t={0}",
                                      m_system->t());
        throw std::runtime_error("Computing eigendecomposition of effective total mass matrix failed");
    }

    // verify M is positive semi-definite
#if !defined(NDEBUG)
    if (!M_eff.isApprox(M_eff.transpose()) || (eigensolver.eigenvalues().array() < 0.0).any())
    {
        throw std::runtime_error("M_eff is non semi-positive definite matrix!");
    }
#endif

    // calculate M^{1/2} & M^{-1/2}
    const Eigen::MatrixXd M_eff_halfPower         = eigensolver.operatorSqrt();        // (7m, 7m)
    const Eigen::MatrixXd M_eff_negativeHalfPower = eigensolver.operatorInverseSqrt(); // (7m, 7m)

    // calculate K
    const Eigen::MatrixXd AM_nHalf      = m_system->udwadiaA() * M_eff_negativeHalfPower; // (c, 7m) = (c, 7m) (7m, 7m)
    const Eigen::MatrixXd AM_nHalf_pInv = AM_nHalf.completeOrthogonalDecomposition().pseudoInverse(); // (7m, c)
    const Eigen::MatrixXd K             = M_eff_halfPower * AM_nHalf_pInv;                            // (7m, c)

    // calculate Q_con
    const Eigen::MatrixXd M_eff_inv   = M_eff.inverse();                   // (7m, 7m)
    const Eigen::MatrixXd M_eff_invQ  = M_eff_inv * Q;                     // (7m, 1)
    const Eigen::VectorXd AM_eff_invQ = m_system->udwadiaA() * M_eff_invQ; // (c, 1) = (c, 7m) (7m, 1)
    Eigen::VectorXd       b_tilde     = m_system->udwadiaB();              // (c, 1)
    b_tilde.noalias() -= AM_eff_invQ;
    const Eigen::VectorXd Q_con = K * b_tilde; // (7m, 1) = (7m, c) (c, 1)

    // calculate accelerations
    Eigen::VectorXd Q_total = Q; // (7m, 1)
    Q_total.noalias() += Q_con;
    acc.noalias() = M_eff_inv * Q_total; // (7m, 1)

    // debugging print statements
    Eigen::IOFormat CleanFmt(4, 0, ", ", "\n", "[", "]");
    std::cout << "M_eff:\n" << M_eff.format(CleanFmt) << "\n\n" << std::endl;
    std::cout << "M_eff_inv:\n" << M_eff_inv.format(CleanFmt) << "\n\n" << std::endl;
    std::cout << "K:\n" << K.format(CleanFmt) << "\n\n" << std::endl;
    std::cout << "udwadiaA:\n" << m_system->udwadiaA().format(CleanFmt) << "\n\n" << std::endl;
    std::cout << "udwadiaB:\n" << m_system->udwadiaB().format(CleanFmt) << "\n\n" << std::endl;
    std::cout << "acc:\n" << acc.format(CleanFmt) << "\n\n" << std::endl;
}
