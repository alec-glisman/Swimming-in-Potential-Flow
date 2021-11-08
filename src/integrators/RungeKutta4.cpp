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

    // create temporary Eigen device to use for tensor computations in constructor
    spdlog::get(m_logName)->critical("Setting up temporary single-thread eigen device");
    Eigen::ThreadPool       thread_pool = Eigen::ThreadPool(1);
    Eigen::ThreadPoolDevice single_core_device(&thread_pool, 1);

    // Create integrator
    spdlog::get(m_logName)->info("Creating integrator");

    // initialize member variables
    spdlog::get(m_logName)->info("Setting attributes");

    m_dt = m_system->dt() * m_system->tau();
    spdlog::get(m_logName)->info("dt (dimensional): {0}", m_dt);
    m_c1_2_dt = m_c1_2 * m_dt;
    spdlog::get(m_logName)->info("1/2 * dt (dimensional): {0}", m_c1_2_dt);
    m_c1_6_dt = m_c1_6 * m_dt;
    spdlog::get(m_logName)->info("1/6 * dt (dimensional): {0}", m_c1_6_dt);

    m_7M = 7 * m_system->numBodies();
    spdlog::get(m_logName)->info("7M: {0}", m_7M);

    // degrees of freedom: divide by 2 if using image system
    m_body_dof = m_system->numBodies(); // = m
    if (m_system->imageSystem())
    {
        m_body_dof /= 2;
    }
    m_body_dof_7 = 7 * m_body_dof;

    spdlog::get(m_logName)->info("body_dof: {0}", m_body_dof);
    spdlog::get(m_logName)->info("body_dof_7: {0}", m_body_dof_7);

    // for initial conditions of all locater points, use the PLFT-free algorithm
    spdlog::get(m_logName)->critical("Setting initial conditions using PLFT-free algorithm.");
    m_system->update(single_core_device);   // update system kinematics and rbm tensors
    m_potHydro->update(single_core_device); // update hydrodynamic tensors
    momForceFree();                         // set initial body kinematics

    if (m_system->imageSystem())
    {
        spdlog::get(m_logName)->info("Setting image position, velocity, and acceleration, using real components.");

        Eigen::VectorXd body_pos = m_system->positionsBodies();
        Eigen::VectorXd body_vel = m_system->velocitiesBodies();
        Eigen::VectorXd body_acc = m_system->accelerationsBodies();

        imageBodyPosVel(body_pos, body_vel);
        imageBodyAcc(body_acc);

        m_system->setPositionsBodies(body_pos);
        m_system->setVelocitiesBodies(body_vel);
        m_system->setAccelerationsBodies(body_acc);
    }

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

    if (m_system->imageSystem())
    {
        imageBodyPosVel(pos, vel);
    }
    m_system->setPositionsBodies(pos);
    m_system->setVelocitiesBodies(vel);

    m_system->update(device);
    m_potHydro->update(device);

    if (m_system->imageSystem())
    {
        // only update Udwadia system for "real" system not image dipoles
        const int num_real_bodies{m_system->numBodies() / 2};
        const int num_real_bodies_7{7 * num_real_bodies};

        Eigen::VectorXd acc_real_body = acc.segment(0, num_real_bodies_7);

        udwadiaKalaba(acc_real_body);

        // update acceleration components using constraints
        acc.segment(0, num_real_bodies_7).noalias() = acc_real_body;
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
    const int num_img_bodies{m_system->numBodies() / 2};

    for (int img_body_id = num_img_bodies; img_body_id < m_system->numBodies(); img_body_id++)
    {
        // `Eigen::Vector` indices
        const int real_body_id_7{7 * (img_body_id - num_img_bodies)};
        const int img_body_id_7{7 * img_body_id};

        /* ANCHOR: Position image system: mirror image about xy-plane (transform z --> -z) */
        pos.segment<7>(img_body_id_7).noalias() = pos.segment<7>(real_body_id_7);

        // (linear components) flip z component, leave x-y unchanged
        pos(img_body_id_7 + 2) *= -1;

        // (quaternion components) flip x-y components of quaternion vector part (see note in header file: look for "C1
        // and C2 are mirrored along the z-axis")
        pos.segment<2>(img_body_id_7 + 4) *= -1;

        /* ANCHOR: Velocity image system: mirror image about xy-plane (transform z --> -z)  */
        vel.segment<7>(img_body_id_7).noalias() = vel.segment<7>(real_body_id_7);

        vel(img_body_id_7 + 2) *= -1;
        vel.segment<2>(img_body_id_7 + 4) *= -1; // FIXME: Is this true for first time derivative?
    }
}

void
RungeKutta4::imageBodyAcc(Eigen::VectorXd& acc)
{
    const int num_img_bodies{m_system->numBodies() / 2};

    for (int img_body_id = num_img_bodies; img_body_id < m_system->numBodies(); img_body_id++)
    {
        // `Eigen::Vector` indices
        const int real_body_id_7{7 * (img_body_id - num_img_bodies)};
        const int img_body_id_7{7 * img_body_id};

        /* ANCHOR: Acceleration image system: mirror image about xy-plane (transform z --> -z)   */
        acc.segment<7>(img_body_id_7).noalias() = acc.segment<7>(real_body_id_7);

        acc(img_body_id_7 + 2) *= -1;
        acc.segment<2>(img_body_id_7 + 4) *= -1; // FIXME: Is this true for second time derivative?
    }
}

void
RungeKutta4::udwadiaKalaba(Eigen::VectorXd& acc)
{
    /* NOTE: Following the formalism developed in Udwadia & Kalaba (1992) Proc. R. Soc. Lond. A
     * Solve system of the form M_eff * acc = Q + Q_con
     * Q is the forces present in unconstrained system
     * Q_con is the generalized constraint forces */
    // calculate Q
    Eigen::VectorXd Q = Eigen::VectorXd::Zero(m_body_dof_7); // (7m, 1)

    // hydrodynamic force
    if (m_system->fluidDensity() > 0)
    {
        Q.noalias() += m_potHydro->fHydroNoInertia().segment(0, m_body_dof_7);
    }

    /* calculate Q_con = K (b - A * M_eff^{-1} * Q)
     * Linear constraint system: A * acc = b
     * Linear proportionality: K = M_eff^{1/2} * (A * M_eff^{-1/2})^{+};
     * + is Moore-Penrose inverse */

    const Eigen::MatrixXd M_eff = m_potHydro->mTotalBodyCoords().block(0, 0, m_body_dof_7, m_body_dof_7); // (7m, 7m)

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
}

void
RungeKutta4::momForceFree()
{
    for (int body_id = 0; body_id < m_body_dof; body_id++)
    {
        const int row_first{7 * body_id};
        const int row_last{};
    }

    // TODO: segment off by bodies

    // calculate M_tilde = Sigma * M_total * Sigma^T;
    Eigen::MatrixXd             M_tilde_hold = m_potHydro->mTotal() * m_system->rbmConn().transpose();
    Eigen::Matrix<double, 7, 7> M_tilde      = m_system->rbmConn() * M_tilde_hold;

    /* ANCHOR: Solve for rigid body motion velocity components */
    // calculate P_script = Sigma * M_total * V_articulation;
    Eigen::VectorXd             P_script_hold = m_potHydro->mTotal() * m_system->velocitiesParticlesArticulation();
    Eigen::Matrix<double, 7, 1> P_script      = m_system->rbmConn() * P_script_hold;

    // calculate U_swim = - M_tilde_inv * P_script; U_swim has translation and rotation components
    const Eigen::Matrix<double, 7, 1> U_swim = -M_tilde.fullPivLu().solve(P_script);

    /* ANCHOR: Output velocity data back to m_system */
    // calculate U = Sigma^T * U_swim + v_artic
    Eigen::Matrix<double, 7, 1> U_out = m_system->rbmConn().transpose() * U_swim;
    U_out.noalias() += m_system->velocitiesParticlesArticulation();
    m_system->setVelocitiesParticles(U_out);

    // // update hydrodynamic force terms for acceleration components
    // m_potHydro->updateForcesOnly();

    // /* ANCHOR: Solve for rigid body motion acceleration components */
    // // calculate b = a_artic + W * r + [v_artic ^ ]^T * Omega_C;  Omega_C = U_swim(3::5)
    // Eigen::Vector3d Omega_C = m_systemParam.U_swim.segment<3>(3);

    // Eigen::Matrix3d W = Omega_C * Omega_C.transpose();
    // W.noalias() -= Omega_C.squaredNorm() * m_system->i3();

    // Eigen::Matrix3d n_v_artic_cross;
    // Eigen::VectorXd b = m_system->accelerationsParticlesArticulation();

    // for (int i = 0; i < m_system->numParticles(); i++)
    // {
    //     int i3{3 * i};

    //     Eigen::Vector3d dr = m_system->positions().segment<3>(i3);
    //     dr.noalias() -= m_RLoc;
    //     b.segment<3>(i3).noalias() += W * dr;

    //     crossProdMat(-m_system->velocitiesParticlesArticulation().segment<3>(i3), n_v_artic_cross);
    //     b.segment<3>(i3).noalias() += n_v_artic_cross * Omega_C;
    // }

    // // calculate gMUU = \nabla M_added : U jdU
    // Eigen::VectorXd gMUU = m_potHydro->t2VelGrad();

    // // calculate F_script = Sigma * (M_total * b + gMUU)
    // Eigen::VectorXd F_script_hold = m_potHydro->mTotal() * b;
    // F_script_hold.noalias() += gMUU;
    // Eigen::VectorXd F_script = m_system->rbmConn() * F_script_hold;

    // // calculate A_swim = - M_tilde_inv * F_script; A_swim has translation and rotation
    // // components
    // m_systemParam.A_swim.noalias() = -M_tilde.fullPivLu().solve(F_script);

    // /* ANCHOR: Output acceleration data back to m_system */
    // // calculate A = Sigma^T * A_swim + b
    // Eigen::VectorXd A_out = m_system->rbmConn().transpose() * m_systemParam.A_swim;
    // A_out.noalias() += b;
    // m_system->setAccelerations(A_out);

    // /* ANCHOR: Output acceleration data back to input variable */
    // acc = A_out;

    // TODO: output body acceleration data too
}