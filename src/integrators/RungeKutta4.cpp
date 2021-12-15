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

    // for initial conditions of all locater points, use the PF-free algorithm
    spdlog::get(m_logName)->critical("Setting initial conditions using PF-free algorithm.");

    /// @review_swimmer: internal dynamics turned off
    const bool internal_dyn_off{abs(m_system->sysSpecU0()) < 1e-12};

    if (internal_dyn_off)
    {
        spdlog::get(m_logName)->warn("Setting initial conditions using constant as internal dynamics off");
        Eigen::VectorXd          vel_body = m_system->velocitiesBodies();
        const double             vel_init_mag{1.7040237270147573e-07}; // from isolated swimmer with R_avg = 4.0
        const Eigen::Quaterniond quat_vel_body(0.0, vel_init_mag, 0.0, 0.0);

        for (int body_id = 0; body_id < m_system->numBodies(); body_id++)
        {
            const int body_id_7{7 * body_id};

            const Eigen::Quaterniond quat_body(
                m_system->positionsBodies()(body_id_7 + 3), m_system->positionsBodies()(body_id_7 + 4),
                m_system->positionsBodies()(body_id_7 + 5), m_system->positionsBodies()(body_id_7 + 6));

            const Eigen::Quaterniond rot_quat_vel_body = quat_body * quat_vel_body * quat_body.inverse();

            // rotate vel mag by unit quaternion
            vel_body.segment<3>(body_id_7).noalias() = rot_quat_vel_body.vec();

            // account for image system
            if ((m_system->imageSystem()) && (body_id >= m_system->numBodies() / 2))
            {
                vel_body(body_id_7 + 2) *= -1.0;
            }
        }

        m_system->setVelocitiesBodies(vel_body);
    }

    spdlog::get(m_logName)->critical("Updating SystemData and PotentialHydrodynamics classes for initial conditions.");
    m_system->update(single_core_device);   // update system kinematics and rbm tensors
    m_potHydro->update(single_core_device); // update hydrodynamic tensors

    if (!internal_dyn_off)
    {
        momForceFree(single_core_device); // set initial body kinematics
    }

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

    spdlog::get(m_logName)->critical("Constructor complete");
    spdlog::get(m_logName)->flush();
}

RungeKutta4::~RungeKutta4()
{
    spdlog::get(m_logName)->critical("Destructing Runge-Kutta 4th order");
    spdlog::get(m_logName)->flush();
    spdlog::drop(m_logName);
}

void
RungeKutta4::integrate(const Eigen::ThreadPoolDevice& device)
{
    integrateSecondOrder(device); // Udwadia-Kalaba method only gives acceleration components
}

void
RungeKutta4::integrateSecondOrder(const Eigen::ThreadPoolDevice& device)
{
    /* Step 1: k1 = f( y(t_0),  t_0 )
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
                                const Eigen::ThreadPoolDevice& device)
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
        Eigen::VectorXd acc_real_body = acc.segment(0, m_body_dof_7);
        udwadiaKalaba(acc_real_body);

        // update acceleration components using constraints
        acc.segment(0, m_body_dof_7).noalias() = acc_real_body;
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
    pos.segment(m_body_dof_7, m_body_dof_7) = pos.segment(0, m_body_dof_7);
    vel.segment(m_body_dof_7, m_body_dof_7) = vel.segment(0, m_body_dof_7);

    for (int img_body_id = m_body_dof; img_body_id < m_system->numBodies(); img_body_id++)
    {
        // `Eigen::Vector` indices
        const int img_body_id_7{7 * img_body_id};

        /* ANCHOR: Position image system: mirror image about xy-plane (transform z --> -z) */
        // (linear components) flip z component, leave x-y unchanged
        pos(img_body_id_7 + 2) *= -1;
        // (quaternion components) flip x-y components of quaternion vector part (see note in header file: look for "C1
        // and C2 are mirrored along the z-axis")
        pos.segment<2>(img_body_id_7 + 4) *= -1;

        /* ANCHOR: Velocity image system: mirror image about xy-plane (transform z --> -z)  */
        vel(img_body_id_7 + 2) *= -1;
        vel.segment<2>(img_body_id_7 + 4) *= -1;
    }
}

void
RungeKutta4::imageBodyAcc(Eigen::VectorXd& acc)
{
    acc.segment(m_body_dof_7, m_body_dof_7) = acc.segment(0, m_body_dof_7);

    for (int img_body_id = m_body_dof; img_body_id < m_system->numBodies(); img_body_id++)
    {
        // `Eigen::Vector` indices
        const int img_body_id_7{7 * img_body_id};

        /* ANCHOR: Acceleration image system: mirror image about xy-plane (transform z --> -z)   */
        acc(img_body_id_7 + 2) *= -1;
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
RungeKutta4::momForceFree(const Eigen::ThreadPoolDevice& device)
{
    const int num_particles{3};                   // Nb
    const int num_particles_3{3 * num_particles}; // 3 * Nb

    const Eigen::array<Eigen::Index, 3> extents_333 = {3, 3, 3};

    /* ANCHOR: Assemble the tensor quantities */
    Eigen::Matrix<double, 3, -1> rbmconn = Eigen::MatrixXd::Zero(3, num_particles_3); // (3 x 3 Nb)
    Eigen::Matrix3d              Id      = Eigen::MatrixXd::Identity(3, 3);

    Eigen::MatrixXd M_eff = Eigen::MatrixXd::Zero(num_particles_3, num_particles_3); // (3 Nb x 3 Nb)

    Eigen::Tensor<double, 3> grad_M_eff =
        Eigen::Tensor<double, 3>(num_particles_3, num_particles_3, num_particles_3); // (3 Nb x 3 Nb x 3 Nb)
    grad_M_eff.setZero();

    Eigen::VectorXd vel_artic = Eigen::VectorXd::Zero(num_particles_3); // (3 Nb x 1)
    Eigen::VectorXd acc_artic = Eigen::VectorXd::Zero(num_particles_3); // (3 Nb x 1)

    for (int particle_id = 0; particle_id < num_particles; particle_id++)
    {
        const int particle_id_3{3 * particle_id};
        const int particle_id_7{7 * particle_id};

        rbmconn.block<3, 3>(0, particle_id_3).noalias() = Id; // translation-translation couple

        for (int j = 0; j < num_particles; j++)
        {
            M_eff.block<3, 3>(particle_id_3, 3 * j).noalias() =
                m_potHydro->mTotal().block<3, 3>(particle_id_7, 7 * j); // translation-translation couple

            for (int k = 0; k < (0 + num_particles); k++)
            {
                const Eigen::array<Eigen::Index, 3> offsets_3   = {particle_id_3, 3 * j, 3 * k};
                const Eigen::array<Eigen::Index, 3> offsets_773 = {particle_id_7, 7 * j, 3 * k};

                grad_M_eff.slice(offsets_3, extents_333).device(device) =
                    m_potHydro->gradMAdded().slice(offsets_773, extents_333);
            }
        }

        vel_artic.segment<3>(particle_id_3).noalias() =
            m_system->velocitiesParticlesArticulation().segment<3>(particle_id_7); // linear components

        acc_artic.segment<3>(particle_id_3).noalias() =
            m_system->accelerationsParticlesArticulation().segment<3>(particle_id_7); // linear components
    }

    const Eigen::Matrix<double, -1, 3> rbmconn_T = rbmconn.transpose();

    /* ANCHOR: Calculate velocity components */
    // calculate M_tilde = Sigma * M_total * Sigma^T;  (3 x 3)
    const Eigen::Matrix<double, -1, 3> M_tilde_hold = M_eff * rbmconn_T;
    const Eigen::Matrix3d              M_tilde      = rbmconn * M_tilde_hold;

    /* STUB: Solve for rigid body motion velocity components */
    // calculate P_script = Sigma * M_total * V_articulation;  (3 x 1)
    const Eigen::VectorXd P_script_hold = M_eff * vel_artic;
    const Eigen::Vector3d P_script      = rbmconn * P_script_hold;

    // calculate U_swim = - M_tilde_inv * P_script; U_swim has translation components
    const Eigen::Vector3d U_swim = -M_tilde.fullPivLu().solve(P_script); // linear components
    // convert U_swim to 4d vector for quaternion rotation
    const Eigen::Quaterniond U_swim_4d(0.0, U_swim(0), U_swim(1), U_swim(2));

    /* ANCHOR: Output velocity data back to m_system */
    Eigen::VectorXd vel_body = Eigen::VectorXd::Zero(m_7M);

    for (int body_id = 0; body_id < (m_body_dof); body_id++)
    {
        const int body_id_7{7 * body_id};

        // body unit quaternion
        const Eigen::Vector4d    theta = m_system->positionsBodies().segment<4>(body_id_7 + 3);
        const Eigen::Quaterniond theta_body(theta(0), theta(1), theta(2), theta(3));

        const Eigen::Quaterniond U_swim_rot_4d = theta_body * U_swim_4d * theta_body.inverse();
        const Eigen::Vector3d    U_swim_rot_3d = U_swim_rot_4d.vec();

        vel_body.segment<3>(body_id_7).noalias() = U_swim_rot_3d;
    }

    if (m_system->imageSystem())
    {
        Eigen::VectorXd pos_body = m_system->positionsBodies();
        imageBodyPosVel(pos_body, vel_body);
    }

    m_system->setVelocitiesBodies(vel_body);

    /* ANCHOR: Calculate acceleration components */
    m_system->update(device); // update velocity components with locater motion

    Eigen::VectorXd vel_part = Eigen::VectorXd::Zero(num_particles_3); // (3 Nb x 1)

    for (int particle_id = 0; particle_id < (0 + num_particles); particle_id++)
    {
        const int particle_id_3{3 * particle_id};
        const int particle_id_7{7 * particle_id};

        vel_part.segment<3>(particle_id_3).noalias() =
            m_system->velocitiesParticles().segment<3>(particle_id_7); // linear components
    }

    // convert particle velocities to `Eigen::Tensor`
    Eigen::Tensor<double, 1> tens_vel_part = Eigen::Tensor<double, 1>(num_particles_3);
    tens_vel_part.device(device)           = TensorCast(vel_part, num_particles_3);

    // `Eigen::Tensor` contract indices
    const Eigen::array<Eigen::IndexPair<int>, 1> contract_ijl_l = {
        Eigen::IndexPair<int>(2, 0)}; // {i, j, l} . {l} --> {i, j}

    // calculate gradM U
    Eigen::Tensor<double, 2> tens_gradM_U = Eigen::Tensor<double, 2>(num_particles_3, num_particles_3);
    tens_gradM_U.device(device)           = grad_M_eff.contract(tens_vel_part, contract_ijl_l);
    const Eigen::MatrixXd gradM_U         = MatrixCast(tens_gradM_U, num_particles_3, num_particles_3, device);

    /* STUB: Solve for rigid body motion velocity components */
    // calculate F_script = Sigma * (M_total * A_articulation + (gradM U) U);  (3 x 1)
    Eigen::VectorXd F_script_hold = M_eff * acc_artic;
    F_script_hold.noalias() += gradM_U * vel_part;
    const Eigen::Vector3d F_script = rbmconn * F_script_hold;

    // calculate A_swim = - M_tilde_inv * P_script; A_swim has translation components
    const Eigen::Vector3d A_swim = -M_tilde.fullPivLu().solve(F_script); // linear components

    // convert A_swim to 4d vector for quaternion rotation
    Eigen::Quaterniond A_swim_4d(0.0, A_swim(0), A_swim(1), A_swim(2));

    /* ANCHOR: Output velocity data back to m_system */
    Eigen::VectorXd acc_body = Eigen::VectorXd::Zero(m_7M);

    for (int body_id = 0; body_id < (m_body_dof); body_id++)
    {
        const int body_id_7{7 * body_id};

        // body unit quaternion
        const Eigen::Vector4d    theta = m_system->positionsBodies().segment<4>(body_id_7 + 3);
        const Eigen::Quaterniond theta_body(theta(0), theta(1), theta(2), theta(3));

        const Eigen::Quaterniond A_swim_rot_4d = theta_body * A_swim_4d * theta_body.inverse();
        const Eigen::Vector3d    A_swim_rot_3d = A_swim_rot_4d.vec();

        acc_body.segment<3>(body_id_7).noalias() = A_swim_rot_3d;
    }

    if (m_system->imageSystem())
    {
        imageBodyAcc(acc_body);
    }

    m_system->setAccelerationsBodies(acc_body);
    m_system->update(device); // update acceleration components with locater motion
}