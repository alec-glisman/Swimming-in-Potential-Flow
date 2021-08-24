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

    m_I_tilde       = m_I;
    m_I_tilde(2, 2) = -1.0;

    // Set specific variables to each system
    initializeSpecificVars();
    initializeConstraintLinearSystem();
}

rungeKutta4::~rungeKutta4()
{
    spdlog::get(m_logName)->info("Destructing Runge-Kutta 4th order");
}

void
rungeKutta4::integrate()
{
    double time{m_system->t() * m_system->tau()};

    /* ANCHOR: Solve system of form: y'(t) = f( y(t),  t )
     * @REFERENCE:
     * https://www.physicsforums.com/threads/using-runge-kutta-method-for-position-calc.553663/post-3634957
     */

    /* Step 1: k1 = f( y(t_0),  t_0 ),
     * initial conditions at current step */
    Eigen::VectorXd v1 = m_system->velocities();
    Eigen::VectorXd x1 = m_system->positions();
    Eigen::VectorXd a1 = Eigen::VectorXd::Zero(3 * m_system->numParticles());
    accelerationUpdate(a1, time);

    /* Step 2: k2 = f( y(t_0) + k1 * dt/2,  t_0 + dt/2 )
     * time rate-of-change k1 evaluated halfway through time step (midpoint) */
    Eigen::VectorXd v2 = v1;
    v2.noalias() += m_c1_2_dt * a1;
    Eigen::VectorXd x2 = x1;
    v2.noalias() += m_c1_2_dt * v2;

    m_system->setVelocities(v2);
    m_system->setPositions(x2);
    m_potHydro->update();

    Eigen::VectorXd a2 = Eigen::VectorXd::Zero(3 * m_system->numParticles());
    accelerationUpdate(a2, time + m_c1_2_dt);

    /* Step 3: k3 = f( y(t_0) + k2 * dt/2,  t_0 + dt/2 )
     * time rate-of-change k2 evaluated halfway through time step (midpoint) */
    Eigen::VectorXd v3 = v1;
    v3.noalias() += m_c1_2_dt * a2;
    Eigen::VectorXd x3 = x1;
    x3.noalias() += m_c1_2_dt * v3;

    m_system->setVelocities(v3);
    m_system->setPositions(x3);
    m_potHydro->update();

    Eigen::VectorXd a3 = Eigen::VectorXd::Zero(3 * m_system->numParticles());
    accelerationUpdate(a3, time + m_c1_2_dt);

    /* Step 4: k4 = f( y(t_0) + k3 * dt,  t_0 + dt )
     * time rate-of-change k3 evaluated at end of step (endpoint) */
    Eigen::VectorXd v4 = v1;
    v4.noalias() += m_dt * a3;
    Eigen::VectorXd x4 = x1;
    x4.noalias() += m_dt * v3;

    m_system->setVelocities(v4);
    m_system->setPositions(x4);
    m_potHydro->update();

    Eigen::VectorXd a4 = Eigen::VectorXd::Zero(3 * m_system->numParticles());
    accelerationUpdate(a4, time + m_dt);

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

    m_system->setVelocities(v_out);
    m_system->setPositions(x_out);
    m_potHydro->update();

    Eigen::VectorXd a_out = Eigen::VectorXd::Zero(3 * m_system->numParticles());
    accelerationUpdate(a_out, time + m_dt);
    m_system->setAccelerations(a_out);
}

/* REVIEW[epic=Change,order=1]: Change initializeSpecificVars() for different systems */
void
rungeKutta4::initializeSpecificVars()
{
    spdlog::get(m_logName)->info("Running initializeSpecificVars()");

    auto gsdParser = m_system->gsdUtil();

    /* ANCHOR: set specific parameters */

    // oscillation velocity amplitude
    spdlog::get(m_logName)->info("GSD parsing U0");
    double U0{-1.0};
    auto   return_bool = gsdParser->readChunk(&U0, gsdParser->frame(), "log/swimmer/U0", 8);
    m_system->setReturnBool(return_bool);
    m_system->check_gsd_return();
    m_systemParam.U0 = U0;
    spdlog::get(m_logName)->info("U_0 : {0}", U0);
    assert(m_systemParam.U0 == U0 && "U0 not properly set");

    // oscillation frequency
    spdlog::get(m_logName)->info("GSD parsing omega");
    double omega{-1.0};
    return_bool = gsdParser->readChunk(&omega, gsdParser->frame(), "log/swimmer/omega", 8);
    m_system->setReturnBool(return_bool);
    m_system->check_gsd_return();
    m_systemParam.omega = omega;
    spdlog::get(m_logName)->info("omega : {0}", omega);
    assert(m_systemParam.omega == omega && "omega not properly set");

    // phase shift between oscillators
    spdlog::get(m_logName)->info("GSD parsing phaseShift");
    double phaseShift{-1.0};
    return_bool =
        gsdParser->readChunk(&phaseShift, gsdParser->frame(), "log/swimmer/phase_shift", 8);
    m_system->setReturnBool(return_bool);
    m_system->check_gsd_return();
    m_systemParam.phaseShift = phaseShift;
    spdlog::get(m_logName)->info("phaseShift : {0}", phaseShift);
    assert(m_systemParam.phaseShift == phaseShift && "phaseShift not properly set");

    // average separation
    spdlog::get(m_logName)->info("GSD parsing Ravg");
    double RAvg{-1.0};
    return_bool = gsdParser->readChunk(&RAvg, gsdParser->frame(), "log/swimmer/R_avg", 8);
    m_system->setReturnBool(return_bool);
    m_system->check_gsd_return();
    m_systemParam.RAvg = RAvg;
    spdlog::get(m_logName)->info("RAvg : {0}", RAvg);
    assert(m_systemParam.RAvg == RAvg && "Ravg not properly set");

    /* ANCHOR: set initial conditions */
    spdlog::get(m_logName)->info("Setting initial conditions");

    // set articulation velocities
    spdlog::get(m_logName)->info("Calling articulationVel()");
    articulationVel(0.0);
    spdlog::get(m_logName)->info("Calling articulationAcc()");
    articulationAcc(0.0);
    // calculate locater point motion via linear and angular momentum free conditions
    spdlog::get(m_logName)->info("Calling rLoc()");
    rLoc();
    spdlog::get(m_logName)->info("Updating (for first time) hydrodynamic tensors");
    m_potHydro->update();
    spdlog::get(m_logName)->info("Calling momentumLinAngFree()");
    momentumLinAngFree();

    // write updated kinematics to original frame (appending current file)
    spdlog::get(m_logName)->info("Updating input GSD with updated kinematic initial conditions");
    spdlog::get(m_logName)->critical(
        "Frame 0 contains header information but incorrect kinematics");
    spdlog::get(m_logName)->critical(
        "Frame 1 should be treated as correct starting frame for analysis");
    gsdParser->writeFrame();
}

/* REVIEW[epic=Change,order=1]: Change initializeConstraintLinearSystem() for different systems */
void
rungeKutta4::initializeConstraintLinearSystem()
{
    // calculate A, not a function of time
    m_A = Eigen::MatrixXd::Zero(15, 18);

    m_A.block<3, 3>(0, 0).noalias() = m_I;
    m_A.block<3, 3>(6, 0).noalias() = m_I;

    m_A.block<3, 3>(0, 3).noalias()  = -m_I;
    m_A.block<3, 3>(3, 3).noalias()  = -m_I;
    m_A.block<3, 3>(12, 3).noalias() = m_I;

    m_A.block<3, 3>(3, 6).noalias()  = m_I;
    m_A.block<3, 3>(12, 6).noalias() = m_I;

    m_A.block<3, 3>(6, 9).noalias() = -m_I_tilde;

    m_A.block<3, 3>(9, 12).noalias() = -m_I_tilde;

    m_A.block<3, 3>(12, 15).noalias() = -m_I_tilde;
}

void
rungeKutta4::constraintLinearSystem(double dimensional_time)
{
    // update articulation velocities for constraint calculation
    articulationAcc(dimensional_time);

    // calculate b, function of time
    m_b                         = Eigen::VectorXd::Zero(15, 1);
    m_b.segment<3>(0).noalias() = m_accArtic.segment<3>(0);
    m_b.segment<3>(3).noalias() = m_accArtic.segment<3>(6);
}

/* REVIEW[epic=Change,order=2]: Change constraint linear system (A, b) for each system */
void
rungeKutta4::accelerationUpdate(Eigen::VectorXd& acc, double dimensional_time)
{
    /* NOTE: Following the formalism developed in Udwadia & Kalaba (1992) Proc. R. Soc. Lond. A
     * Solve system of the form M_total * acc = Q + Q_con
     * Q is the forces present in unconstrained system
     * Q_con is the generalized constraint forces */

    constraintLinearSystem(dimensional_time);

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
    Eigen::MatrixXd M_total_halfPower         = eigensolver.operatorSqrt();
    Eigen::MatrixXd M_total_negativeHalfPower = eigensolver.operatorInverseSqrt();

    // calculate K
    Eigen::MatrixXd AM_nHalf      = m_A * M_total_negativeHalfPower;
    Eigen::MatrixXd AM_nHalf_pInv = AM_nHalf.completeOrthogonalDecomposition().pseudoInverse();
    Eigen::MatrixXd K             = M_total_halfPower * AM_nHalf_pInv;

    // calculate Q_con
    Eigen::MatrixXd M_total_inv   = m_potHydro->mTotal().inverse();
    Eigen::MatrixXd M_total_invQ  = M_total_inv * Q;
    Eigen::VectorXd AM_total_invQ = m_A * M_total_invQ;
    Eigen::VectorXd b_tilde       = m_b;
    b_tilde.noalias() -= AM_total_invQ;
    Eigen::VectorXd Q_con = K * b_tilde;

    // calculate accelerations
    Eigen::VectorXd Q_total = Q;
    Q_total.noalias() += Q_con;
    acc.noalias() = m_potHydro->mTotal().llt().solve(Q_total);
}

void
rungeKutta4::momentumLinAngFree()
{
    /* ANCHOR: Solve for rigid body motion (rbm) tensors */
    // initialize variables
    Eigen::MatrixXd rbmconn = Eigen::MatrixXd::Zero(6, 3 * m_system->numParticles()); // [6 x 3N]

    // assemble the rigid body motion connectivity tensor (Sigma);  [6 x 3N]
    for (int i = 0; i < m_system->numParticles(); i++)
    {
        int i3{3 * i};

        Eigen::Vector3d n_dr = -m_system->positions().segment<3>(i3);
        n_dr.noalias() += m_RLoc;
        Eigen::Matrix3d n_dr_cross;
        crossProdMat(n_dr, n_dr_cross);

        rbmconn.block<3, 3>(0, i3).noalias() = m_I;        // translation-translation couple
        rbmconn.block<3, 3>(3, i3).noalias() = n_dr_cross; // translation-rotation couple
    }
    Eigen::MatrixXd rbmconn_T = rbmconn.transpose();

    // calculate M_tilde = Sigma * M_total * Sigma^T;  [6 x 6]
    Eigen::MatrixXd M_tilde_hold = m_potHydro->mTotal() * rbmconn_T;
    Eigen::MatrixXd M_tilde      = rbmconn * M_tilde_hold;

    /* ANCHOR: Solve for rigid body motion velocity components */
    // calculate P_script = Sigma * M_total * V_articulation;  [6 x 1]
    Eigen::VectorXd P_script_hold = m_potHydro->mTotal() * m_velArtic;
    Eigen::VectorXd P_script      = rbmconn * P_script_hold;

    // calculate U_swim = - M_tilde_inv * P_script; U_swim has translation and rotation
    // components
    Eigen::VectorXd U_swim = -M_tilde.fullPivLu().solve(P_script);

    /* ANCHOR: Output velocity data back to m_system */
    // calculate U = Sigma^T * U_swim + v_artic
    Eigen::VectorXd U_out = rbmconn_T * U_swim;
    U_out.noalias() += m_velArtic;
    m_system->setVelocities(U_out);

    // update hydrodynamic force terms for acceleration components
    m_potHydro->updateForcesOnly();

    /* ANCHOR: Solve for rigid body motion acceleration components */
    // calculate b = a_artic + W * r + [v_artic ^ ]^T * Omega_C;  Omega_C = U_swim(3::5)
    Eigen::Vector3d Omega_C = U_swim.segment<3>(3);

    Eigen::Matrix3d W = Omega_C * Omega_C.transpose();
    W.noalias() -= Omega_C.squaredNorm() * m_I;

    Eigen::Matrix3d n_v_artic_cross;
    Eigen::VectorXd b = m_accArtic;

    for (int i = 0; i < m_system->numParticles(); i++)
    {
        int i3{3 * i};

        Eigen::Vector3d dr = m_system->positions().segment<3>(i3);
        dr.noalias() -= m_RLoc;
        b.segment<3>(i3).noalias() += W * dr;

        crossProdMat(-m_velArtic.segment<3>(i3), n_v_artic_cross);
        b.segment<3>(i3).noalias() += n_v_artic_cross * Omega_C;
    }

    // calculate gMUU = \nabla M_added : U jdU
    Eigen::VectorXd gMUU = m_potHydro->t2VelGrad();

    // calculate F_script = Sigma * (M_total * b + gMUU)
    Eigen::VectorXd F_script_hold = m_potHydro->mTotal() * b;
    F_script_hold.noalias() += gMUU;
    Eigen::VectorXd F_script = rbmconn * F_script_hold;

    // calculate A_swim = - M_tilde_inv * F_script; A_swim has translation and rotation
    // components
    Eigen::VectorXd A_swim = -M_tilde.fullPivLu().solve(F_script);

    /* ANCHOR: Output acceleration data back to m_system */
    // calculate A = Sigma^T * A_swim + b
    Eigen::VectorXd A_out = rbmconn_T * A_swim;
    A_out.noalias() += b;
    m_system->setAccelerations(A_out);
}

/* REVIEW[epic=Change,order=3]: Change assignment of m_velArtic for different systems */
void
rungeKutta4::articulationVel(double dimensional_time)
{
    // ANCHOR: Orientation vectors, q = R_1 - R_3
    Eigen::Vector3d q = m_system->positions().segment<3>(0) - m_system->positions().segment<3>(6);
    q.normalize();
    Eigen::Vector3d q_tilde = m_I_tilde * q;

    // articulation velocity magnitudes
    double v1_mag = m_systemParam.U0 * cos(m_systemParam.omega * dimensional_time);
    double v3_mag =
        m_systemParam.U0 * cos(m_systemParam.omega * dimensional_time + m_systemParam.phaseShift);

    // Zero and then calculate m_velArtic
    m_velArtic                             = Eigen::VectorXd::Zero(3 * m_system->numParticles());
    m_velArtic.segment<3>(3 * 0).noalias() = v1_mag * q;
    m_velArtic.segment<3>(3 * 2).noalias() = v3_mag * q;
    m_velArtic.segment<3>(3 * 3).noalias() = v1_mag * q_tilde;
    m_velArtic.segment<3>(3 * 5).noalias() = v3_mag * q_tilde;
}

/* REVIEW[epic=Change,order=4]: Change assignment of m_accArtic for different systems */
void
rungeKutta4::articulationAcc(double dimensional_time)
{
    // ANCHOR: Orientation vectors, q = R_1 - R_3
    Eigen::Vector3d q = m_system->positions().segment<3>(0) - m_system->positions().segment<3>(6);
    q.normalize();
    Eigen::Vector3d q_tilde = m_I_tilde * q;

    // articulation acceleration magnitudes
    double a1_mag =
        -m_systemParam.U0 * m_systemParam.omega * sin(m_systemParam.omega * dimensional_time);
    double a3_mag = -m_systemParam.U0 * m_systemParam.omega *
                    sin(m_systemParam.omega * dimensional_time + m_systemParam.phaseShift);

    // Zero and then calculate m_accArtic
    m_accArtic                             = Eigen::VectorXd::Zero(3 * m_system->numParticles());
    m_accArtic.segment<3>(3 * 0).noalias() = a1_mag * q;
    m_accArtic.segment<3>(3 * 2).noalias() = a3_mag * q;
    m_accArtic.segment<3>(3 * 3).noalias() = a1_mag * q_tilde;
    m_accArtic.segment<3>(3 * 5).noalias() = a3_mag * q_tilde;

    m_accArtic = Eigen::VectorXd::Zero(3 * m_system->numParticles());
}

/* REVIEW[epic=Change,order=5]: Change assignment of m_RLoc for different systems */
void
rungeKutta4::rLoc()
{
    m_RLoc = m_system->positions().segment<3>(3 * 1);
}
