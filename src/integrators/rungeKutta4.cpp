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

    // Set specific variables to each system
    initializeSpecificVars();
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
    Eigen::VectorXd v1 = m_system->velocities;
    Eigen::VectorXd x1 = m_system->positions;

    m_potHydro->update();

    Eigen::VectorXd a1 = Eigen::VectorXd::Zero(3 * m_system->numParticles());
    accelerationUpdate(a1, time);

    /* Step 2: k2 = dt * f ( y(t_0) + 1/2 * k1, t_0 + 1/2 dt ) */
    Eigen::VectorXd v2 = v1 + m_c1_2_dt * a1;
    Eigen::VectorXd x2 = x1 + m_c1_2_dt * v2;

    m_system->accelerations.noalias() = a1;
    m_system->velocities.noalias()    = v2;
    m_system->positions.noalias()     = x2;

    m_potHydro->update();

    Eigen::VectorXd a2 = Eigen::VectorXd::Zero(3 * m_system->numParticles());
    accelerationUpdate(a2, time + m_c1_2_dt);

    /* Step 3 */
    Eigen::VectorXd v3 = v1 + m_c1_2_dt * a2;
    Eigen::VectorXd x3 = x1 + m_c1_2_dt * v3;

    m_system->accelerations.noalias() = a2;
    m_system->velocities.noalias()    = v3;
    m_system->positions.noalias()     = x3;

    m_potHydro->update();

    Eigen::VectorXd a3 = Eigen::VectorXd::Zero(3 * m_system->numParticles());
    accelerationUpdate(a3, time + m_c1_2_dt);

    /* Step 4 */
    Eigen::VectorXd v4 = v1 + m_dt * a3;
    Eigen::VectorXd x4 = x1 + m_dt * v3;

    m_system->accelerations.noalias() = a3;
    m_system->velocities.noalias()    = v4;
    m_system->positions.noalias()     = x4;

    m_potHydro->update();

    Eigen::VectorXd a4 = Eigen::VectorXd::Zero(3 * m_system->numParticles());
    accelerationUpdate(a4, time + m_dt);

    /* Output data */
    m_system->positions.noalias()  = x1 + m_c1_6_dt * (v1 + 2.0 * v2 + 2.0 * v3 + v4);
    m_system->velocities.noalias() = v1 + m_c1_6_dt * (a1 + 2.0 * a2 + 2.0 * a3 + a4);

    m_potHydro->update();

    Eigen::VectorXd a_out = Eigen::VectorXd::Zero(3 * m_system->numParticles());
    accelerationUpdate(a_out, time + m_dt);
    m_system->accelerations.noalias() = a_out;
}

void
rungeKutta4::accelerationUpdate(Eigen::VectorXd& acc, double dimensional_time)
{
    /* NOTE: Following the formalism developed in Udwadia & Kalaba (1992) Proc. R. Soc. Lond. A
     * Solve system of the form M_total * acc = Q + Q_con
     * Q is the forces present in unconstrained system
     * Q_con is the generalized constraint forces */

    // calculate Q
    Eigen::VectorXd Q = Eigen::VectorXd::Zero(3 * m_system->numParticles());
    if (m_system->particleDensity() >= 0) // hydrodynamic force
    {
        Q.noalias() += m_potHydro->fHydroNoInertia();
    }

    /* calculate Q_con = K (b - A * M_total^{-1} * Q)
     * Linear constraint system: A * acc = b
     * Linear proportionality: K = M_total^{1/2} * (A * M_total^{-1/2})^{+};
     * + is Moore-Penrose inverse */

    // REVIEW[epic=Change,order=5]: alter constraint linear system for each system
    // calculate A
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(2, 9);
    A(0, 0)           = 1.0;
    A(0, 3)           = -1.0;
    A(1, 3)           = -1.0;
    A(1, 6)           = 1.0;
    // calculate B
    Eigen::Vector2d b = Eigen::Vector2d::Zero(2, 1);
    b(0)              = -m_U0 * m_omega * sin(m_omega * dimensional_time);
    b(1)              = -m_U0 * m_omega * sin(m_omega * dimensional_time + m_phase_shift);

    // calculate M^{1/2} & M^{-1/2}
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(m_potHydro->mTotal());
    if (eigensolver.info() != Eigen::Success)
    {
        spdlog::get(m_logName)->error("Computing eigendecomposition of M_total failed at t={0}",
                                      m_system->t());
        throw std::runtime_error("Computing eigendecomposition of M_total failed");
    }
    Eigen::MatrixXd M_total_halfPower    = eigensolver.operatorSqrt();
    Eigen::MatrixXd M_total_negativeHalf = eigensolver.operatorInverseSqrt();

    // calculate K
    Eigen::MatrixXd AM_nHalf      = A * M_total_negativeHalf;
    Eigen::MatrixXd AM_nHalf_pInv = AM_nHalf.completeOrthogonalDecomposition().pseudoInverse();
    Eigen::MatrixXd K             = M_total_halfPower * AM_nHalf_pInv;

    // calculate Q_con
    Eigen::MatrixXd M_total_inv   = m_potHydro->mTotal().inverse();
    Eigen::MatrixXd M_total_invQ  = M_total_inv * Q;
    Eigen::VectorXd AM_total_invQ = A * M_total_invQ;
    Eigen::VectorXd b_tilde       = b;
    b_tilde.noalias() -= AM_total_invQ;
    Eigen::VectorXd Q_con = K * b_tilde;

    // calculate accelerations
    Eigen::VectorXd Q_total = Q;
    Q_total.noalias() += Q_con;
    acc.noalias() = m_potHydro->mTotal().llt().solve(Q_total);
}

/* REVIEW[epic=Change,order=1]: Change initializeSpecificVars() for different systems*/
void
rungeKutta4::initializeSpecificVars()
{
    spdlog::get(m_logName)->info("Running initializeSpecificVars()");

    auto gsdParser = m_system->gsdUtil();

    /* ANCHOR: set specific parameters */

    // oscillation velocity amplitude
    spdlog::get(m_logName)->info("GSD parsing U0");
    float U0{-1.0};
    auto  return_bool = gsdParser->readChunk(&U0, gsdParser->frame(), "log/swimmer/U0", 4);
    m_system->setReturnBool(return_bool);
    m_system->check_gsd_return();
    m_U0 = double(U0);
    spdlog::get(m_logName)->info("U_0 : {0}", m_U0);
    assert(double(U0) == m_U0 && "U0 not properly set");

    // oscillation frequency
    spdlog::get(m_logName)->info("GSD parsing omega");
    float omega{-1.0};
    return_bool = gsdParser->readChunk(&omega, gsdParser->frame(), "log/swimmer/omega", 4);
    m_system->setReturnBool(return_bool);
    m_system->check_gsd_return();
    m_omega = double(omega);
    spdlog::get(m_logName)->info("omega : {0}", m_omega);
    assert(double(omega) == m_omega && "omega not properly set");

    // phase shift between oscillators
    spdlog::get(m_logName)->info("GSD parsing phase_shift");
    float phase_shift{-1.0};
    return_bool =
        gsdParser->readChunk(&phase_shift, gsdParser->frame(), "log/swimmer/phase_shift", 4);
    m_system->setReturnBool(return_bool);
    m_system->check_gsd_return();
    m_phase_shift = double(phase_shift);
    spdlog::get(m_logName)->info("phase_shift : {0}", m_phase_shift);
    assert(double(phase_shift) == m_phase_shift && "phase_shift not properly set");

    // average separation
    spdlog::get(m_logName)->info("GSD parsing R_avg");
    float R_avg{-1.0};
    return_bool = gsdParser->readChunk(&R_avg, gsdParser->frame(), "log/swimmer/R_avg", 4);
    m_system->setReturnBool(return_bool);
    m_system->check_gsd_return();
    m_R_avg = double(R_avg);
    spdlog::get(m_logName)->info("R_avg : {0}", m_R_avg);
    assert(double(R_avg) == m_R_avg && "R_avg not properly set");

    /* ANCHOR: set initial conditions */
    spdlog::get(m_logName)->info("Setting initial conditions");

    // reset initial positions
    spdlog::get(m_logName)->warn("Setting initial positions");
    m_system->positions.setZero();
    m_system->positions(0) = m_R_avg;
    m_system->positions(6) = -m_R_avg + m_U0 / m_omega * sin(m_phase_shift);

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

/* REVIEW[epic=Change,order=2]: Change articulationVel() for different systems*/
void
rungeKutta4::articulationVel(double dimensional_time)
{
    m_velArtic = Eigen::VectorXd::Zero(3 * m_system->numParticles());

    m_velArtic(0) = m_U0 * cos(m_omega * dimensional_time);
    m_velArtic(6) = m_U0 * cos(m_omega * dimensional_time + m_phase_shift);
}

/* REVIEW[epic=Change,order=3]: Change articulationAcc() for different systems*/
void
rungeKutta4::articulationAcc(double dimensional_time)
{
    m_accArtic = Eigen::VectorXd::Zero(3 * m_system->numParticles());

    m_accArtic(0) = -m_U0 * m_omega * sin(m_omega * dimensional_time);
    m_accArtic(6) = -m_U0 * m_omega * sin(m_omega * dimensional_time + m_phase_shift);
}

/* REVIEW[epic=Change,order=4]: Change rLoc() for different systems*/
void
rungeKutta4::rLoc()
{
    m_RLoc = m_system->positions.segment<3>(3);
}

void
rungeKutta4::momentumLinAngFree()
{
    /* ANCHOR: Solve for rigid body motion (rbm) tensors */
    // initialize variables
    Eigen::Matrix3d I       = Eigen::Matrix3d::Identity(3, 3);                        // [3 x 3]
    Eigen::MatrixXd rbmconn = Eigen::MatrixXd::Zero(6, 3 * m_system->numParticles()); // [6 x 3N]

    // assemble the rigid body motion connectivity tensor (Sigma);  [6 x 3N]
    for (int i = 0; i < m_system->numParticles(); i++)
    {
        int i3{3 * i};

        Eigen::Vector3d n_dr = -m_system->positions.segment<3>(i3);
        n_dr.noalias() += m_RLoc;
        Eigen::Matrix3d n_dr_cross;
        crossProdMat(n_dr, n_dr_cross);

        rbmconn.block<3, 3>(0, i3).noalias() = I;          // translation-translation couple
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
    m_system->velocities.noalias() = rbmconn_T * U_swim;
    m_system->velocities.noalias() += m_velArtic;

    // update hydrodynamic force terms for acceleration components
    m_potHydro->updateForcesOnly();

    /* ANCHOR: Solve for rigid body motion acceleration components */
    // calculate b = a_artic + W * r + [v_artic ^ ]^T * Omega_C;  Omega_C = U_swim(3::5)
    Eigen::Vector3d Omega_C = U_swim.segment<3>(3);

    Eigen::Matrix3d W = Omega_C * Omega_C.transpose();
    W.noalias() -= Omega_C.squaredNorm() * I;

    Eigen::Matrix3d n_v_artic_cross;
    Eigen::VectorXd b = m_accArtic;

    for (int i = 0; i < m_system->numParticles(); i++)
    {
        int i3{3 * i};

        Eigen::Vector3d dr = m_system->positions.segment<3>(i3);
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
    m_system->accelerations.noalias() = rbmconn_T * A_swim;
    m_system->accelerations.noalias() += b;
}