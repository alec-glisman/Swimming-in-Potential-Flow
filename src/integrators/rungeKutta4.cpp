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

    // // FIXME: debugging
    // std::setprecision(10);
    // Eigen::IOFormat CleanFmt(Eigen::FullPrecision, 0, ", ", ";\n", "", "", "[", "]");
    // std::string     sep = "\n----------------------------------------\n";
    // std::cout << "a_out at t = " << time + m_dt << std::endl;
    // std::cout << a_out.format(CleanFmt) << sep << std::endl;
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

    /* ANCHOR: Other member variables */
    int numRealPart{m_system->numParticles() / 2};
    int halfMat{3 * numRealPart};
    int fullMat{3 * m_system->numParticles()};

    spdlog::get(m_logName)->info("Setting sigma");
    m_systemParam.sigma = Eigen::MatrixXd::Zero(fullMat, halfMat);

    for (int i = 0; i < numRealPart; i++)
    {
        int i3{3 * i};

        m_systemParam.sigma.block<3, 3>(i3, i3).noalias()           = m_I;
        m_systemParam.sigma.block<3, 3>(i3 + halfMat, i3).noalias() = m_I_tilde;
    }

    m_systemParam.sigma_T = m_systemParam.sigma.transpose();

    /* ANCHOR: set initial conditions */
    spdlog::get(m_logName)->info("Setting initial conditions");

    // set articulation velocities
    spdlog::get(m_logName)->info("Updating (for first time) hydrodynamic tensors");
    m_potHydro->update();
    spdlog::get(m_logName)->info("Calling accelerationUpdate()");
    Eigen::VectorXd acc = m_system->accelerations();
    accelerationUpdate(acc, 0.0);

    /* ANCHOR: write updated kinematics to original frame (appending current file) */
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
    // Initialize matrices
    m_A = Eigen::MatrixXd::Zero(m_systemParam.num_constraints, m_systemParam.num_DoF);
    m_b = Eigen::VectorXd::Zero(m_systemParam.num_constraints);

    updateConstraintLinearSystem(m_system->t());
}

/* REVIEW[epic=Change,order=2]: Change constraint linear system (A, b) for each system */
void
rungeKutta4::updateConstraintLinearSystem(double dimensional_time)
{
    /* ANCHOR: Calculate quantities for linear system */
    // orientation vector, q = R_1 - R_3
    Eigen::Vector3d q = m_system->positions().segment<3>(3 * 0);
    q.noalias() -= m_system->positions().segment<3>(3 * 2);
    q.normalize();

    /* ANCHOR: Calculate A, function of time (Indexing: (constraint #, particle DoF #)) */
    m_A.setZero(m_systemParam.num_constraints, m_systemParam.num_DoF);

    // (1): q^T * U_1 - q^T * U_2 = || V_1 ||
    m_A.block<1, 3>(0, 0).noalias() = q;
    m_A.block<1, 3>(0, 3).noalias() = -q;

    // (2): q^T * U_3 - q^T * U_2 = || V_3 ||
    m_A.block<1, 3>(1, 6).noalias() = q;
    m_A.block<1, 3>(1, 3).noalias() = -q;

    // (3-5): U_1 - I_tilde * U_4 = 0
    m_A.block<3, 3>(2, 0).noalias() = m_I;
    m_A.block<3, 3>(2, 9).noalias() = -m_I_tilde;

    //  (6-8): U_2 - I_tilde * U_5 = 0
    m_A.block<3, 3>(5, 3).noalias()  = m_I;
    m_A.block<3, 3>(5, 12).noalias() = -m_I_tilde;

    //  (9-11): U_3 - I_tilde * U_6 = 0
    m_A.block<3, 3>(8, 6).noalias()  = m_I;
    m_A.block<3, 3>(8, 15).noalias() = -m_I_tilde;

    // (12): Collinear body constraint FIXME

    // relative positions
    Eigen::Vector3d rLoc = m_system->positions().segment<3>(3 * 1);
    Eigen::Vector3d r1   = m_system->positions().segment<3>(3 * 0);
    r1.noalias() -= rLoc;
    Eigen::Vector3d r3 = m_system->positions().segment<3>(3 * 2);
    r3.noalias() -= rLoc;

    // relative velocities
    Eigen::Vector3d vLoc = m_system->velocities().segment<3>(3 * 1);
    Eigen::Vector3d v1   = m_system->velocities().segment<3>(3 * 0);
    v1.noalias() -= vLoc;
    Eigen::Vector3d v3 = m_system->velocities().segment<3>(3 * 2);
    v3.noalias() -= vLoc;

    // Constants
    const double r1norm = r1.norm();
    const double r3norm = r3.norm();

    const double b1 = 1.0 / r1norm;
    const double b3 = 1.0 / r3norm;

    const double b1_sqr = 1.0 / std::pow(r1norm, 2.0);
    const double b1_cub = 1.0 / std::pow(r1norm, 3.0);
    const double b3_sqr = 1.0 / std::pow(r3norm, 2.0);
    const double b3_cub = 1.0 / std::pow(r3norm, 3.0);

    const double c1 = -b1_cub * r1.dot(v1);
    const double c3 = -b3_cub * r3.dot(v3);

    const Eigen::Matrix3d d1 = -b1_cub * r1 * r1.transpose() + b1 * m_I;
    const Eigen::Matrix3d d3 = -b3_cub * r3 * r3.transpose() + b3 * m_I;

    const Eigen::Vector3d e1 =
        2.0 * c1 * v1 - (3.0 * b1_sqr * c1 * r1.dot(v1) + b1_cub * v1.dot(v1)) * r1;
    const Eigen::Vector3d e3 =
        2.0 * c3 * v3 - (3.0 * b3_sqr * c3 * r3.dot(v3) + b3_cub * v3.dot(v3)) * r3;

    // kinematic vectors used in constraint
    const Eigen::Vector3d r1_hat = b1 * r1;
    const Eigen::Vector3d r3_hat = b3 * r3;

    const Eigen::Vector3d v1_hat = c1 * r1 + b1 * v1;
    const Eigen::Vector3d v3_hat = c3 * r3 + b3 * v3;

    // output quantities
    const double beta = r3_hat.dot(e1) + 2.0 * v1_hat.dot(v3_hat) + r1_hat.dot(e3);

    const Eigen::Vector3d alpha1 = r3_hat.transpose() * d1;
    const Eigen::Vector3d alpha3 = r1_hat.transpose() * d3;

    // output values
    m_A.block<1, 3>(11, 0).noalias() = alpha1;
    m_A.block<1, 3>(11, 3).noalias() = -alpha1;
    m_A.block<1, 3>(11, 3).noalias() -= alpha3;
    m_A.block<1, 3>(11, 6).noalias() = alpha3;

    /* ANCHOR: Calculate b, function of time */
    // articulation acceleration magnitudes
    double a1_mag =
        -m_systemParam.U0 * m_systemParam.omega * sin(m_systemParam.omega * dimensional_time);
    double a3_mag = -m_systemParam.U0 * m_systemParam.omega *
                    sin(m_systemParam.omega * dimensional_time + m_systemParam.phaseShift);

    // output results
    m_b.setZero(m_systemParam.num_constraints);
    m_b(0)  = a1_mag; // (1)
    m_b(1)  = a3_mag; // (2)
    m_b(11) = -beta;  // (12)
}

void
rungeKutta4::accelerationUpdate(Eigen::VectorXd& acc, double dimensional_time)
{
    /* REVIEW[epic=Change,seq=1]: Change integration schemes for different systems */
    bool useUdwadiaMethod{true};
    bool useImageMethod{true}; // can be used for initial conditions or total time integration

    if (useUdwadiaMethod && dimensional_time > 0.0)
    {
        udwadiaKalaba(acc, dimensional_time);
    }
    else if ((useImageMethod && !useUdwadiaMethod) || (useImageMethod && dimensional_time == 0.0))
    {
        momentumLinAngFreeImageSystem(acc, dimensional_time);
    }
    else
    {
        momentumLinAngFree(acc, dimensional_time);
    }
}

void
rungeKutta4::udwadiaKalaba(Eigen::VectorXd& acc, double dimensional_time)
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

    // update constraint linear system
    updateConstraintLinearSystem(dimensional_time);

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
rungeKutta4::momentumLinAngFree(Eigen::VectorXd& acc, double dimensional_time)
{
    /* ANCHOR: Compute articulation data */
    articulationAcc(dimensional_time);
    articulationVel(dimensional_time);
    rLoc();

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
    m_systemParam.U_swim.noalias() = -M_tilde.fullPivLu().solve(P_script);

    /* ANCHOR: Output velocity data back to m_system */
    // calculate U = Sigma^T * U_swim + v_artic
    Eigen::VectorXd U_out = rbmconn_T * m_systemParam.U_swim;
    U_out.noalias() += m_velArtic;
    m_system->setVelocities(U_out);

    // update hydrodynamic force terms for acceleration components
    m_potHydro->updateForcesOnly();

    /* ANCHOR: Solve for rigid body motion acceleration components */
    // calculate b = a_artic + W * r + [v_artic ^ ]^T * Omega_C;  Omega_C = U_swim(3::5)
    Eigen::Vector3d Omega_C = m_systemParam.U_swim.segment<3>(3);

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
    m_systemParam.A_swim.noalias() = -M_tilde.fullPivLu().solve(F_script);

    /* ANCHOR: Output acceleration data back to m_system */
    // calculate A = Sigma^T * A_swim + b
    Eigen::VectorXd A_out = rbmconn_T * m_systemParam.A_swim;
    A_out.noalias() += b;
    m_system->setAccelerations(A_out);

    /* ANCHOR: Output acceleration data back to input variable */
    acc = A_out;
}

void
rungeKutta4::momentumLinAngFreeImageSystem(Eigen::VectorXd& acc, double dimensional_time)
{
    /* ANCHOR: Method variables */
    int numRealPart{m_system->numParticles() / 2};

    /* ANCHOR: Compute articulation data */
    articulationAcc(dimensional_time);
    articulationVel(dimensional_time);
    rLoc();

    /* ANCHOR: Solve for rigid body motion (rbm) tensors */
    // initialize variables
    Eigen::MatrixXd rbmconn = Eigen::MatrixXd::Zero(6, 3 * m_system->numParticles()); // [6 x 3N]
    int             halfMat{3 * numRealPart};

    // assemble the rigid body motion connectivity tensor (Sigma);  [6 x 3N]
    for (int i = 0; i < numRealPart; i++)
    {
        int i3{3 * i};

        Eigen::Vector3d n_dr = -m_system->positions().segment<3>(i3);
        n_dr.noalias() += m_RLoc;
        Eigen::Matrix3d n_dr_cross;
        crossProdMat(n_dr, n_dr_cross);

        // "real" particles
        rbmconn.block<3, 3>(0, i3).noalias() = m_I;        // translation-translation couple
        rbmconn.block<3, 3>(3, i3).noalias() = n_dr_cross; // translation-rotation couple

        // "image" particles
        rbmconn.block<3, 3>(0, i3 + halfMat).noalias() = m_I; // translation-translation couple
        rbmconn.block<3, 3>(3, i3 + halfMat).noalias() =
            m_I_tilde * n_dr_cross; // translation-rotation couple
    }
    Eigen::MatrixXd rbmconn_T = rbmconn.transpose();

    // assemble rbm_conn for only "real" (not "image") particles
    Eigen::MatrixXd rbmconn_hat   = rbmconn.block(0, 0, 6, halfMat);
    Eigen::MatrixXd rbmconn_hat_T = rbmconn_hat.transpose();

    // calculate M_tilde = Sigma * M_total * sigma * Sigma_hat^T;  [6 x 6]
    Eigen::MatrixXd M_sigma      = m_potHydro->mTotal() * m_systemParam.sigma;
    Eigen::MatrixXd M_tilde_hold = M_sigma * rbmconn_hat_T;
    Eigen::MatrixXd M_tilde      = rbmconn * M_tilde_hold;

    /* ANCHOR: Solve for rigid body motion velocity components */
    // calculate P_script = Sigma * M_total * V;  [6 x 1]
    Eigen::VectorXd P_script_hold = m_potHydro->mTotal() * m_velArtic;
    Eigen::VectorXd P_script      = rbmconn * P_script_hold;

    // calculate U_swim = - M_tilde_inv * P_script;
    // U_swim has translation and rotation components
    m_systemParam.U_swim.noalias() = -M_tilde.fullPivLu().solve(P_script);

    /* ANCHOR: Output velocity data back to m_system */
    // calculate U = (sigma * Sigma_hat^T) * U_swim + V
    Eigen::VectorXd U_out_hold = rbmconn_hat_T * m_systemParam.U_swim;
    Eigen::VectorXd U_out      = m_systemParam.sigma * U_out_hold;
    U_out.noalias() += m_velArtic;
    m_system->setVelocities(U_out);

    // update hydrodynamic force terms for acceleration components
    m_potHydro->updateForcesOnly();

    /* ANCHOR: Solve for rigid body motion acceleration components */
    // calculate b_hat = [ (U_c - U_hat_alpha) ^ ] * Omega_C
    Eigen::Vector3d U_C     = m_systemParam.U_swim.segment<3>(0);
    Eigen::Vector3d Omega_C = m_systemParam.U_swim.segment<3>(3);

    Eigen::Matrix3d Uc_minus_Ualpha_cross;
    Eigen::VectorXd b_hat = Eigen::VectorXd::Zero(halfMat);

    for (int i = 0; i < numRealPart; i++)
    {
        int i3{3 * i};

        Eigen::Vector3d Uc_minus_Ualpha = U_C;
        Uc_minus_Ualpha.noalias() -= m_system->velocities().segment<3>(i3);

        crossProdMat(Uc_minus_Ualpha, Uc_minus_Ualpha_cross);
        b_hat.segment<3>(i3).noalias() += Uc_minus_Ualpha_cross * Omega_C;
    }

    // calculate c = sigma * b_hat + A_artic
    Eigen::VectorXd c = m_systemParam.sigma * b_hat;
    c.noalias() += m_accArtic;

    // calculate gMUU = \nabla M_added : ( U U )
    Eigen::VectorXd gMUU = m_potHydro->t2VelGrad();

    // FIXME
    // calculate F_script = Sigma * (M_total * c + gMUU)
    Eigen::VectorXd F_script_hold = m_potHydro->mTotal() * c;
    F_script_hold.noalias() += gMUU;
    Eigen::VectorXd F_script = rbmconn * F_script_hold;

    // calculate A_swim = - M_tilde_inv * F_script;
    // A_swim has translation and rotation components
    m_systemParam.A_swim.noalias() = -M_tilde.fullPivLu().solve(F_script);

    /* ANCHOR: Output acceleration data back to m_system */
    // calculate A = sigma * ( Sigma_hat^T * A_swim + b_hat )
    Eigen::VectorXd A_out_hold = rbmconn_hat_T * m_systemParam.A_swim;
    A_out_hold.noalias() += b_hat;
    Eigen::VectorXd A_out = m_systemParam.sigma * A_out_hold;
    m_system->setAccelerations(A_out);

    /* ANCHOR: Output acceleration data back to input variable */
    acc = A_out;
}

/* ANCHOR: Output data to GSD */
void
rungeKutta4::GSDOutput()
{
    spdlog::get(m_logName)->info("Writing momentumLinAngFree() data to GSD");

    std::vector<double> d_U_swim(6);
    d_U_swim.reserve(1); //! make sure we allocate
    for (uint32_t i = 0; i < 6; i++)
    {
        d_U_swim[i] = m_systemParam.U_swim(i);
    }
    spdlog::get(m_logName)->info("Writing swimmer/U_swim data to GSD");
    spdlog::get(m_logName)->info("[{0}, {1}, {2}, {3}, {4}, {5}]", d_U_swim[0], d_U_swim[1],
                                 d_U_swim[2], d_U_swim[3], d_U_swim[4], d_U_swim[5]);
    auto return_val = gsd_write_chunk(m_system->handle().get(), "swimmer/U_swim", GSD_TYPE_DOUBLE,
                                      6, 1, 0, (void*)&d_U_swim[0]);
    m_system->setReturnVal(return_val);
    m_system->check_gsd_return();

    std::vector<double> d_A_swim(6);
    d_A_swim.reserve(1); //! make sure we allocate
    for (uint32_t i = 0; i < 6; i++)
    {
        d_A_swim[i] = m_systemParam.A_swim(i);
    }
    spdlog::get(m_logName)->info("Writing swimmer/A_swim data to GSD");
    spdlog::get(m_logName)->info("[{0}, {1}, {2}, {3}, {4}, {5}]", d_A_swim[0], d_A_swim[1],
                                 d_A_swim[2], d_A_swim[3], d_A_swim[4], d_A_swim[5]);
    return_val = gsd_write_chunk(m_system->handle().get(), "swimmer/A_swim", GSD_TYPE_DOUBLE, 6, 1,
                                 0, (void*)&d_A_swim[0]);
    m_system->setReturnVal(return_val);
    m_system->check_gsd_return();

    spdlog::get(m_logName)->flush();
}

/* REVIEW[epic=Change,order=3]: Change assignment of m_velArtic for different systems */
void
rungeKutta4::articulationVel(double dimensional_time)
{
    // ANCHOR: Orientation vectors, q = R_1 - R_3
    Eigen::Vector3d q =
        m_system->positions().segment<3>(3 * 0) - m_system->positions().segment<3>(3 * 2);
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
    Eigen::Vector3d q =
        m_system->positions().segment<3>(3 * 0) - m_system->positions().segment<3>(3 * 2);
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
}

/* REVIEW[epic=Change,order=5]: Change assignment of m_RLoc for different systems */
void
rungeKutta4::rLoc()
{
    m_RLoc = m_system->positions().segment<3>(3 * 1);
}
