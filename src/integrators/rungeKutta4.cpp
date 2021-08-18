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
    /* ANCHOR: Solve system of form: y'(t) = f( y(t),  t ) */

    /* Step 1: k1 = f( y(t_0),  t_0 ),
     * initial conditions at current step */
    Eigen::VectorXd x1 = m_system->positions;
    Eigen::VectorXd v1 = m_system->velocities;
    Eigen::VectorXd a1 = m_system->accelerations;

    /* Step 2: k2 = f( y(t_0) + k1 * dt/2,  t_0 + dt/2 )
     * time rate-of-change k1 evaluated halfway through time step (midpoint) */
    Eigen::VectorXd x2 = m_c1_2_dt * v1;
    x2.noalias() += x1;
    m_system->positions.noalias() =
        x2; // temporarily change positions to get correct hydro tensors during acceleration update

    Eigen::VectorXd v2 = m_c1_2_dt * a1;
    v2.noalias() += v1;

    m_potHydro->update(); // update hydrodynamics tensors
    Eigen::VectorXd a2 = Eigen::VectorXd::Zero(3 * m_system->numParticles());
    accelerationUpdate(a2);

    /* Step 3: k3 = f( y(t_0) + k2 * dt/2,  t_0 + dt/2 )
     * time rate-of-change k2 evaluated halfway through time step (midpoint) */
    Eigen::VectorXd x3 = m_c1_2_dt * v2;
    x3.noalias() += x1;
    m_system->positions.noalias() = x3;

    Eigen::VectorXd v3 = m_c1_2_dt * a2;
    v3.noalias() += v1;

    m_potHydro->update();
    Eigen::VectorXd a3 = Eigen::VectorXd::Zero(3 * m_system->numParticles());
    accelerationUpdate(a3);

    /* Step 4: k4 = f( y(t_0) + k3 * dt,  t_0 + dt )
     * time rate-of-change k3 evaluated at end of step (endpoint) */
    Eigen::VectorXd x4 = m_dt * v3;
    x4.noalias() += x1;
    m_system->positions.noalias() = x4;

    Eigen::VectorXd v4 = m_dt * a3;
    v4.noalias() += v1;

    m_potHydro->update();
    Eigen::VectorXd a4 = Eigen::VectorXd::Zero(3 * m_system->numParticles());
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
    Eigen::VectorXd f = Eigen::VectorXd::Zero(3 * m_system->numParticles());

    if (m_system->particleDensity() >= 0) // hydrodynamic force
    {
        f.noalias() += m_potHydro->fHydroNoInertia();
    }

    acc.noalias() = m_potHydro->mTotal().llt().solve(f);
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

    /* ANCHOR: set initial conditions */
    spdlog::get(m_logName)->info("Setting initial conditions");

    // set articulation velocities
    spdlog::get(m_logName)->info("Calling articulationVel()");
    Eigen::VectorXd v_artic = articulationVel(0.0);
    spdlog::get(m_logName)->info("Calling articulationAcc()");
    Eigen::VectorXd a_artic = articulationAcc(0.0);
    // calculate locater point motion via linear and angular momentum free conditions
    spdlog::get(m_logName)->info("Calling rLoc()");
    Eigen::Vector3d rloc = rLoc();
    spdlog::get(m_logName)->info("Updating (for first time) hydrodynamic tensors");
    m_potHydro->update();
    spdlog::get(m_logName)->info("Calling momentumLinAngFree()");
    momentumLinAngFree(rloc, v_artic, a_artic);
}

/* REVIEW[epic=Change,order=2]: Change articulationVel() for different systems*/
const Eigen::VectorXd&
rungeKutta4::articulationVel(double dimensional_time)
{
    m_velArtic = Eigen::VectorXd::Zero(3 * m_system->numParticles());

    double part0_x_vel = -m_U0 * cos(m_omega * dimensional_time);
    double part2_x_vel = m_U0 * cos(m_omega * dimensional_time + m_phase_shift);

    m_velArtic(0)     = part0_x_vel;
    m_velArtic(3 * 2) = part2_x_vel;

    return m_velArtic;
}

/* REVIEW[epic=Change,order=3]: Change articulationAcc() for different systems*/
const Eigen::VectorXd&
rungeKutta4::articulationAcc(double dimensional_time)
{
    m_accArtic = Eigen::VectorXd::Zero(3 * m_system->numParticles());

    double part0_x_acc = m_U0 * m_omega * sin(m_omega * dimensional_time);
    double part2_x_acc = -m_U0 * m_omega * sin(m_omega * dimensional_time + m_phase_shift);

    m_accArtic(0)     = part0_x_acc;
    m_accArtic(3 * 2) = part2_x_acc;

    return m_accArtic;
}

/* REVIEW[epic=Change,order=4]: Change rLoc() for different systems*/
const Eigen::Vector3d&
rungeKutta4::rLoc()
{
    m_RLoc = m_system->positions.segment<3>(3 * 1);
    return m_RLoc;
}

void
rungeKutta4::momentumLinAngFree(Eigen::Vector3d& r_loc, Eigen::VectorXd& v_artic,
                                Eigen::VectorXd a_artic)
{ // TODO
    /* ANCHOR: Solve for rigid body motion (rbm) tensors */
    // initialize variables
    Eigen::Matrix3d I       = Eigen::Matrix3d::Identity(3, 3);
    Eigen::Matrix3d rbmconn = Eigen::MatrixXd::Zero(6, 3 * m_system->numParticles());

    // assemble the rigid body motion connectivity tensor (Sigma)
    for (int i = 0; i < m_system->numParticles(); i++)
    {
        int i3{3 * i};

        Eigen::Vector3d n_dr = m_system->positions.segment<3>(i3);
        n_dr.noalias() -= r_loc;
        n_dr *= -1;
        Eigen::Matrix3d n_dr_cross;
        crossProdMat(n_dr, n_dr_cross);

        rbmconn.block<3, 3>(0, i3).noalias() = I;          // translation-translation couple
        rbmconn.block<3, 3>(3, i3).noalias() = n_dr_cross; // translation-rotation couple
    }

    // calculate M_tilde = Sigma * M_total * Sigma^T
    Eigen::MatrixXd M_tilde_hold = m_potHydro->mTotal() * rbmconn.transpose();
    Eigen::MatrixXd M_tilde      = rbmconn * M_tilde_hold;
    Eigen::MatrixXd M_tilde_inv  = M_tilde.inverse();

    /* ANCHOR: Solve for rigid body motion velocity components */
    // calculate P_script = Sigma * M_total * V_articulation
    Eigen::VectorXd P_script_hold = m_potHydro->mTotal() * v_artic;
    Eigen::VectorXd P_script      = rbmconn * P_script_hold;

    // calculate U_swim = - M_tilde_inv * P_script; U_swim has translation and rotation components
    Eigen::VectorXd U_swim = -M_tilde_inv * P_script;

    /* ANCHOR: Output velocity data back to m_system */
    // calculate U = Sigma^T * U_swim + v_artic
    m_system->velocities.noalias() = rbmconn.transpose() * U_swim;
    m_system->velocities.noalias() += v_artic;

    // update hydrodynamic force terms for acceleration components
    m_potHydro->updateForcesOnly();

    /* ANCHOR: Solve for rigid body motion acceleration components */
    // calculate b = a_artic + W * r + [v_artic ^ ]^T * Omega_C;  Omega_C = U_swim(3::5)
    Eigen::Vector3d Omega_C = U_swim.segment<3>(3);

    Eigen::Matrix3d W = Omega_C * Omega_C.transpose();
    W.noalias() -= Omega_C.squaredNorm() * I;

    Eigen::Matrix3d n_v_artic_cross;
    Eigen::VectorXd b = a_artic;

    for (int i = 0; i < m_system->numParticles(); i++)
    {
        int i3{3 * i};

        Eigen::Vector3d dr = m_system->positions.segment<3>(i3);
        dr.noalias() -= r_loc;
        b.segment<3>(i3).noalias() += W * dr;

        crossProdMat(-v_artic.segment<3>(i3), n_v_artic_cross);
        b.segment<3>(i3).noalias() += n_v_artic_cross * Omega_C;
    }

    // calculate gMUU = \nabla M_added : U U
    Eigen::VectorXd gMUU = m_potHydro->t2VelGrad();

    // calculate F_script = Sigma * (M_total * b + gMUU)
    Eigen::VectorXd F_script_hold = m_potHydro->mTotal() * b;
    F_script_hold.noalias() += gMUU;
    Eigen::VectorXd F_script = rbmconn * F_script_hold;

    // calculate A_swim = - M_tilde_inv * F_script; A_swim has translation and rotation components
    Eigen::VectorXd A_swim = -M_tilde_inv * F_script;

    /* ANCHOR: Output acceleration data back to m_system */
    // calculate A = Sigma^T * A_swim + b
    m_system->accelerations.noalias() = rbmconn.transpose() * A_swim;
    m_system->accelerations.noalias() += b;
}
