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
    m_system->velocities    = articulationVel(0.0);
    m_system->accelerations = articulationAcc(0.0);
    // calculate locater point motion via linear and angular momentum free conditions
    int             first_idx{0};
    int             last_idx{3 * m_system->numParticles() - 1};
    Eigen::Vector3d rloc = rLoc();
    momentumLinAngFree(rloc, first_idx, last_idx);
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
rungeKutta4::momentumLinAngFree(Eigen::Vector3d& r_loc, int first_idx, int last_idx)
{
    // TODO

    /* ANCHOR: Generate linear system for locater point motion */

    /* ANCHOR: Solve linear system for locater point motion */

    /* ANCHOR: Output kinematic data back to m_system */
}
