//
// Created by Alec Glisman on 07/30/21
//

#include <systemData.hpp>

systemData::systemData(std::string inputGSDFile, std::string outputDir)
    : m_inputGSDFile(std::move(inputGSDFile)), m_outputDir(std::move(outputDir))
{
    // Initialize logger
    m_logFile   = m_outputDir + "/logs/" + m_logName + "-log.txt";
    auto logger = spdlog::basic_logger_mt(m_logName, m_logFile);
    spdlog::get(m_logName)->info("Initializing system data");
    spdlog::get(m_logName)->info("Output path: {0}", m_outputDir);

    // set "identity" tensors
    spdlog::get(m_logName)->info("Setting general-use tensors");
    m_I_tilde       = m_I;
    m_I_tilde(2, 2) = -1.0;

    spdlog::get(m_logName)->info("Constructor complete");
    spdlog::get(m_logName)->flush();
}

systemData::~systemData()
{
    spdlog::get(m_logName)->info("systemData destructor called");
    gsd_close(m_handle.get());
    spdlog::get(m_logName)->flush();
}

void
systemData::initializeData()
{
    spdlog::get(m_logName)->info("Running initializeData()");

    // load and check data from GSD
    parseGSD();

    // initialize D.o.F. parameters
    m_num_DoF         = 6 * m_num_bodies; // D.o.F. are linear and angular positions of body centers
    m_num_constraints = m_num_bodies;     // 1 unit quaternion constraint per body

    // initialize general-use tensors
    levi_cevita.setZero();
    levi_cevita(1, 2, 0) = 1;
    levi_cevita(2, 1, 0) = -1;

    levi_cevita(0, 2, 1) = -1;
    levi_cevita(2, 0, 1) = 1;

    levi_cevita(0, 1, 2) = 1;
    levi_cevita(1, 0, 2) = -1;

    kappa_tilde.setZero();
    kappa_tilde(0, 1, 3) = 1;
    kappa_tilde(1, 2, 3) = 1;
    kappa_tilde(2, 3, 3) = 1;

    kappa_tilde(0, 0, 4) = -1;

    // initialize constraints
    spdlog::get(m_logName)->info("Initializing constraints");
    updateConstraints(0.0);
}

void
systemData::updateConstraints(double time)
{
    velocitiesArticulation(time);
    accelerationsArticulation(time);
    locaterPointLocations();

    udwadiaLinearSystem(time);
}

void
systemData::parseGSD()
{
    spdlog::get(m_logName)->info("Running parseGSD()");

    // parse GSD file and load data into *this
    m_gsdUtil    = std::make_shared<GSDUtil>(shared_from_this());
    m_GSD_parsed = true;

    // verify data is not a-physical
    checkInput();
}

void
systemData::checkInput()
{
    spdlog::get(m_logName)->info("Input checking assertions");

    assert(m_positions_particles.size() == 4 * m_num_particles &&
           "Particle orientation (unit quaternions) vector has incorrect length, not 4N.");
    assert(m_positions_particles.size() == 3 * m_num_particles &&
           "Particle position vector has incorrect length, not 3N.");
    assert(m_velocities_particles.size() == 3 * m_num_particles &&
           "Particle velocity vector has incorrect length, not 3N.");
    assert(m_accelerations_particles.size() == 3 * m_num_particles &&
           "Particle acceleration vector has incorrect length, not 3N.");

    assert(m_positions_bodies.size() == 7 * m_num_particles &&
           "Body position vector has incorrect length, not 7N.");
    assert(m_velocities_bodies.size() == 7 * m_num_particles &&
           "Body velocity vector has incorrect length, not 7N.");
    assert(m_accelerations_bodies.size() == 7 * m_num_particles &&
           "Body acceleration vector has incorrect length, not 7N.");

    assert(m_particle_type_id.size() == m_num_particles &&
           "Particle type ID vector has incorrect length, not N.");

    assert(m_num_particles > 0 && "Must have at least one particle to simulate.");
    assert(m_num_bodies > 0 && "Must have at least one body to simulate.");
    assert(m_num_spatial_dim == 3 &&
           "Simulation framwork currently requires 3 spatial dimensions.");

    assert(m_tf > 0. && "Must have positive integration time.");
    assert(m_dt > 0.0 && "Integration time step must be positive.");
    assert(m_tf > m_t && "Final time must be greater than current time.");
    assert(m_tau > 0.0 && "Characteristic time must be greater than zero.");

    assert(m_fluid_density >= 0.0 && "Fluid density must be non-negative.");
    assert(m_particle_density >= 0.0 && "Particle density must be non-negative");

    assert(m_wca_epsilon >= 0.0 && "WCA_epsilon must be non-negative.");
    assert(m_wca_sigma >= 0.0 && "WCA_sigma must be non-negative.");
}

/* REVIEW[epic=Change,order=3]: Change assignment of m_velocities_particles_articulation for
 * different systems */
void
systemData::velocitiesArticulation(double time)
{
    double dimensional_time{m_tau * time};

    // ANCHOR: Orientation vectors, q = R_1 - R_3
    Eigen::Vector3d q =
        m_positions_particles.segment<3>(3 * 0) - m_positions_particles.segment<3>(3 * 2);
    q.normalize();
    Eigen::Vector3d q_tilde = m_I_tilde * q;

    // articulation velocity magnitudes
    const double v1_mag = m_sys_spec_U0 * cos(m_sys_spec_omega * dimensional_time);
    const double v3_mag =
        m_sys_spec_U0 * cos(m_sys_spec_omega * dimensional_time + m_sys_spec_phase_shift);

    // Zero and then calculate  m_velocities_particles_articulation
    m_velocities_particles_articulation = Eigen::VectorXd::Zero(3 * m_num_particles);
    m_velocities_particles_articulation.segment<3>(3 * 0).noalias() = v1_mag * q;
    m_velocities_particles_articulation.segment<3>(3 * 2).noalias() = v3_mag * q;
    m_velocities_particles_articulation.segment<3>(3 * 3).noalias() = v1_mag * q_tilde;
    m_velocities_particles_articulation.segment<3>(3 * 5).noalias() = v3_mag * q_tilde;
}

/* REVIEW[epic=Change,order=4]: Change assignment of m_accelerations_particles_articulation for
 * different systems */
void
systemData::accelerationsArticulation(double time)
{
    double dimensional_time{m_tau * time};

    // ANCHOR: Orientation vectors, q = R_1 - R_3
    Eigen::Vector3d q =
        m_positions_particles.segment<3>(3 * 0) - m_positions_particles.segment<3>(3 * 2);
    q.normalize();
    Eigen::Vector3d q_tilde = m_I_tilde * q;

    // articulation acceleration magnitudes
    const double a1_mag =
        -m_sys_spec_U0 * m_sys_spec_omega * sin(m_sys_spec_omega * dimensional_time);
    const double a3_mag = -m_sys_spec_U0 * m_sys_spec_omega *
                          sin(m_sys_spec_omega * dimensional_time + m_sys_spec_phase_shift);

    // Zero and then calculate m_accelerations_particles_articulation
    m_accelerations_particles_articulation = Eigen::VectorXd::Zero(3 * m_num_particles);
    m_accelerations_particles_articulation.segment<3>(3 * 0).noalias() = a1_mag * q;
    m_accelerations_particles_articulation.segment<3>(3 * 2).noalias() = a3_mag * q;
    m_accelerations_particles_articulation.segment<3>(3 * 3).noalias() = a1_mag * q_tilde;
    m_accelerations_particles_articulation.segment<3>(3 * 5).noalias() = a3_mag * q_tilde;
}

/* REVIEW[epic=Change,order=5]: Change assignment of m_positions_locater_particles for different
 * systems */
void
systemData::locaterPointLocations()
{
    m_positions_locater_particles = Eigen::VectorXd::Zero(3 * m_num_bodies);
    int body_count{0};

    for (int i = 0; i < m_num_particles; i++)
    {
        // only fill for locater particles
        if (m_particle_type_id(i) == 1)
        {
            m_positions_bodies.segment<3>(3 * body_count).noalias() =
                m_positions_particles.segment<3>(3 * i);

            body_count += 1;
        }
    }
    assert(body_count == m_num_bodies && "Incorrect number of bodies filled");
}

/* REVIEW[epic=Change,order=1]: Change udwadiaLinearSystem() for different systems
 */
void
systemData::udwadiaLinearSystem(double time)
{
    // Initialize matrices
    m_Udwadia_A = Eigen::MatrixXd::Zero(m_num_constraints, m_num_DoF);
    m_Udwadia_b = Eigen::VectorXd::Zero(m_num_constraints);

    // TODO: Implement the constraint linear system for the unit quaternions
}