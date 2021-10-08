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
    kappa_tilde(1, 3, 4) = -1;
    kappa_tilde(2, 2, 4) = 1;

    kappa_tilde(0, 3, 5) = 1;
    kappa_tilde(1, 0, 5) = -1;
    kappa_tilde(2, 1, 5) = -1;

    kappa_tilde(0, 2, 6) = -1;
    kappa_tilde(1, 1, 6) = 1;
    kappa_tilde(2, 0, 6) = -1;

    // initialize constraints
    spdlog::get(m_logName)->info("Initializing constraints");
    updateConstraints(0.0);
}

void
systemData::updateConstraints(double time)
{
    locaterPointLocations();

    velocitiesArticulation(time);
    accelerationsArticulation(time);

    rigidBodyMotionTensors();
    gradientChangeOfVariableTensors();

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

void
systemData::rigidBodyMotionTensors()
{
    /* ANCHOR: Compute m_rbm_conn */
    m_rbm_conn = Eigen::MatrixXd::Zero(6 * m_num_bodies, 3 * m_num_particles);

    for (int i = 0; i < m_num_particles; i++)
    {
        int i3{3 * i};

        const Eigen::Vector3d R_i = m_positions_particles.segment<3>(i3);

        // Eigen3 is column-major by default so loop over rows (bodies) in inner loop
        for (int j = 0; j < m_num_bodies; j++)
        {
            int j6{6 * j};
            int j7{7 * j};

            // (negative) moment arm of particle about its body's locater position
            Eigen::Vector3d n_dr = m_positions_bodies.segment<3>(j7);
            n_dr.noalias() -= R_i;

            // skew-symmetric matrix representation of cross product
            Eigen::Matrix3d n_dr_cross;
            crossProdMat(n_dr, n_dr_cross);

            // rigid body motion connectivity tensor elements
            m_rbm_conn.block<3, 3>(j6, i3).noalias() = m_I; // translation-translation couple
            m_rbm_conn.block<3, 3>(j6 + 3, i3).noalias() =
                n_dr_cross; // translation-rotation couple
        }
    }

    /* ANCHOR: Compute m_psi_conv_quat_ang */
    m_psi_conv_quat_ang = Eigen::MatrixXd::Zero(6 * m_num_bodies, 7 * m_num_bodies);

    for (int k = 0; k < m_num_bodies; k++)
    {
        int k6{6 * k};
        int k7{7 * k};

        // matrix E from quaterion of body k
        Eigen::Matrix<double, 3, 4> E_theta_k;
        eMatrix(m_positions_bodies.segment<4>(k7 + 3), E_theta_k);

        // matrix elements of Psi
        m_psi_conv_quat_ang.block<3, 3>(k6, k7).noalias() =
            m_I; // no conversion from linear components
        m_psi_conv_quat_ang.block<3, 4>(k6 + 3, k7 + 3).noalias() =
            2 * E_theta_k; // angular-quaternion velocity couple
    }

    /* ANCHOR: Compute m_C_conv_quat_part */
    m_C_conv_quat_part.noalias() = m_rbm_conn.transpose() * m_psi_conv_quat_ang;
}

void
systemData::gradientChangeOfVariableTensors()
{
    /* ANCHOR: Compute m_D_conv_quat_part */
    m_D_conv_quat_part = Eigen::MatrixXd::Zero(7 * m_num_bodies, 3 * m_num_particles);

    for (int i = 0; i < m_num_particles; i++)
    {
        int i3{3 * i};

        // 4-vector version of particle i position (prepend zero element)
        Eigen::Vector4d R_i         = Eigen::Vector4d::Zero(4, 1);
        R_i.segment<3>(1).noalias() = m_positions_particles.segment<3>(i3);

        // Eigen3 is column-major by default so loop over rows (bodies) in inner loop
        for (int j = 0; j < m_num_bodies; j++)
        {
            int j7{7 * j};

            // moment arm of particle i about its body j locater position
            Eigen::Vector4d dr = R_i;
            dr.segment<3>(1).noalias() -= m_positions_bodies.segment<3>(j7);

            // matrix E from quaterion of body j
            Eigen::Matrix<double, 3, 4> E_theta_j;
            eMatrix(m_positions_bodies.segment<4>(j7 + 3), E_theta_j);

            // matrix P_tilde from moment arm
            Eigen::Matrix<double, 3, 4> P_i_tilde_hold;
            eMatrix(dr, P_i_tilde_hold);

            // change of variables gradient tensor elements
            m_D_conv_quat_part.block<3, 3>(j7, i3).noalias() =
                -m_I; // translation-translation couple
            m_D_conv_quat_part.block<3, 3>(j7 + 4, i3).noalias() =
                2 * E_theta_j *
                P_i_tilde_hold.transpose(); // quaternion-rotation couple (first row is zero)
        }
    }
}