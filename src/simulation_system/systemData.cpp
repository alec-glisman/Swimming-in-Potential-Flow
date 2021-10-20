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

    spdlog::get(m_logName)->info("Setting general-use tensors");

    // Eigen3 parallelization
    int m_num_physical_cores = std::thread::hardware_concurrency();
    Eigen::setNbThreads(m_num_physical_cores); // for Eigen OpenMP parallelization

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
    spdlog::get(m_logName)->info("Setting D.o.F. parameters");
    if (m_image_system)
    {
        m_num_DoF         = 6 * (m_num_bodies / 2);
        m_num_constraints = (m_num_bodies / 2);
    }
    else
    {
        m_num_DoF         = 6 * m_num_bodies; // D.o.F. are linear and angular positions of body centers
        m_num_constraints = m_num_bodies;     // 1 unit quaternion constraint per body
    }

    spdlog::get(m_logName)->info("Degrees of freedom: {0}", m_num_DoF);
    spdlog::get(m_logName)->info("Number of constraints: {0}", m_num_constraints);

    // initialize other integrator parameters
    spdlog::get(m_logName)->info("Setting integrator parameters");
    m_t0 = m_t;
    spdlog::get(m_logName)->info("Initial integration time: {0}", m_t0);

    // initialize general-use tensors
    m_levi_cevita.setZero();

    m_levi_cevita(1, 2, 0) = 1;
    m_levi_cevita(2, 1, 0) = -1;

    m_levi_cevita(0, 2, 1) = -1;
    m_levi_cevita(2, 0, 1) = 1;

    m_levi_cevita(0, 1, 2) = 1;
    m_levi_cevita(1, 0, 2) = -1;

    m_kappa_tilde.setZero();

    m_kappa_tilde(0, 1, 3) = 1;
    m_kappa_tilde(1, 2, 3) = 1;
    m_kappa_tilde(2, 3, 3) = 1;

    m_kappa_tilde(0, 0, 4) = -1;
    m_kappa_tilde(1, 3, 4) = -1;
    m_kappa_tilde(2, 2, 4) = 1;

    m_kappa_tilde(0, 3, 5) = 1;
    m_kappa_tilde(1, 0, 5) = -1;
    m_kappa_tilde(2, 1, 5) = -1;

    m_kappa_tilde(0, 2, 6) = -1;
    m_kappa_tilde(1, 1, 6) = 1;
    m_kappa_tilde(2, 0, 6) = -1;

    // initialize particle vectors
    spdlog::get(m_logName)->info("Setting particle group (swimmer) ids");
    m_particle_group_id = Eigen::VectorXi::Zero(m_num_particles);

    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        const double body_num            = (m_particle_type_id.segment(0, particle_id + 1).array() == 1).count() - 1;
        m_particle_group_id(particle_id) = std::round(body_num); // convert data type

        spdlog::get(m_logName)->info("Particle {0} group id: {1}", particle_id + 1, m_particle_group_id(particle_id));
    }
    assert(m_particle_group_id(0) == 0 && "Particle 0 must belong to group 0");
    assert(m_particle_group_id(m_num_particles - 1) == m_num_bodies - 1 && "Particle N must belong to group N");

    // initialize kinematic vectors
    m_orientations_particles                     = Eigen::VectorXd::Zero(3 * m_num_particles);
    m_positions_particles_articulation_init_norm = Eigen::VectorXd::Zero(3 * m_num_particles);
    m_positions_particles_articulation           = Eigen::VectorXd::Zero(3 * m_num_particles);
    m_velocities_particles_articulation          = Eigen::VectorXd::Zero(3 * m_num_particles);
    m_accelerations_particles_articulation       = Eigen::VectorXd::Zero(3 * m_num_particles);

    // initialize constraint matrices
    m_Udwadia_A = Eigen::MatrixXd::Zero(m_num_constraints, m_num_DoF + m_num_constraints);
    m_Udwadia_b = Eigen::VectorXd::Zero(m_num_constraints);

    // initialize rigid body motion matrices
    m_rbm_conn             = Eigen::MatrixXd::Zero(6 * m_num_bodies, 3 * m_num_particles);
    m_psi_conv_quat_ang    = Eigen::MatrixXd::Zero(6 * m_num_bodies, 7 * m_num_bodies);
    m_rbm_conn_T_quat      = Eigen::MatrixXd::Zero(3 * m_num_particles, 7 * m_num_bodies);
    m_tens_rbm_conn_T_quat = Eigen::Tensor<double, 2>(3 * m_num_particles, 7 * m_num_bodies);
    m_tens_rbm_conn_T_quat.setZero();

    // initialize gradient matrices
    m_conv_body_2_part_dof      = Eigen::MatrixXd::Zero(7 * m_num_bodies, 3 * m_num_particles);
    m_tens_conv_body_2_part_dof = Eigen::Tensor<double, 2>(7 * m_num_bodies, 3 * m_num_particles);
    m_tens_conv_body_2_part_dof.setZero();
    m_rbm_conn_T_quat_grad = Eigen::Tensor<double, 3>(3 * m_num_particles, 7 * m_num_bodies, 7 * m_num_bodies);
    m_rbm_conn_T_quat_grad.setZero();

    // Set initial configuration orientation
    spdlog::get(m_logName)->info("Setting initial configuration orientation");
    m_positions_particles_articulation_init_norm.setZero();

    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        if (m_particle_type_id(particle_id) == 1)
        {
            continue; // continue to next loop as all elements are zero
        }

        const int particle_id_3{3 * particle_id};
        const int body_id_7{7 * m_particle_group_id(particle_id)};

        const Eigen::Vector3d disp =
            m_positions_particles.segment<3>(particle_id_3) - m_positions_bodies.segment<3>(body_id_7);

        m_positions_particles_articulation_init_norm.segment<3>(particle_id_3).noalias() = disp.normalized();

        spdlog::get(m_logName)->info("Particle {0} initial orientation: [{1:03.3f}, {2:03.3f}, {3:03.3f}]",
                                     particle_id + 1, m_positions_particles_articulation_init_norm(particle_id_3),
                                     m_positions_particles_articulation_init_norm(particle_id_3 + 1),
                                     m_positions_particles_articulation_init_norm(particle_id_3 + 2));
    }

    // initialize constraints
    spdlog::get(m_logName)->critical("Setting up temporary single-thread eigen device");
    Eigen::ThreadPool       thread_pool = Eigen::ThreadPool(1);
    Eigen::ThreadPoolDevice single_core_device(&thread_pool, 1);

    spdlog::get(m_logName)->info("Initializing constraints");
    update(single_core_device);

    // output data
    logData();

    spdlog::get(m_logName)->info("Initialization complete");
    spdlog::get(m_logName)->flush();
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

    assert(m_quaternions_particles.size() == 4 * m_num_particles &&
           "Particle orientation (unit quaternions) vector has incorrect length, not 4N.");
    assert(m_positions_particles.size() == 3 * m_num_particles &&
           "Particle position vector has incorrect length, not 3N.");
    assert(m_velocities_particles.size() == 3 * m_num_particles &&
           "Particle velocity vector has incorrect length, not 3N.");
    assert(m_accelerations_particles.size() == 3 * m_num_particles &&
           "Particle acceleration vector has incorrect length, not 3N.");

    assert(m_positions_bodies.size() == 7 * m_num_bodies && "Body position vector has incorrect length, not 7M.");
    assert(m_velocities_bodies.size() == 7 * m_num_bodies && "Body velocity vector has incorrect length, not 7M.");
    assert(m_accelerations_bodies.size() == 7 * m_num_bodies &&
           "Body acceleration vector has incorrect length, not 7M.");

    assert(m_particle_type_id.size() == m_num_particles && "Particle type ID vector has incorrect length, not N.");

    assert(m_num_particles > 0 && "Must have at least one particle to simulate.");
    assert(m_num_bodies > 0 && "Must have at least one body to simulate.");
    assert(m_num_spatial_dim == 3 && "Simulation framwork currently requires 3 spatial dimensions.");

    assert(m_tf > 0. && "Must have positive integration time.");
    assert(m_dt > 0.0 && "Integration time step must be positive.");
    assert(m_tf > m_t && "Final time must be greater than current time.");
    assert(m_tau > 0.0 && "Characteristic time must be greater than zero.");

    assert(m_fluid_density >= 0.0 && "Fluid density must be non-negative.");
    assert(m_particle_density >= 0.0 && "Particle density must be non-negative");

    assert(m_wca_epsilon >= 0.0 && "WCA_epsilon must be non-negative.");
    assert(m_wca_sigma >= 0.0 && "WCA_sigma must be non-negative.");
}

void
systemData::logData()
{
    /* ANCHOR: Output simulation data */
    spdlog::get(m_logName)->info("Starting logdata()");
    spdlog::get(m_logName)->info("time: {0}", m_t);

    /* ANCHOR: Output particle data */
    spdlog::get(m_logName)->info("Particle orientations:");
    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        const int particle_id_3{3 * particle_id};
        spdlog::get(m_logName)->info("\tParticle {0}: [{1:03.3f}, {2:03.3f}, {3:03.3f}]", particle_id + 1,
                                     m_orientations_particles(particle_id_3),
                                     m_orientations_particles(particle_id_3 + 1),
                                     m_orientations_particles(particle_id_3 + 2));
    }

    spdlog::get(m_logName)->info("Particle positions:");
    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        const int particle_id_3{3 * particle_id};
        spdlog::get(m_logName)->info("\tParticle {0}: [{1:03.3f}, {2:03.3f}, {3:03.3f}]", particle_id + 1,
                                     m_positions_particles(particle_id_3), m_positions_particles(particle_id_3 + 1),
                                     m_positions_particles(particle_id_3 + 2));
    }

    spdlog::get(m_logName)->info("Particle velocities:");
    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        const int particle_id_3{3 * particle_id};
        spdlog::get(m_logName)->info("\tParticle {0}: [{1:03.3f}, {2:03.3f}, {3:03.3f}]", particle_id + 1,
                                     m_velocities_particles(particle_id_3), m_velocities_particles(particle_id_3 + 1),
                                     m_velocities_particles(particle_id_3 + 2));
    }

    spdlog::get(m_logName)->info("Particle accelerations:");
    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        const int particle_id_3{3 * particle_id};
        spdlog::get(m_logName)->info("\tParticle {0}: [{1:03.3f}, {2:03.3f}, {3:03.3f}]", particle_id + 1,
                                     m_accelerations_particles(particle_id_3),
                                     m_accelerations_particles(particle_id_3 + 1),
                                     m_accelerations_particles(particle_id_3 + 2));
    }

    /* ANCHOR: Output particle *articulation* data */
    spdlog::get(m_logName)->info("Particle articulation positions:");
    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        const int particle_id_3{3 * particle_id};
        spdlog::get(m_logName)->info("\tParticle {0}: [{1:03.3f}, {2:03.3f}, {3:03.3f}]", particle_id + 1,
                                     m_positions_particles_articulation(particle_id_3),
                                     m_positions_particles_articulation(particle_id_3 + 1),
                                     m_positions_particles_articulation(particle_id_3 + 2));
    }

    spdlog::get(m_logName)->info("Particle articulation velocities:");
    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        const int particle_id_3{3 * particle_id};
        spdlog::get(m_logName)->info("\tParticle {0}: [{1:03.3f}, {2:03.3f}, {3:03.3f}]", particle_id + 1,
                                     m_velocities_particles_articulation(particle_id_3),
                                     m_velocities_particles_articulation(particle_id_3 + 1),
                                     m_velocities_particles_articulation(particle_id_3 + 2));
    }

    spdlog::get(m_logName)->info("Particle articulation accelerations:");
    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        const int particle_id_3{3 * particle_id};
        spdlog::get(m_logName)->info("\tParticle {0}: [{1:03.3f}, {2:03.3f}, {3:03.3f}]", particle_id + 1,
                                     m_accelerations_particles_articulation(particle_id_3),
                                     m_accelerations_particles_articulation(particle_id_3 + 1),
                                     m_accelerations_particles_articulation(particle_id_3 + 2));
    }

    /* ANCHOR: Output body data */
    spdlog::get(m_logName)->info("(linear) Body positions:");
    for (int body_id = 0; body_id < m_num_bodies; body_id++)
    {
        const int body_id_7{7 * body_id};
        spdlog::get(m_logName)->info("\tBody {0}: [{1:03.3f}, {2:03.3f}, {3:03.3f}]", body_id + 1,
                                     m_positions_bodies(body_id_7), m_positions_bodies(body_id_7 + 1),
                                     m_positions_bodies(body_id_7 + 2));
    }
    spdlog::get(m_logName)->info("(linear) Body velocities:");
    for (int body_id = 0; body_id < m_num_bodies; body_id++)
    {
        const int body_id_7{7 * body_id};
        spdlog::get(m_logName)->info("\tBody {0}: [{1:03.3f}, {2:03.3f}, {3:03.3f}]", body_id + 1,
                                     m_velocities_bodies(body_id_7), m_velocities_bodies(body_id_7 + 1),
                                     m_velocities_bodies(body_id_7 + 2));
    }
    spdlog::get(m_logName)->info("(linear) Body accelerations:");
    for (int body_id = 0; body_id < m_num_bodies; body_id++)
    {
        const int body_id_7{7 * body_id};
        spdlog::get(m_logName)->info("\tBody {0}: [{1:03.3f}, {2:03.3f}, {3:03.3f}]", body_id + 1,
                                     m_accelerations_bodies(body_id_7), m_accelerations_bodies(body_id_7 + 1),
                                     m_accelerations_bodies(body_id_7 + 2));
    }

    spdlog::get(m_logName)->info("Ending logdata()");
}

void
systemData::update(Eigen::ThreadPoolDevice& device)
{
    // NOTE: Internal particle orientation D.o.F. calculated 1st
    orientationsArticulation();

    // NOTE: Articulation functions calculated 2nd (need m_orientations_particles)
    positionsArticulation();
    velocitiesArticulation();
    accelerationsArticulation();

    // NOTE: Rigid body motion tensors calculated 3rd (need m_positions_particles_articulation)
    positionsParticlesfromBodies();
    rigidBodyMotionTensors(device);
    gradientChangeOfVariableTensors(device);

    // NOTE: Particle degrees of freedom calculated 4th (need rbm tensors)
    convertBody2ParticleDoF(device);

    // NOTE: Udwadia linear system calculated 5th
    udwadiaLinearSystem();
}

void
systemData::orientationsArticulation()
{
    m_orientations_particles.setZero();

    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        if (m_particle_type_id(particle_id) == 1)
        {
            continue; // continue to next loop as all elements are zero
        }

        const int particle_id_3{3 * particle_id};
        const int body_id_7{7 * m_particle_group_id(particle_id)};

        // Get unit quaternion for body that particle_id is a member of
        const Eigen::Quaterniond particle_quat(m_positions_bodies.segment<4>(body_id_7 + 3));
        assert(abs(particle_quat.norm() - 1.0) < 1e-12 && "Quaternion is not unitary");

        // Rotate particle from initial configuration using unit quaternion
        Eigen::Quaterniond orient_init;
        orient_init.w()   = 0.0;
        orient_init.vec() = m_positions_particles_articulation_init_norm.segment<3>(particle_id_3);

        const Eigen::Quaterniond orient_rot    = particle_quat * orient_init * particle_quat.inverse();
        const Eigen::Vector3d    orient_rot_3d = orient_rot.vec();
        assert(abs(orient_rot.w()) < 1e-12 && "0th component of particle orientation vector is not zero");

        // Output rotated orientation
        m_orientations_particles.segment<3>(particle_id_3).noalias() = orient_rot_3d;
    }
}

void
systemData::positionsArticulation()
{
    double t_dimensional{m_tau * m_t};

    // articulation acceleration magnitudes
    const double r1_mag = m_sys_spec_R_avg + (m_sys_spec_U0 / m_sys_spec_omega) * sin(m_sys_spec_omega * t_dimensional);
    const double r3_mag = m_sys_spec_R_avg - (m_sys_spec_U0 / m_sys_spec_omega) *
                                                 sin(m_sys_spec_omega * t_dimensional + m_sys_spec_phase_shift);

    Eigen::VectorXd particle_distances(m_num_particles);
    particle_distances << 0.0, r1_mag, r3_mag, 0.0, r1_mag, r3_mag;

    m_positions_particles_articulation.setZero();

    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        if (m_particle_type_id(particle_id) == 1)
        {
            continue; // continue to next loop as all elements are zero
        }

        const int particle_id_3{3 * particle_id};
        const int body_id_7{7 * m_particle_group_id(particle_id)};

        m_positions_particles_articulation.segment<3>(particle_id_3).noalias() =
            particle_distances(particle_id) * m_orientations_particles.segment<3>(particle_id_3);
    }
}

void
systemData::velocitiesArticulation()
{
    double t_dimensional{m_tau * m_t};

    // articulation velocity magnitudes
    const double v1_mag = m_sys_spec_U0 * cos(m_sys_spec_omega * t_dimensional);
    const double v3_mag = -m_sys_spec_U0 * cos(m_sys_spec_omega * t_dimensional + m_sys_spec_phase_shift);

    Eigen::VectorXd particle_velocities(m_num_particles);
    particle_velocities << 0.0, v1_mag, v3_mag, 0.0, v1_mag, v3_mag;

    // Zero and then calculate  m_velocities_particles_articulation
    m_velocities_particles_articulation.setZero();

    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        if (m_particle_type_id(particle_id) == 1)
        {
            continue; // continue to next loop as all elements are zero
        }

        const int particle_id_3{3 * particle_id};
        const int body_id_7{7 * m_particle_group_id(particle_id)};

        m_velocities_particles_articulation.segment<3>(particle_id_3).noalias() =
            particle_velocities(particle_id) * m_orientations_particles.segment<3>(particle_id_3);

        // Image system: flip x-y components, leave z unchanged
        if (particle_id >= (m_num_particles / 2))
        {
            m_velocities_particles_articulation.segment<2>(particle_id_3) *= -1;
        }
    }
}

void
systemData::accelerationsArticulation()
{
    double t_dimensional{m_tau * m_t};

    // articulation acceleration magnitudes
    const double a1_mag = -m_sys_spec_U0 * m_sys_spec_omega * sin(m_sys_spec_omega * t_dimensional);
    const double a3_mag =
        m_sys_spec_U0 * m_sys_spec_omega * sin(m_sys_spec_omega * t_dimensional + m_sys_spec_phase_shift);

    Eigen::VectorXd particle_accelerations(m_num_particles);
    particle_accelerations << 0.0, a1_mag, a3_mag, 0.0, a1_mag, a3_mag;

    // Zero and then calculate m_accelerations_particles_articulation
    m_accelerations_particles_articulation.setZero();

    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        if (m_particle_type_id(particle_id) == 1)
        {
            continue; // continue to next loop as all elements are zero
        }

        const int particle_id_3{3 * particle_id};
        const int body_id_7{7 * m_particle_group_id(particle_id)};

        m_accelerations_particles_articulation.segment<3>(particle_id_3).noalias() =
            particle_accelerations(particle_id) * m_orientations_particles.segment<3>(particle_id_3);

        // Image system: flip x-y components, leave z unchanged
        if (particle_id >= (m_num_particles / 2))
        {
            m_accelerations_particles_articulation.segment<2>(particle_id_3) *= -1;
        }
    }
}

void
systemData::positionsParticlesfromBodies()
{
    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        const int particle_id_3{3 * particle_id};
        const int body_id_7{7 * m_particle_group_id(particle_id)};

        m_positions_particles.segment<3>(particle_id_3).noalias() =
            m_positions_bodies.segment<3>(body_id_7) + m_positions_particles_articulation.segment<3>(particle_id_3);
    }
}

void
systemData::udwadiaLinearSystem()
{
    m_Udwadia_A.setZero();
    m_Udwadia_b.setZero();

    for (int body_id = 0; body_id < m_num_constraints; body_id++)
    {
        const int body_id_7{7 * body_id};

        m_Udwadia_A.block<1, 4>(body_id, body_id_7 + 3).noalias() = m_positions_bodies.segment<4>(body_id_7 + 3);
        m_Udwadia_b(body_id) = -m_velocities_bodies.segment<4>(body_id_7 + 3).squaredNorm();
    }
}

void
systemData::rigidBodyMotionTensors(Eigen::ThreadPoolDevice& device)
{
    /* ANCHOR: Compute m_rbm_conn */
    m_rbm_conn.setZero();

    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        if (m_particle_type_id(particle_id) == 1)
        {
            continue; // continue to next loop as all elements are zero
        }

        const int particle_id_3{3 * particle_id};
        const int body_id_6{6 * m_particle_group_id(particle_id)};

        // skew-symmetric matrix representation of (negative) cross product
        Eigen::Matrix3d n_dr_cross;
        crossProdMat(-m_positions_particles_articulation.segment<3>(particle_id_3), n_dr_cross);

        // rigid body motion connectivity tensor elements
        m_rbm_conn.block<3, 3>(body_id_6, particle_id_3).noalias()     = m_I3;       // translation-translation couple
        m_rbm_conn.block<3, 3>(body_id_6 + 3, particle_id_3).noalias() = n_dr_cross; // translation-rotation couple
    }

    /* ANCHOR: Compute m_psi_conv_quat_ang */
    m_psi_conv_quat_ang.setZero();

    for (int body_id = 0; body_id < m_num_bodies; body_id++)
    {
        const int body_id_6{6 * body_id};
        const int body_id_7{7 * body_id};

        // matrix E from quaterion of body k
        Eigen::Matrix<double, 3, 4> E_theta_k;
        eMatrix(m_positions_bodies.segment<4>(body_id_7 + 3), E_theta_k);

        // matrix elements of Psi
        m_psi_conv_quat_ang.block<3, 3>(body_id_6, body_id_7).noalias() = m_I3; // no conversion from linear components
        m_psi_conv_quat_ang.block<3, 4>(body_id_6 + 3, body_id_7 + 3).noalias() =
            2 * E_theta_k; // angular-quaternion velocity couple
    }

    /* ANCHOR: Compute m_rbm_conn_T_quat & m_tens_rbm_conn_T_quat */
    m_rbm_conn_T_quat.noalias() = m_rbm_conn.transpose() * m_psi_conv_quat_ang;
    m_tens_rbm_conn_T_quat      = TensorCast(m_rbm_conn_T_quat);
}

void
systemData::gradientChangeOfVariableTensors(Eigen::ThreadPoolDevice& device)
{
    /* ANCHOR: Compute m_conv_body_2_part_dof and m_tens_conv_body_2_part_dof */
    m_conv_body_2_part_dof.setZero();

    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        if (m_particle_type_id(particle_id) == 1)
        {
            continue; // continue to next loop as all elements are zero
        }

        const int particle_id_3{3 * particle_id};
        const int body_id_6{6 * m_particle_group_id(particle_id)};
        const int body_id_7{7 * m_particle_group_id(particle_id)};

        // (4 vector) displacement of particle from locater point moment arm
        Eigen::Vector4d dr         = Eigen::Vector4d::Zero(4, 1);
        dr.segment<3>(1).noalias() = m_positions_particles_articulation.segment<3>(particle_id_3);

        // quaternion product representation of particle displacement (transpose)
        Eigen::Matrix<double, 3, 4> P_j_tilde_T;
        eMatrix(dr, P_j_tilde_T);
        Eigen::Matrix<double, 4, 3> P_j_tilde = P_j_tilde_T.transpose();

        // change of variables gradient tensor elements
        m_conv_body_2_part_dof.block<3, 3>(body_id_7, particle_id_3).noalias() =
            -m_I3; // translation-translation couple
        m_conv_body_2_part_dof.block<3, 3>(body_id_7 + 4, particle_id_3).noalias() =
            m_psi_conv_quat_ang.block<3, 4>(body_id_6 + 3, body_id_7 + 3) *
            P_j_tilde; // quaternion-rotation couple (first row is zero)
    }

    m_tens_conv_body_2_part_dof = TensorCast(m_conv_body_2_part_dof);

    /* ANCHOR : Compute m_rbm_conn_T_quat_grad */
    m_rbm_conn_T_quat_grad.setZero();

    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        if (m_particle_type_id(particle_id) == 1)
        {
            // Continue to next loop as all elements are zero
            continue;
        }

        const int particle_id_3{3 * particle_id};
        const int body_id_6{6 * m_particle_group_id(particle_id)};
        const int body_id_7{7 * m_particle_group_id(particle_id)};

        // get moment arm to locater point from C
        const Eigen::Matrix3d two_r_cross_mat = -2 * m_rbm_conn.block<3, 3>(body_id_6 + 3, particle_id_3);

        // get E matrix representation of body quaternion from m_psi_conv_quat_ang
        const Eigen::Matrix<double, 3, 4> two_E_body = m_psi_conv_quat_ang.block<3, 4>(body_id_6 + 3, body_id_7 + 3);

        // get change of variable matrix element D_{\alpha}
        const Eigen::Matrix<double, 7, 3> n_D_alpha = -m_conv_body_2_part_dof.block<7, 3>(body_id_7, particle_id_3);

        /* ANCHOR: tensor contractions to produce result */

        // convert Eigen::Matrix --> Eigen::Tensor
        const Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3>> tens_two_r_cross_mat = TensorCast(two_r_cross_mat);
        const Eigen::TensorFixedSize<double, Eigen::Sizes<3, 4>> tens_two_E_body      = TensorCast(two_E_body);
        const Eigen::TensorFixedSize<double, Eigen::Sizes<7, 3>> tens_n_D_alpha       = TensorCast(n_D_alpha);

        // 1) Compute result_{i m k} = levi_cevita{l i m} n_D_alpha_{k l}
        const Eigen::array<Eigen::IndexPair<int>, 1> contract_lim_kl = {Eigen::IndexPair<int>(0, 1)}; // {i m k}
        Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3, 7>> d_rCrossMat_d_xi;
        d_rCrossMat_d_xi.device(device) = m_levi_cevita.contract(tens_n_D_alpha, contract_lim_kl);

        // 2) Compute result_{i j k} = d_rCrossMat_d_xi_{i m k} two_E_body{m j}
        const Eigen::array<Eigen::IndexPair<int>, 1> contract_imk_mj = {
            Eigen::IndexPair<int>(1, 0)}; // {i k j}, must shuffle
        const Eigen::array<int, 3> swap_last_two_indices({0, 2, 1});

        Eigen::TensorFixedSize<double, Eigen::Sizes<3, 7, 4>> d_rCrossMat_d_xi_times_two_E_body_preshuffle;
        d_rCrossMat_d_xi_times_two_E_body_preshuffle.device(device) =
            d_rCrossMat_d_xi.contract(tens_two_E_body, contract_imk_mj);

        Eigen::TensorFixedSize<double, Eigen::Sizes<3, 4, 7>> d_rCrossMat_d_xi_times_two_E_body;
        d_rCrossMat_d_xi_times_two_E_body.device(device) =
            d_rCrossMat_d_xi_times_two_E_body_preshuffle.shuffle(swap_last_two_indices);

        // 3) Compute result_{i j k} = two_r_cross_mat_{i m} Kappa_tilde_{m j k} (same indices as (2))
        const Eigen::array<Eigen::IndexPair<int>, 1> contract_im_mjk = {Eigen::IndexPair<int>(1, 0)}; // {i j k}
        Eigen::TensorFixedSize<double, Eigen::Sizes<3, 4, 7>> rCrossMat_times_d_E_body_d_xi;
        rCrossMat_times_d_E_body_d_xi.device(device) = tens_two_r_cross_mat.contract(m_kappa_tilde, contract_im_mjk);

        // Compute and output matrix element
        Eigen::TensorFixedSize<double, Eigen::Sizes<3, 4, 7>> angular_gradient;
        angular_gradient.device(device) = -d_rCrossMat_d_xi_times_two_E_body - rCrossMat_times_d_E_body_d_xi;

        // output matrix indices
        const Eigen::array<Eigen::Index, 3> offsets    = {particle_id_3, body_id_7 + 3, body_id_7};
        const Eigen::array<Eigen::Index, 3> extents    = {3, 4, 7};
        m_rbm_conn_T_quat_grad.slice(offsets, extents) = angular_gradient;
    }
}

void
systemData::convertBody2ParticleDoF(Eigen::ThreadPoolDevice& device)
{
    /* ANCHOR: Convert velocity D.o.F. */
    m_velocities_particles.noalias() = m_rbm_conn_T_quat * m_velocities_bodies;
    m_velocities_particles.noalias() += m_velocities_particles_articulation;

    /* ANCHOR: Convert acceleration D.o.F. */

    // Eigen::Tensor contraction indices
    const Eigen::array<Eigen::IndexPair<int>, 0> outer_product   = {}; // no contractions
    const Eigen::array<Eigen::IndexPair<int>, 2> contract_ijk_jk = {
        Eigen::IndexPair<int>(1, 0), Eigen::IndexPair<int>(2, 1)}; // {i j k}, {j k} --> {i}

    const Eigen::Tensor<double, 1> xi    = TensorCast(m_velocities_bodies);
    const Eigen::Tensor<double, 2> xi_xi = xi.contract(xi, outer_product);

    Eigen::Tensor<double, 1> velocities_rbm_particles = Eigen::Tensor<double, 1>(3 * m_num_particles);
    velocities_rbm_particles.device(device)           = m_rbm_conn_T_quat_grad.contract(xi_xi, contract_ijk_jk);

    m_accelerations_bodies.noalias() = MatrixCast(velocities_rbm_particles, 3 * m_num_particles, 1, device);
    m_accelerations_bodies.noalias() += m_accelerations_particles_articulation;
}