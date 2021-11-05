//
// Created by Alec Glisman on 07/30/21
//

#include <systemData.hpp>

systemData::systemData(std::string inputGSDFile, std::string outputDir)
    : m_inputGSDFile(inputGSDFile), m_outputDir(outputDir)
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
    spdlog::drop(m_logName);
}

void
systemData::initializeData()
{
    spdlog::get(m_logName)->info("Running initializeData()");

    // load and check data from GSD
    parseGSD();

    // initialize D.o.F. parameters
    spdlog::get(m_logName)->info("Setting D.o.F. parameters");
    m_num_DoF         = 6 * m_num_bodies; // D.o.F. are linear and angular positions of body centers
    m_num_constraints = m_num_bodies;     // 1 unit quaternion constraint per body

    if (m_image_system)
    {
        m_num_DoF /= 2;
        m_num_constraints /= 2;
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

    m_kappa.setZero();

    m_kappa(0, 0, 3) = 1;
    m_kappa(1, 1, 3) = 1;
    m_kappa(2, 2, 3) = 1;
    m_kappa(3, 3, 3) = 1;

    m_kappa(0, 1, 4) = -1;
    m_kappa(1, 0, 4) = 1;
    m_kappa(2, 3, 4) = -1;
    m_kappa(3, 2, 4) = 1;

    m_kappa(0, 2, 5) = -1;
    m_kappa(1, 3, 5) = 1;
    m_kappa(2, 0, 5) = 1;
    m_kappa(3, 1, 5) = -1;

    m_kappa(0, 3, 6) = -1;
    m_kappa(1, 2, 6) = -1;
    m_kappa(2, 1, 6) = 1;
    m_kappa(3, 0, 6) = 1;

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

    // tensor lengths
    const int n3{3 * m_num_particles};
    const int n4{4 * m_num_particles};
    const int n7{7 * m_num_particles};

    const int m7{7 * m_num_bodies};

    // initialize kinematic vectors
    m_orientations_particles                     = Eigen::VectorXd::Zero(n3);
    m_positions_particles_articulation_init_norm = Eigen::VectorXd::Zero(n3);
    m_positions_particles_articulation           = Eigen::VectorXd::Zero(n3);
    m_velocities_particles_articulation          = Eigen::VectorXd::Zero(n7);
    m_accelerations_particles_articulation       = Eigen::VectorXd::Zero(n7);

    // initialize constraint matrices
    m_Udwadia_A = Eigen::MatrixXd::Zero(m_num_constraints, m_num_DoF + m_num_constraints);
    m_Udwadia_b = Eigen::VectorXd::Zero(m_num_constraints);

    // initialize rigid body motion matrices
    m_rbm_conn          = Eigen::MatrixXd::Zero(m7, n7);
    m_psi_conv_quat_ang = Eigen::MatrixXd::Zero(m7, m7);
    m_zeta              = Eigen::MatrixXd::Zero(m7, n7);
    m_tens_zeta         = Eigen::Tensor<double, 2>(m7, n7);
    m_tens_zeta.setZero();

    // initialize gradient matrices
    m_chi      = Eigen::MatrixXd::Zero(m7, n7);
    m_tens_chi = Eigen::Tensor<double, 2>(m7, n7);
    m_tens_chi.setZero();
    m_tens_grad_zeta = Eigen::Tensor<double, 3>(m7, n7, m7);
    m_tens_grad_zeta.setZero();

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

        /// @todo: also load data from frame0 of GSD for initial positions so that simulations can continue from
        /// previous GSD
        /// @todo: add a check that the first frame's unit quaternions are all identity quaternions (1, 0, 0, 0)^T
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
    spdlog::get(m_logName)->info("Starting parseGSD()");

    // parse GSD file and load data into *this
    m_gsdUtil    = std::make_shared<GSDUtil>(shared_from_this());
    m_GSD_parsed = true;

    // verify data is not a-physical
    checkInput();

    spdlog::get(m_logName)->info("Ending parseGSD()");
    spdlog::get(m_logName)->flush();
}

void
systemData::checkInput()
{
    spdlog::get(m_logName)->info("Input checking assertions");

    assert(m_quaternions_particles.size() == 4 * m_num_particles &&
           "Particle orientation (unit quaternions) vector has incorrect length, not 4N.");
    assert(m_positions_particles.size() == 3 * m_num_particles &&
           "Particle position vector has incorrect length, not 3N.");
    assert(m_velocities_particles.size() == 7 * m_num_particles &&
           "Particle velocity vector has incorrect length, not 7N.");
    assert(m_accelerations_particles.size() == 7 * m_num_particles &&
           "Particle acceleration vector has incorrect length, not 7N.");

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
        const int particle_id_7{7 * particle_id};
        spdlog::get(m_logName)->info(
            "\tParticle {0}: [{1:03.3f}, {2:03.3f}, {3:03.3f}, {4:03.3f}, {5:03.3f}, {6:03.3f}, {7:03.3f}]",
            particle_id + 1, m_velocities_particles_articulation(particle_id_7),
            m_velocities_particles_articulation(particle_id_7 + 1),
            m_velocities_particles_articulation(particle_id_7 + 2),
            m_velocities_particles_articulation(particle_id_7 + 3),
            m_velocities_particles_articulation(particle_id_7 + 4),
            m_velocities_particles_articulation(particle_id_7 + 5),
            m_velocities_particles_articulation(particle_id_7 + 6));
    }

    spdlog::get(m_logName)->info("Particle articulation accelerations:");
    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        const int particle_id_7{7 * particle_id};
        spdlog::get(m_logName)->info(
            "\tParticle {0}: [{1:03.3f}, {2:03.3f}, {3:03.3f}, {4:03.3f}, {5:03.3f}, {6:03.3f}, {7:03.3f}]",
            particle_id + 1, m_accelerations_particles_articulation(particle_id_7),
            m_accelerations_particles_articulation(particle_id_7 + 1),
            m_accelerations_particles_articulation(particle_id_7 + 2),
            m_accelerations_particles_articulation(particle_id_7 + 3),
            m_accelerations_particles_articulation(particle_id_7 + 4),
            m_accelerations_particles_articulation(particle_id_7 + 5),
            m_accelerations_particles_articulation(particle_id_7 + 6));
    }

    /* ANCHOR: Output body data */
    spdlog::get(m_logName)->info("Body positions:");
    for (int body_id = 0; body_id < m_num_bodies; body_id++)
    {
        const int body_id_7{7 * body_id};
        spdlog::get(m_logName)->info(
            "\tBody {0}: [{1:03.3f}, {2:03.3f}, {3:03.3f}, {4:03.3f}, {5:03.3f}, {6:03.3f}, {7:03.3f}]", body_id + 1,
            m_positions_bodies(body_id_7), m_positions_bodies(body_id_7 + 1), m_positions_bodies(body_id_7 + 2),
            m_positions_bodies(body_id_7 + 3), m_positions_bodies(body_id_7 + 4), m_positions_bodies(body_id_7 + 5),
            m_positions_bodies(body_id_7 + 6));
    }
    spdlog::get(m_logName)->info("Body velocities:");
    for (int body_id = 0; body_id < m_num_bodies; body_id++)
    {
        const int body_id_7{7 * body_id};
        spdlog::get(m_logName)->info(
            "\tBody {0}: [{1:03.3f}, {2:03.3f}, {3:03.3f}, {4:03.3f}, {5:03.3f}, {6:03.3f}, {7:03.3f}]", body_id + 1,
            m_velocities_bodies(body_id_7), m_velocities_bodies(body_id_7 + 1), m_velocities_bodies(body_id_7 + 2),
            m_velocities_bodies(body_id_7 + 3), m_velocities_bodies(body_id_7 + 4), m_velocities_bodies(body_id_7 + 5),
            m_velocities_bodies(body_id_7 + 6));
    }
    spdlog::get(m_logName)->info("Body accelerations:");
    for (int body_id = 0; body_id < m_num_bodies; body_id++)
    {
        const int body_id_7{7 * body_id};
        spdlog::get(m_logName)->info(
            "\tBody {0}: [{1:03.3f}, {2:03.3f}, {3:03.3f}, {4:03.3f}, {5:03.3f}, {6:03.3f}, {7:03.3f}]", body_id + 1,
            m_accelerations_bodies(body_id_7), m_accelerations_bodies(body_id_7 + 1),
            m_accelerations_bodies(body_id_7 + 2), m_accelerations_bodies(body_id_7 + 3),
            m_accelerations_bodies(body_id_7 + 4), m_accelerations_bodies(body_id_7 + 5),
            m_accelerations_bodies(body_id_7 + 6));
    }

    spdlog::get(m_logName)->info("Ending logdata()");
    spdlog::get(m_logName)->flush();
}

void
systemData::update(Eigen::ThreadPoolDevice& device)
{
    // NOTE: Internal particle orientation D.o.F. calculated 1st
    convertBody2ParticleOrient();

    // NOTE: Articulation functions calculated 2nd (need m_orientations_particles)
    positionsArticulation();
    velocitiesArticulation();
    accelerationsArticulation();

    // NOTE: Rigid body motion tensors calculated 3rd (need m_positions_particles_articulation)
    convertBody2ParticlePos();
    rigidBodyMotionTensors(device);
    gradientChangeOfVariableTensors(device);

    // NOTE: Particle degrees of freedom calculated 4th (need rbm tensors)
    convertBody2ParticleVelAcc(device);

    // NOTE: Udwadia linear system calculated 5th
    udwadiaLinearSystem();
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
        const int particle_id_7{7 * particle_id};
        const int body_id_7{7 * m_particle_group_id(particle_id)};

        m_velocities_particles_articulation.segment<3>(particle_id_7).noalias() =
            particle_velocities(particle_id) * m_orientations_particles.segment<3>(particle_id_3);

        // Image system: flip x-y components, leave z unchanged
        // TODO: Related to Issue #3
        if (particle_id >= (m_num_particles / 2))
        {
            m_velocities_particles_articulation.segment<2>(particle_id_7) *= -1;
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
        const int particle_id_7{7 * particle_id};
        const int body_id_7{7 * m_particle_group_id(particle_id)};

        m_accelerations_particles_articulation.segment<3>(particle_id_7).noalias() =
            particle_accelerations(particle_id) * m_orientations_particles.segment<3>(particle_id_3);

        // Image system: flip x-y components, leave z unchanged
        // TODO: Related to Issue #3
        if (particle_id >= (m_num_particles / 2))
        {
            m_accelerations_particles_articulation.segment<2>(particle_id_7) *= -1;
        }
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

        const int particle_id_3{3 * particle_id};
        const int particle_id_7{7 * particle_id};
        const int body_id_7{7 * m_particle_group_id(particle_id)};

        // skew-symmetric matrix representation of cross product
        Eigen::Matrix3d mat_dr_cross;
        crossProdMat(m_positions_particles_articulation.segment<3>(particle_id_3), mat_dr_cross);

        // rigid body motion connectivity tensor elements
        m_rbm_conn.block<3, 3>(body_id_7, particle_id_7).noalias() = m_I3; // translation-translation couple

        /// @FIXME: The teaching SD to swim paper had a negative sign here, but my derivations do not have that
        m_rbm_conn.block<3, 3>(body_id_7 + 3, particle_id_7).noalias() = mat_dr_cross; // translation-rotation couple

        m_rbm_conn.block<3, 3>(body_id_7 + 3, particle_id_7 + 3).noalias() = m_I3; // rotation-rotation couple
    }

    /* ANCHOR: Compute m_psi_conv_quat_ang */
    m_psi_conv_quat_ang.setZero();

    for (int body_id = 0; body_id < m_num_bodies; body_id++)
    {
        const int body_id_7{7 * body_id};

        // matrix E from quaterion of body k
        Eigen::Matrix<double, 3, 4> E_theta_k;
        eMatrix(m_positions_bodies.segment<4>(body_id_7 + 3), E_theta_k);

        // matrix elements of Psi
        m_psi_conv_quat_ang.block<3, 3>(body_id_7, body_id_7).noalias() = m_I3; // no conversion from linear components
        m_psi_conv_quat_ang.block<3, 4>(body_id_7 + 3, body_id_7 + 3).noalias() =
            2 * E_theta_k; // angular-quaternion velocity couple
    }

    /* ANCHOR: Compute m_zeta & m_tens_zeta */
    m_zeta.noalias() = m_psi_conv_quat_ang.transpose() * m_rbm_conn;
    m_tens_zeta      = TensorCast(m_zeta);
}

void
systemData::gradientChangeOfVariableTensors(Eigen::ThreadPoolDevice& device)
{
    /* ANCHOR: Compute m_chi and m_tens_chi */
    m_chi.setZero();

    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        const int particle_id_3{3 * particle_id};
        const int particle_id_7{7 * particle_id};
        const int body_id_7{7 * m_particle_group_id(particle_id)};

        const bool is_locater{m_particle_group_id(particle_id) == 1};

        /// S_alpha = {-1 for non-locater particles, +1 for locater particles}
        const int s_part{-1 + 2 * is_locater};

        // (4 vector) displacement of particle from locater point moment arm
        Eigen::Vector4d dr_init         = Eigen::Vector4d::Zero(4, 1);
        dr_init.segment<3>(1).noalias() = m_positions_particles_articulation_init_norm.segment<3>(particle_id_3);

        const double r_alpha_norm{m_positions_particles_articulation.segment<3>(particle_id_3).norm()};
        dr_init *= r_alpha_norm; // incorporate norm factor here rather than in chi calculation for numerical efficiency

        // quaternion product representation of particle displacement (transpose)
        Eigen::Matrix<double, 3, 4> g_T;
        eMatrix(dr_init, g_T);

        // change of variables gradient tensor elements
        m_chi.block<3, 3>(body_id_7, particle_id_7).noalias() = s_part * m_I3; // translation-translation couple
        m_chi.block<3, 3>(body_id_7 + 4, particle_id_7).noalias() =
            m_psi_conv_quat_ang.block<3, 4>(body_id_7 + 3, body_id_7 + 3) *
            g_T.transpose(); // quaternion-rotation couple (first row is zero)
    }

    m_tens_chi = TensorCast(m_chi);

    /* ANCHOR : Compute m_tens_grad_zeta */
    m_tens_grad_zeta.setZero();

    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        /* ANCHOR: Tensor indices */
        const int particle_id_7{7 * particle_id};
        const int body_id_7{7 * m_particle_group_id(particle_id)};

        const bool is_locater{m_particle_group_id(particle_id) == 1};

        // `Eigen::Tensor` contract indices
        const Eigen::array<Eigen::IndexPair<int>, 1> contract_ilk_lj = {
            Eigen::IndexPair<int>(1, 0)}; // {i, l, k} . {l, j} --> {i, k, j}
        const Eigen::array<Eigen::IndexPair<int>, 1> contract_li_ljk = {
            Eigen::IndexPair<int>(0, 0)}; // {l, i} . {l, j, k} --> {i, j, k}
        const Eigen::array<Eigen::IndexPair<int>, 1> contract_ljm_km = {
            Eigen::IndexPair<int>(2, 1)}; // {l, j, m} . {k, m} --> {l, j, k}

        // `Eigen::Tensor` permute indices
        const Eigen::array<int, 3> permute_ikj_ijk({0, 2, 1}); // {i, k, j} --> {i, j, k}

        // `Eigen::Tensor` output tensor indices start location (offset) and extent
        const Eigen::array<Eigen::Index, 3> offsets_left  = {body_id_7 + 3, particle_id_7, body_id_7};
        const Eigen::array<Eigen::Index, 3> offsets_right = {body_id_7 + 3, particle_id_7 + 3, body_id_7};
        const Eigen::array<Eigen::Index, 3> extents       = {4, 3, 7};

        // get moment arm to locater point from particle
        const Eigen::Matrix3d two_r_cross_mat = 2 * m_rbm_conn.block<3, 3>(body_id_7 + 3, particle_id_7);
        const Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3>> tens_two_r_cross_mat = TensorCast(two_r_cross_mat);

        // get E matrix representation of body quaternion from m_psi_conv_quat_ang
        const Eigen::Matrix<double, 3, 4> two_E_body = m_psi_conv_quat_ang.block<3, 4>(body_id_7 + 3, body_id_7 + 3);
        const Eigen::TensorFixedSize<double, Eigen::Sizes<3, 4>> tens_two_E_body = TensorCast(two_E_body);

        // get change of variable matrix element m_chi_tilde
        const Eigen::Matrix<double, 7, 3> n_chi_tilde = -m_chi.block<7, 3>(body_id_7, particle_id_7);
        const Eigen::TensorFixedSize<double, Eigen::Sizes<7, 3>> tens_n_chi_tilde = TensorCast(n_chi_tilde);

        /* ANCHOR: tensor contractions for left-half of gradient */
        // compute grad_r_cross{l, j, k} = - m_levi_cevita{l, j, m} chi_tilde{k, m}
        Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3, 7>> grad_r_cross; // (3, 3, 7)  {l, j, k}
        grad_r_cross.device(device) = m_levi_cevita.contract(tens_n_chi_tilde, contract_ljm_km);

        // compute 2E_T_grad_r{i, j, k} = 2E_{l, i} grad_r_cross{l, j, k}
        Eigen::TensorFixedSize<double, Eigen::Sizes<4, 3, 7>> two_E_grad_r_cross; // (4, 3, 7)  {i, j, k}
        two_E_grad_r_cross.device(device) = tens_two_E_body.contract(grad_r_cross, contract_li_ljk);

        // compute 2grad_E_r_cross_{i, k, j} = Kappa{i, l, k} r_cross{l, j}
        Eigen::TensorFixedSize<double, Eigen::Sizes<4, 7, 3>> two_grad_E_r_cross_preshuffle; // (4, 3, 7)  {i, k, j}
        two_grad_E_r_cross_preshuffle.device(device) = m_kappa.contract(tens_two_r_cross_mat, contract_ilk_lj);

        // shuffle two_grad_E_r_cross_preshuffle
        Eigen::TensorFixedSize<double, Eigen::Sizes<4, 3, 7>> two_grad_E_r_cross; // (4, 3, 7)  {i, j, k}
        two_grad_E_r_cross.device(device) = two_grad_E_r_cross_preshuffle.shuffle(permute_ikj_ijk);

        // left term
        const Eigen::TensorFixedSize<double, Eigen::Sizes<4, 3, 7>> mixed_gradient =
            two_grad_E_r_cross + two_E_grad_r_cross;
        m_tens_grad_zeta.slice(offsets_left, extents) = mixed_gradient;

        /* ANCHOR: tensor contractions for right-half of gradient */
        m_tens_grad_zeta.slice(offsets_right, extents) = 2 * m_kappa;
    }
}

void
systemData::convertBody2ParticleOrient()
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
systemData::convertBody2ParticlePos()
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
systemData::convertBody2ParticleVelAcc(Eigen::ThreadPoolDevice& device)
{
    /* ANCHOR: Convert velocity D.o.F. */
    m_velocities_particles.noalias() = m_zeta.transpose() * m_velocities_bodies;
    m_velocities_particles.noalias() += m_velocities_particles_articulation;

    /* ANCHOR: Convert acceleration D.o.F. */

    // Eigen::Tensor contraction indices
    const Eigen::array<Eigen::IndexPair<int>, 0> outer_product   = {}; // no contractions
    const Eigen::array<Eigen::IndexPair<int>, 2> contract_jik_jk = {
        Eigen::IndexPair<int>(0, 0), Eigen::IndexPair<int>(2, 1)}; // {j i k}, {j k} --> {i}

    const Eigen::Tensor<double, 1> xi    = TensorCast(m_velocities_bodies); // (7M x 1)
    const Eigen::Tensor<double, 2> xi_xi = xi.contract(xi, outer_product);  // (7M x 7M)

    Eigen::Tensor<double, 1> velocities_rbm_particles = Eigen::Tensor<double, 1>(7 * m_num_particles); // (7N x 1)
    velocities_rbm_particles.device(device)           = m_tens_grad_zeta.contract(xi_xi, contract_jik_jk);

    m_accelerations_bodies.noalias() = MatrixCast(velocities_rbm_particles, 7 * m_num_particles, 1, device); // (7N x 1)
    m_accelerations_bodies.noalias() += m_accelerations_particles_articulation;                              // (7N x 1)
}