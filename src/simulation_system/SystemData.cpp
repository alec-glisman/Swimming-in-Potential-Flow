//
// Created by Alec Glisman on 07/30/21
//

#include <SystemData.hpp>

SystemData::SystemData(std::string inputGSDFile, std::string outputDir)
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

SystemData::~SystemData()
{
    spdlog::get(m_logName)->info("SystemData destructor called");
    gsd_close(m_handle.get());
    spdlog::get(m_logName)->flush();
    spdlog::drop(m_logName);
}

void
SystemData::initializeData()
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
    m_rbm_conn      = Eigen::MatrixXd::Zero(m7, n7);
    m_tens_rbm_conn = Eigen::Tensor<double, 2>(m7, n7);
    m_tens_rbm_conn.setZero();

    // initialize gradient matrices
    m_chi      = Eigen::MatrixXd::Zero(m7, n3);
    m_tens_chi = Eigen::Tensor<double, 2>(m7, n3);
    m_tens_chi.setZero();
    m_tens_grad_rbm_conn = Eigen::Tensor<double, 3>(m7, n7, m7);
    m_tens_grad_rbm_conn.setZero();

    // Set initial configuration orientation
    spdlog::get(m_logName)->info("Setting initial configuration orientation");
    m_positions_particles_articulation_init_norm.setZero();

    /// @review_swimmer: change initial orientations
    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        const int particle_id_3{3 * particle_id};

        double orient_dir{0.0};

        if (particle_id % 3 == 1)
        {
            orient_dir = 1.0; // +x placement
        }
        else if (particle_id % 3 == 2)
        {
            orient_dir = -1.0; // -x placement
        }

        m_positions_particles_articulation_init_norm(particle_id_3) =
            orient_dir; // REVIEW: All of body coordinate in +x-axis, quaternion rotates

        spdlog::get(m_logName)->info("Particle {0} initial orientation: [{1:03.14f}, {2:03.14f}, {3:03.14f}]",
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
SystemData::parseGSD()
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
SystemData::checkInput()
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
SystemData::logData()
{
    /* ANCHOR: Output simulation data */
    spdlog::get(m_logName)->info("Starting logdata()");
    spdlog::get(m_logName)->info("time: {0}", m_t);

    /* ANCHOR: Output body data */
    spdlog::get(m_logName)->info("Body positions:");
    for (int body_id = 0; body_id < m_num_bodies; body_id++)
    {
        const int body_id_7{7 * body_id};
        spdlog::get(m_logName)->info(
            "\tBody {0}: [{1:03.14f}, {2:03.14f}, {3:03.14f}, {4:03.14f}, {5:03.14f}, {6:03.14f}, {7:03.14f}]",
            body_id + 1, m_positions_bodies(body_id_7), m_positions_bodies(body_id_7 + 1),
            m_positions_bodies(body_id_7 + 2), m_positions_bodies(body_id_7 + 3), m_positions_bodies(body_id_7 + 4),
            m_positions_bodies(body_id_7 + 5), m_positions_bodies(body_id_7 + 6));
    }
    spdlog::get(m_logName)->info("Body velocities:");
    for (int body_id = 0; body_id < m_num_bodies; body_id++)
    {
        const int body_id_7{7 * body_id};
        spdlog::get(m_logName)->info(
            "\tBody {0}: [{1:03.14f}, {2:03.14f}, {3:03.14f}, {4:03.14f}, {5:03.14f}, {6:03.14f}, {7:03.14f}]",
            body_id + 1, m_velocities_bodies(body_id_7), m_velocities_bodies(body_id_7 + 1),
            m_velocities_bodies(body_id_7 + 2), m_velocities_bodies(body_id_7 + 3), m_velocities_bodies(body_id_7 + 4),
            m_velocities_bodies(body_id_7 + 5), m_velocities_bodies(body_id_7 + 6));
    }
    spdlog::get(m_logName)->info("Body accelerations:");
    for (int body_id = 0; body_id < m_num_bodies; body_id++)
    {
        const int body_id_7{7 * body_id};
        spdlog::get(m_logName)->info(
            "\tBody {0}: [{1:03.14f}, {2:03.14f}, {3:03.14f}, {4:03.14f}, {5:03.14f}, {6:03.14f}, {7:03.14f}]",
            body_id + 1, m_accelerations_bodies(body_id_7), m_accelerations_bodies(body_id_7 + 1),
            m_accelerations_bodies(body_id_7 + 2), m_accelerations_bodies(body_id_7 + 3),
            m_accelerations_bodies(body_id_7 + 4), m_accelerations_bodies(body_id_7 + 5),
            m_accelerations_bodies(body_id_7 + 6));
    }

    /* ANCHOR: Output particle data */
    spdlog::get(m_logName)->info("Particle orientations:");
    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        const int particle_id_3{3 * particle_id};
        spdlog::get(m_logName)->info("\tParticle {0}: [{1:03.14f}, {2:03.14f}, {3:03.14f}]", particle_id + 1,
                                     m_orientations_particles(particle_id_3),
                                     m_orientations_particles(particle_id_3 + 1),
                                     m_orientations_particles(particle_id_3 + 2));
    }

    spdlog::get(m_logName)->info("Particle positions:");
    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        const int particle_id_3{3 * particle_id};
        spdlog::get(m_logName)->info("\tParticle {0}: [{1:03.14f}, {2:03.14f}, {3:03.14f}]", particle_id + 1,
                                     m_positions_particles(particle_id_3), m_positions_particles(particle_id_3 + 1),
                                     m_positions_particles(particle_id_3 + 2));
    }

    spdlog::get(m_logName)->info("Particle velocities:");
    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        const int particle_id_7{7 * particle_id};
        spdlog::get(m_logName)->info(
            "\tParticle {0}: [{1:03.14f}, {2:03.14f}, {3:03.14f}]", particle_id + 1,
            m_velocities_particles(particle_id_7), m_velocities_particles(particle_id_7 + 1),
            m_velocities_particles(particle_id_7 + 2), m_velocities_particles(particle_id_7 + 3),
            m_velocities_particles(particle_id_7 + 4), m_velocities_particles(particle_id_7 + 5),
            m_velocities_particles(particle_id_7 + 6));
    }

    spdlog::get(m_logName)->info("Particle accelerations:");
    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        const int particle_id_7{7 * particle_id};
        spdlog::get(m_logName)->info(
            "\tParticle {0}: [{1:03.14f}, {2:03.14f}, {3:03.14f}]", particle_id + 1,
            m_accelerations_particles(particle_id_7), m_accelerations_particles(particle_id_7 + 1),
            m_accelerations_particles(particle_id_7 + 2), m_accelerations_particles(particle_id_7 + 3),
            m_accelerations_particles(particle_id_7 + 4), m_accelerations_particles(particle_id_7 + 5),
            m_accelerations_particles(particle_id_7 + 6));
    }

    /* ANCHOR: Output particle *articulation* data */
    spdlog::get(m_logName)->info("Particle articulation positions:");
    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        const int particle_id_3{3 * particle_id};
        spdlog::get(m_logName)->info("\tParticle {0}: [{1:03.14f}, {2:03.14f}, {3:03.14f}]", particle_id + 1,
                                     m_positions_particles_articulation(particle_id_3),
                                     m_positions_particles_articulation(particle_id_3 + 1),
                                     m_positions_particles_articulation(particle_id_3 + 2));
    }

    spdlog::get(m_logName)->info("Particle articulation velocities:");
    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        const int particle_id_7{7 * particle_id};
        spdlog::get(m_logName)->info(
            "\tParticle {0}: [{1:03.14f}, {2:03.14f}, {3:03.14f}, {4:03.14f}, {5:03.14f}, {6:03.14f}, {7:03.14f}]",
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
            "\tParticle {0}: [{1:03.14f}, {2:03.14f}, {3:03.14f}, {4:03.14f}, {5:03.14f}, {6:03.14f}, {7:03.14f}]",
            particle_id + 1, m_accelerations_particles_articulation(particle_id_7),
            m_accelerations_particles_articulation(particle_id_7 + 1),
            m_accelerations_particles_articulation(particle_id_7 + 2),
            m_accelerations_particles_articulation(particle_id_7 + 3),
            m_accelerations_particles_articulation(particle_id_7 + 4),
            m_accelerations_particles_articulation(particle_id_7 + 5),
            m_accelerations_particles_articulation(particle_id_7 + 6));
    }

    spdlog::get(m_logName)->info("Ending logdata()");
    spdlog::get(m_logName)->flush();
}

void
SystemData::update(const Eigen::ThreadPoolDevice& device)
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
SystemData::normalizeQuaternions()
{
    for (int body_id = 0; body_id < m_num_bodies; body_id++)
    {
        const int body_id_7{7 * body_id};

        const double quaternion_norm{m_positions_bodies.segment<4>(body_id_7 + 3).norm()};

        m_positions_bodies.segment<4>(body_id_7 + 3) /= quaternion_norm;
    }
}

void
SystemData::positionsArticulation()
{

    // articulation acceleration magnitudes
    const double t_dimensional{m_tau * m_t};
    const double r1_mag = m_sys_spec_R_avg + (m_sys_spec_U0 / m_sys_spec_omega) * sin(m_sys_spec_omega * t_dimensional);
    const double r3_mag = m_sys_spec_R_avg - (m_sys_spec_U0 / m_sys_spec_omega) *
                                                 sin(m_sys_spec_omega * t_dimensional + m_sys_spec_phase_shift);

    int particle_id{0};
    m_positions_particles_articulation.setZero();

    for (int body_id = 0; body_id < m_num_bodies; body_id++)
    {
        const int constrained1_3{3 * (particle_id + 1)};
        const int constrained2_3{3 * (particle_id + 2)};

        m_positions_particles_articulation.segment<3>(constrained1_3).noalias() =
            r1_mag * m_orientations_particles.segment<3>(constrained1_3);

        m_positions_particles_articulation.segment<3>(constrained2_3).noalias() =
            r3_mag * m_orientations_particles.segment<3>(constrained2_3);

        particle_id += 3;
    }
}

void
SystemData::velocitiesArticulation()
{

    // articulation velocity magnitudes
    const double t_dimensional{m_tau * m_t};
    const double v1_mag = m_sys_spec_U0 * cos(m_sys_spec_omega * t_dimensional);
    const double v3_mag = -m_sys_spec_U0 * cos(m_sys_spec_omega * t_dimensional + m_sys_spec_phase_shift);

    // Zero and then calculate  m_velocities_particles_articulation
    int particle_id{0};
    m_velocities_particles_articulation.setZero();

    for (int body_id = 0; body_id < m_num_bodies; body_id++)
    {
        const int constrained1_3{3 * (particle_id + 1)};
        const int constrained1_7{7 * (particle_id + 1)};

        const int constrained2_3{3 * (particle_id + 2)};
        const int constrained2_7{7 * (particle_id + 2)};

        m_velocities_particles_articulation.segment<3>(constrained1_7).noalias() =
            v1_mag * m_orientations_particles.segment<3>(constrained1_3);

        m_velocities_particles_articulation.segment<3>(constrained2_7).noalias() =
            v3_mag * m_orientations_particles.segment<3>(constrained2_3);

        // Image system: flip z components, leave x,y unchanged
        if (particle_id >= (m_num_particles / 2))
        {
            m_velocities_particles_articulation(constrained1_7 + 2) *= -1;
            m_velocities_particles_articulation(constrained2_7 + 2) *= -1;
        }

        particle_id += 3;
    }
}

void
SystemData::accelerationsArticulation()
{

    // articulation acceleration magnitudes
    const double t_dimensional{m_tau * m_t};
    const double a1_mag = -m_sys_spec_U0 * m_sys_spec_omega * sin(m_sys_spec_omega * t_dimensional);
    const double a3_mag =
        m_sys_spec_U0 * m_sys_spec_omega * sin(m_sys_spec_omega * t_dimensional + m_sys_spec_phase_shift);

    // Zero and then calculate  m_accelerations_particles_articulation
    int particle_id{0};
    m_accelerations_particles_articulation.setZero();

    for (int body_id = 0; body_id < m_num_bodies; body_id++)
    {
        const int constrained1_3{3 * (particle_id + 1)};
        const int constrained1_7{7 * (particle_id + 1)};

        const int constrained2_3{3 * (particle_id + 2)};
        const int constrained2_7{7 * (particle_id + 2)};

        m_accelerations_particles_articulation.segment<3>(constrained1_7).noalias() =
            a1_mag * m_orientations_particles.segment<3>(constrained1_3);

        m_accelerations_particles_articulation.segment<3>(constrained2_7).noalias() =
            a3_mag * m_orientations_particles.segment<3>(constrained2_3);

        // Image system: flip z components, leave x,y unchanged
        if (particle_id >= (m_num_particles / 2))
        {
            m_accelerations_particles_articulation(constrained1_7 + 2) *= -1;
            m_accelerations_particles_articulation(constrained2_7 + 2) *= -1;
        }

        particle_id += 3;
    }
}

void
SystemData::udwadiaLinearSystem()
{
    m_Udwadia_A.setZero();
    m_Udwadia_b.setZero();

    for (int body_id = 0; body_id < m_num_constraints; body_id++)
    {
        const int body_id_7{7 * body_id};
        const int quat_start{body_id_7 + 3};

        const Eigen::Vector4d quat = m_positions_bodies.segment<4>(quat_start);
        const double          d_quat_norm_sqr{m_velocities_bodies.segment<4>(quat_start).squaredNorm()};

        m_Udwadia_A.row(body_id).segment<4>(quat_start) = quat;
        m_Udwadia_b(body_id)                            = -d_quat_norm_sqr;
    }
}

void
SystemData::rbmMatrixElement(const int particle_id)
{
    const int particle_id_3{3 * particle_id};
    const int particle_id_7{7 * particle_id};
    const int body_id_7{7 * m_particle_group_id(particle_id)};

    // \[r_\alpha ^\]: skew-symmetric matrix representation of cross product
    Eigen::Matrix3d mat_two_dr_cross;
    crossProdMat(m_positions_particles_articulation.segment<3>(particle_id_3), mat_two_dr_cross);
    // preprend row of zeros
    Eigen::Matrix<double, 4, 3> mat_dr_cross_43;
    mat_dr_cross_43.setZero();
    mat_dr_cross_43.block<3, 3>(1, 0).noalias() = mat_two_dr_cross;

    // E_{(i)}: matrix representation of left quaternion composition
    Eigen::Matrix4d E_body;
    eMatrix(m_positions_bodies.segment<4>(body_id_7 + 3), E_body);
    const Eigen::Matrix4d two_E_T = 2.0 * E_body.transpose();

    // rigid body motion connectivity tensor elements
    m_rbm_conn.block<3, 3>(body_id_7, particle_id_7).noalias() = m_I3; // translation-translation couple

    /// @FIXME: (check literature) The teaching SD to swim paper had a negative sign here, but my derivations do not
    /// have that
    m_rbm_conn.block<4, 3>(body_id_7 + 3, particle_id_7).noalias() =
        two_E_T * mat_dr_cross_43; // translation-quaternion couple

    m_rbm_conn.block<4, 4>(body_id_7 + 3, particle_id_7 + 3).noalias() = two_E_T; // quaternion-quaternion couple
}

void
SystemData::chiMatrixElement(const int particle_id)
{
    const int particle_id_3{3 * particle_id};                  // particle number
    const int body_id_7{7 * m_particle_group_id(particle_id)}; // body number

    // Scaling prefactor: 2 * || r_particle_id ||
    const double prefactor{2.0 * m_positions_particles_articulation.segment<3>(particle_id_3).norm()};

    // unit quaternion basis directions
    const Eigen::Quaterniond q_i(0.0, 1.0, 0.0, 0.0);
    const Eigen::Quaterniond q_j(0.0, 0.0, 1.0, 0.0);
    const Eigen::Quaterniond q_k(0.0, 0.0, 0.0, 1.0);

    /* ANCHOR: Compute G matrix element */
    // body unit quaternion
    const Eigen::Quaterniond theta_body(m_positions_bodies(body_id_7 + 3), m_positions_bodies(body_id_7 + 4),
                                        m_positions_bodies(body_id_7 + 5), m_positions_bodies(body_id_7 + 6));

    // particle initial configuration unit quaternion
    const Eigen::Quaterniond r_body_quat(0.0, prefactor * m_positions_particles_articulation_init_norm(particle_id_3),
                                         prefactor * m_positions_particles_articulation_init_norm(particle_id_3 + 1),
                                         prefactor * m_positions_particles_articulation_init_norm(particle_id_3 + 2));

    // quaternion product: r_body * theta
    const Eigen::Quaterniond r_theta = r_body_quat * theta_body;
    // quaternion product: r_theta * (unit quaternion basis directions)
    const Eigen::Quaterniond r_theta_i = r_theta * q_i;
    const Eigen::Quaterniond r_theta_j = r_theta * q_j;
    const Eigen::Quaterniond r_theta_k = r_theta * q_k;

    // G matrix element
    Eigen::Matrix<double, 4, 3> g_matrix = Eigen::Matrix<double, 4, 3>::Zero(4, 3);
    g_matrix.row(0).noalias()            = r_theta.vec();
    g_matrix.row(1).noalias()            = -r_theta_i.vec();
    g_matrix.row(2).noalias()            = -r_theta_j.vec();
    g_matrix.row(3).noalias()            = -r_theta_k.vec();

    /* ANCHOR: Compute chi matrix element */
    m_chi.block<3, 3>(body_id_7, particle_id_3).noalias() =
        m_I3; // convert body (linear) position derivatives to particle linear coordinate derivatives
    m_chi.block<4, 3>(body_id_7 + 3, particle_id_3).noalias() =
        g_matrix; // convert body (quaternion) position derivatives to particle linear coordinate derivatives
}

void
SystemData::gradRbmConnTensorElement(const int particle_id, const Eigen::ThreadPoolDevice& device)
{
    /* ANCHOR: Tensor indices */
    const int  particle_id_3{3 * particle_id};
    const int  particle_id_7{7 * particle_id};
    const int  body_id_7{7 * m_particle_group_id(particle_id)};
    const bool is_locater{m_particle_type_id(particle_id) == 1};

    // `Eigen::Tensor` contract indices
    const Eigen::array<Eigen::IndexPair<int>, 1> contract_ilk_lj = {
        Eigen::IndexPair<int>(1, 0)}; // {i, l, k} . {l, j} --> {i, k, j}
    const Eigen::array<Eigen::IndexPair<int>, 1> contract_il_ljk = {
        Eigen::IndexPair<int>(1, 0)}; // {l, i} . {l, j, k} --> {i, j, k}
    const Eigen::array<Eigen::IndexPair<int>, 1> contract_ljm_km = {
        Eigen::IndexPair<int>(2, 1)}; // {l, j, m} . {k, m} --> {l, j, k}

    // `Eigen::Tensor` permute indices
    const Eigen::array<int, 3> permute_ikj_ijk({0, 2, 1}); // {i, k, j} --> {i, j, k}

    // `Eigen::Tensor` output tensor indices start location (offset) and extent
    const Eigen::array<Eigen::Index, 3> offset_grad_Rc_r_cross  = {1, 0, 0};
    const Eigen::array<Eigen::Index, 3> extents_grad_Rc_r_cross = {3, 3, 3};

    const Eigen::array<Eigen::Index, 3> offset_grad_theta_r_cross  = {1, 0, 3};
    const Eigen::array<Eigen::Index, 3> extents_grad_theta_r_cross = {3, 3, 4};

    const Eigen::array<Eigen::Index, 3> offsets_mixed = {body_id_7 + 3, particle_id_7, body_id_7};
    const Eigen::array<Eigen::Index, 3> extents_mixed = {4, 3, 7};

    const Eigen::array<Eigen::Index, 3> offsets_angular = {body_id_7 + 3, particle_id_7 + 3, body_id_7};
    const Eigen::array<Eigen::Index, 3> extents_angular = {4, 4, 7};

    /* ANCHOR: Tensor quantities that will be contracted */
    // moment arm to locater point from particle
    Eigen::Matrix3d mat_two_dr_cross;
    crossProdMat(2.0 * m_positions_particles_articulation.segment<3>(particle_id_3), mat_two_dr_cross);
    // preprend row of zeros
    Eigen::Matrix<double, 4, 3> two_r_tilde_cross_mat = Eigen::MatrixXd::Zero(4, 3);
    two_r_tilde_cross_mat.block<3, 3>(1, 0).noalias() = mat_two_dr_cross;

    const Eigen::TensorFixedSize<double, Eigen::Sizes<4, 3>> tens_two_r_tilde_cross_mat =
        TensorCast(two_r_tilde_cross_mat, 4, 3);

    // E_{(i)}: matrix representation of body quaternion from m_rbm_conn
    const Eigen::TensorFixedSize<double, Eigen::Sizes<4, 4>> tens_two_ET_body =
        TensorCast(m_rbm_conn.block<4, 4>(body_id_7 + 3, particle_id_7 + 3), 4, 4);

    // G_{(i) a}: Jacobian matrix from particle positions to body quaternion
    const Eigen::TensorFixedSize<double, Eigen::Sizes<4, 3>> tens_g_ia =
        TensorCast(m_chi.block<4, 3>(body_id_7 + 3, particle_id_3), 4, 3); // (4, 3)

    /* ANCHOR: tensor contractions for left-half of gradient */
    // compute grad_r_cross{l, j, k}
    Eigen::TensorFixedSize<double, Eigen::Sizes<4, 3, 7>> grad_r_cross; // (4, 3, 7)  {l, j, k}
    grad_r_cross.setZero();

    if (!is_locater)
    {
        grad_r_cross.slice(offset_grad_Rc_r_cross, extents_grad_Rc_r_cross).device(device) =
            m_levi_cevita; // body locater point gradient
        grad_r_cross.slice(offset_grad_theta_r_cross, extents_grad_theta_r_cross).device(device) =
            -m_levi_cevita.contract(tens_g_ia, contract_ljm_km); // body unit quaternion gradient
    }

    // compute 2ET_grad_r{i, j, k} = 2ET_{i, l} grad_r_cross{l, j, k}
    Eigen::TensorFixedSize<double, Eigen::Sizes<4, 3, 7>> two_ET_grad_r_cross; // (4, 3, 7)  {i, j, k}
    two_ET_grad_r_cross.device(device) = tens_two_ET_body.contract(grad_r_cross, contract_il_ljk);

    // compute 2grad_E_r_cross_{i, k, j} = Kappa{i, l, k} r_cross{l, j}
    Eigen::TensorFixedSize<double, Eigen::Sizes<4, 7, 3>> two_grad_ET_r_cross_preshuffle; // (4, 7, 3)  {i, k, j}
    two_grad_ET_r_cross_preshuffle.device(device) = m_kappa.contract(tens_two_r_tilde_cross_mat, contract_ilk_lj);
    // shuffle two_grad_ET_r_cross_preshuffle
    Eigen::TensorFixedSize<double, Eigen::Sizes<4, 3, 7>> two_grad_ET_r_cross_tilde; // (4, 3, 7)  {i, j, k}
    two_grad_ET_r_cross_tilde.device(device) = two_grad_ET_r_cross_preshuffle.shuffle(permute_ikj_ijk);

    // mixed gradient term
    Eigen::TensorFixedSize<double, Eigen::Sizes<4, 3, 7>> mixed_gradient; // (4, 3, 7)  {i, j, k}
    mixed_gradient.device(device) = two_grad_ET_r_cross_tilde;
    mixed_gradient.device(device) += two_ET_grad_r_cross;
    m_tens_grad_rbm_conn.slice(offsets_mixed, extents_mixed).device(device) = mixed_gradient;

    /* ANCHOR: tensor contractions for quaternion-quaternion of gradient term */
    m_tens_grad_rbm_conn.slice(offsets_angular, extents_angular).device(device) = 2.0 * m_kappa;
}

void
SystemData::rigidBodyMotionTensors(const Eigen::ThreadPoolDevice& device)
{
    /* ANCHOR: Compute m_rbm_conn */
    m_rbm_conn.setZero();

    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        rbmMatrixElement(particle_id);
    }

    m_tens_rbm_conn = TensorCast(m_rbm_conn, 7 * m_num_bodies, 7 * m_num_particles);
}

void
SystemData::gradientChangeOfVariableTensors(const Eigen::ThreadPoolDevice& device)
{
    /* ANCHOR: Compute m_chi and m_tens_chi */
    m_chi.setZero();

    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        chiMatrixElement(particle_id);
    }

    m_tens_chi = TensorCast(m_chi, 7 * m_num_bodies, 3 * m_num_particles);

    /* ANCHOR : Compute m_tens_grad_rbm_conn */
    m_tens_grad_rbm_conn.setZero();

    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        gradRbmConnTensorElement(particle_id, device);
    }
}

void
SystemData::convertBody2ParticleOrient()
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

        // body unit quaternion
        const Eigen::Vector4d    theta = m_positions_bodies.segment<4>(body_id_7 + 3);
        const Eigen::Quaterniond theta_body(theta(0), theta(1), theta(2), theta(3));

        // particle initial configuration unit quaternion
        const Eigen::Vector3d    r_hat_init = m_positions_particles_articulation_init_norm.segment<3>(particle_id_3);
        const Eigen::Quaterniond r_body_quat(0.0, r_hat_init(0), r_hat_init(1), r_hat_init(2));

        const Eigen::Quaterniond orient_rot    = theta_body * r_body_quat * theta_body.inverse();
        const Eigen::Vector3d    orient_rot_3d = orient_rot.vec();

        // Output rotated orientation
        m_orientations_particles.segment<3>(particle_id_3).noalias() = orient_rot_3d;

#if !defined(NDEBUG)
        Eigen::IOFormat CleanFmt(16, 0, ", ", "\n", "[", "]");

        const double epsilon_zero{1e-12}; // "small value" close to zero for double comparison
        const double epsilon_norm{1e-6};  // "small value" close to zero for double comparison

        const double particle_quat_norm{theta_body.norm()};
        const double unit_norm_error_mag{particle_quat_norm - 1.0};

        const bool zeroith_comp_cond{abs(orient_rot.w()) >= epsilon_zero};
        const bool unit_norm_cond{abs(unit_norm_error_mag) >= epsilon_norm};

        if (zeroith_comp_cond)
        {
            const std::string zeroith_comp_msg{
                "At t = " + std::to_string(m_t) + ", 0th component of rotated particle # " +
                std::to_string(particle_quat_norm) + " vector is " + std::to_string(orient_rot.w())};

            throw std::logic_error(zeroith_comp_msg);
        }
        else if (unit_norm_cond)
        {
            std::ostringstream stream_error_mag;
            stream_error_mag << unit_norm_error_mag;

            const std::string unit_norm_msg{
                "Quaternion is not unitary at t = " + std::to_string(m_t) +
                ".\n\tParticle #: " + std::to_string(particle_id) + ",\n\tNorm - 1: " + stream_error_mag.str() +
                ",\n\ttheta_body: " + std::to_string(theta_body.w()) + ", " + std::to_string(theta_body.x()) + ", " +
                std::to_string(theta_body.y()) + ", " + std::to_string(theta_body.z())};

            throw std::logic_error(unit_norm_msg);
        }
#endif
    }
}

void
SystemData::convertBody2ParticlePos()
{
    for (int particle_id = 0; particle_id < m_num_particles; particle_id++)
    {
        const int particle_id_3{3 * particle_id};
        const int body_id_7{7 * m_particle_group_id(particle_id)};

        m_positions_particles.segment<3>(particle_id_3).noalias() = m_positions_bodies.segment<3>(body_id_7);
        m_positions_particles.segment<3>(particle_id_3).noalias() +=
            m_positions_particles_articulation.segment<3>(particle_id_3);
    }
}

void
SystemData::convertBody2ParticleVelAcc(const Eigen::ThreadPoolDevice& device)
{
    /* ANCHOR: Convert velocity D.o.F. */
    m_velocities_particles.noalias() = m_rbm_conn.transpose() * m_velocities_bodies;
    m_velocities_particles.noalias() += m_velocities_particles_articulation;

    /* ANCHOR: Convert acceleration D.o.F. */

    // Eigen::Tensor contraction indices
    const Eigen::array<Eigen::IndexPair<int>, 0> outer_product   = {}; // no contractions
    const Eigen::array<Eigen::IndexPair<int>, 2> contract_jik_jk = {
        Eigen::IndexPair<int>(0, 0), Eigen::IndexPair<int>(2, 1)}; // {j i k}, {j k} --> {i}

    const Eigen::Tensor<double, 1> xi    = TensorCast(m_velocities_bodies, 7 * m_num_bodies); // (7M x 1)
    const Eigen::Tensor<double, 2> xi_xi = xi.contract(xi, outer_product);                    // (7M x 7M)

    Eigen::Tensor<double, 1> velocities_rbm_particles = Eigen::Tensor<double, 1>(7 * m_num_particles); // (7N x 1)
    velocities_rbm_particles.device(device)           = m_tens_grad_rbm_conn.contract(xi_xi, contract_jik_jk);
    const Eigen::VectorXd mat_vel_rbm_particles = MatrixCast(velocities_rbm_particles, 7 * m_num_particles, 1, device);

    m_accelerations_particles.noalias() = mat_vel_rbm_particles;                            // (7N x 1)
    m_accelerations_particles.noalias() += m_rbm_conn.transpose() * m_accelerations_bodies; // (7N x 1)
    m_accelerations_particles.noalias() += m_accelerations_particles_articulation;          // (7N x 1)
}
