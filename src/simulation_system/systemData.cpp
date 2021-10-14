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
    m_num_DoF         = 6 * m_num_bodies; // D.o.F. are linear and angular positions of body centers
    m_num_constraints = m_num_bodies;     // 1 unit quaternion constraint per body

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

    // initialize kinematic vectors
    m_positions_locater_particles          = Eigen::VectorXd::Zero(3 * m_num_bodies);
    m_velocities_particles_articulation    = Eigen::VectorXd::Zero(3 * m_num_particles);
    m_accelerations_particles_articulation = Eigen::VectorXd::Zero(3 * m_num_particles);

    // initialize constraint matrices
    m_Udwadia_A = Eigen::MatrixXd::Zero(m_num_constraints, m_num_DoF);
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

    // initialize constraints
    spdlog::get(m_logName)->critical("Setting up temporary single-thread eigen device");
    Eigen::ThreadPool       thread_pool = Eigen::ThreadPool(1);
    Eigen::ThreadPoolDevice single_core_device(&thread_pool, 1);

    spdlog::get(m_logName)->info("Initializing constraints");
    update(single_core_device);
}

void
systemData::update(Eigen::ThreadPoolDevice& device)
{
    // NOTE: Ordering of following functions does not matter
    velocitiesArticulation();
    accelerationsArticulation();

    locaterPointLocations();
    rigidBodyMotionTensors(device);

    udwadiaLinearSystem();

    // NOTE: Call gradientChangeOfVariableTensors() after rigidBodyMotionTensors()
    gradientChangeOfVariableTensors(device);
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

    assert(m_orientations_particles.size() == 4 * m_num_particles &&
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

// FIXME: Redo in light of new particle ordering
/* REVIEW[epic=Change,order=3]: Change assignment of m_velocities_particles_articulation for
 * different systems */
void
systemData::velocitiesArticulation()
{
    double dimensional_time{m_tau * m_t};

    // ANCHOR: Orientation vectors, q = R_1 - R_3
    Eigen::Vector3d q = m_positions_particles.segment<3>(3 * 0) - m_positions_particles.segment<3>(3 * 2);
    q.normalize();
    Eigen::Vector3d q_tilde = q;
    q_tilde(2) *= -1;

    // articulation velocity magnitudes
    const double v1_mag = m_sys_spec_U0 * cos(m_sys_spec_omega * dimensional_time);
    const double v3_mag = m_sys_spec_U0 * cos(m_sys_spec_omega * dimensional_time + m_sys_spec_phase_shift);

    // Zero and then calculate  m_velocities_particles_articulation
    m_velocities_particles_articulation.setZero();
    m_velocities_particles_articulation.segment<3>(3 * 0).noalias() = v1_mag * q;
    m_velocities_particles_articulation.segment<3>(3 * 2).noalias() = v3_mag * q;
    m_velocities_particles_articulation.segment<3>(3 * 3).noalias() = v1_mag * q_tilde;
    m_velocities_particles_articulation.segment<3>(3 * 5).noalias() = v3_mag * q_tilde;
}

// FIXME: Redo in light of new particle ordering
/* REVIEW[epic=Change,order=4]: Change assignment of m_accelerations_particles_articulation for
 * different systems */
void
systemData::accelerationsArticulation()
{
    double dimensional_time{m_tau * m_t};

    // ANCHOR: Orientation vectors, q = R_1 - R_3
    Eigen::Vector3d q = m_positions_particles.segment<3>(3 * 0) - m_positions_particles.segment<3>(3 * 2);
    q.normalize();
    Eigen::Vector3d q_tilde = q;
    q_tilde(2) *= -1;

    // articulation acceleration magnitudes
    const double a1_mag = -m_sys_spec_U0 * m_sys_spec_omega * sin(m_sys_spec_omega * dimensional_time);
    const double a3_mag =
        -m_sys_spec_U0 * m_sys_spec_omega * sin(m_sys_spec_omega * dimensional_time + m_sys_spec_phase_shift);

    // Zero and then calculate m_accelerations_particles_articulation
    m_accelerations_particles_articulation.setZero();
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
            m_positions_bodies.segment<3>(7 * body_count).noalias() = m_positions_particles.segment<3>(3 * i);

            body_count += 1;
        }
    }
    assert(body_count == m_num_bodies && "Incorrect number of bodies filled");
}

/* REVIEW[epic=Change,order=1]: Change udwadiaLinearSystem() for different systems
 */
void
systemData::udwadiaLinearSystem()
{
    m_Udwadia_A.setZero();
    m_Udwadia_b.setZero();

    // TODO: Implement the constraint linear system for the unit quaternions
}

void
systemData::rigidBodyMotionTensors(Eigen::ThreadPoolDevice& device)
{
    /* ANCHOR: Compute m_rbm_conn */
    m_rbm_conn.setZero();

    for (int i = 0; i < m_num_particles; i++)
    {
        const int             i3{3 * i};
        const Eigen::Vector3d R_i = m_positions_particles.segment<3>(i3);

        // Calculate which body j particle i belongs to
        const int body_j = (m_particle_type_id.segment(0, i + 1).array() == 1).count() - 1;
        const int j6{6 * body_j};
        const int j7{7 * body_j};

        // (negative) moment arm of particle about its body's locater position
        Eigen::Vector3d n_dr = m_positions_bodies.segment<3>(j7);
        n_dr.noalias() -= R_i;

        // skew-symmetric matrix representation of cross product
        Eigen::Matrix3d n_dr_cross;
        crossProdMat(n_dr, n_dr_cross);

        // rigid body motion connectivity tensor elements
        m_rbm_conn.block<3, 3>(j6, i3).noalias()     = m_I3;       // translation-translation couple
        m_rbm_conn.block<3, 3>(j6 + 3, i3).noalias() = n_dr_cross; // translation-rotation couple

        if (i == 0)
        {
            assert(body_j == 0 && "First particle must be part of body 0.");
        }
        if (i == m_num_particles - 1)
        {
            assert(body_j + 1 == m_num_bodies && "Last particle must be part of body M.");
        }
    }

    /* ANCHOR: Compute m_psi_conv_quat_ang */
    m_psi_conv_quat_ang.setZero();

    for (int k = 0; k < m_num_bodies; k++)
    {
        const int k6{6 * k};
        const int k7{7 * k};

        // matrix E from quaterion of body k
        Eigen::Matrix<double, 3, 4> E_theta_k;
        eMatrix(m_positions_bodies.segment<4>(k7 + 3), E_theta_k);

        // matrix elements of Psi
        m_psi_conv_quat_ang.block<3, 3>(k6, k7).noalias()         = m_I3; // no conversion from linear components
        m_psi_conv_quat_ang.block<3, 4>(k6 + 3, k7 + 3).noalias() = 2 * E_theta_k; // angular-quaternion velocity couple
    }

    /* ANCHOR: Compute m_rbm_conn_T_quat */
    m_rbm_conn_T_quat.noalias() = m_rbm_conn.transpose() * m_psi_conv_quat_ang;
    m_tens_rbm_conn_T_quat      = TensorCast(m_rbm_conn_T_quat);
}

void
systemData::gradientChangeOfVariableTensors(Eigen::ThreadPoolDevice& device)
{
    /* ANCHOR: Compute m_conv_body_2_part_dof */
    m_conv_body_2_part_dof.setZero();

    int                         body_num{-1}; // -1 as first particle should be locater particle and increment this
    Eigen::Matrix<double, 3, 4> twoE_body;    ///< 2 * E_body
    Eigen::Vector4d             R_c_body;     ///< R_c^{(i)}, locater position

    assert(m_particle_type_id(0) == 1 && "First particle index must be locater by convention");
    for (int j = 0; j < m_num_particles; j++)
    {

        if (m_particle_type_id(j) == 1)
        {
            // Zero out the relevent locater values
            twoE_body.setZero();
            R_c_body.setZero();
            // Increment body count
            body_num += 1;

            // Get locater position of body i
            R_c_body.segment<3>(1).noalias() = m_positions_locater_particles.segment<3>(3 * body_num);
            // Get E matrix representation of body quaternion from m_psi_conv_quat_ang
            twoE_body.noalias() = m_psi_conv_quat_ang.block<3, 4>(6 * body_num + 3, 7 * body_num + 3);

            // Continue to next loop as all elements are zero
            continue;
        }

        // j -> particle number
        const int j3{3 * j};

        // 4-vector version of particle j position (prepend zero element)
        Eigen::Vector4d R_j         = Eigen::Vector4d::Zero(4, 1);
        R_j.segment<3>(1).noalias() = m_positions_particles.segment<3>(j3);

        // matrix P_tilde from moment arm
        const Eigen::Vector4d       dr = R_j - R_c_body;
        Eigen::Matrix<double, 3, 4> P_j_tilde_hold;
        eMatrix(dr, P_j_tilde_hold);
        Eigen::Matrix<double, 4, 3> P_j_tilde = P_j_tilde_hold.transpose();

        // change of variables gradient tensor elements
        m_conv_body_2_part_dof.block<3, 3>(7 * body_num, j3).noalias() = -m_I3; // translation-translation couple
        m_conv_body_2_part_dof.block<3, 3>(7 * body_num + 4, j3).noalias() =
            twoE_body * P_j_tilde; // quaternion-rotation couple (first row is zero)
    }
    assert(body_num + 1 == m_num_bodies && "Not all bodies were indexed correctly");

    /* ANCHOR : Compute m_rbm_conn_T_quat_grad */
    m_rbm_conn_T_quat_grad.setZero();
    body_num = -1;

    assert(m_particle_type_id(0) == 1 && "First particle index must be locater by convention");
    for (int i = 0; i < m_num_particles; i++)
    {
        if (m_particle_type_id(i) == 1)
        {
            // Increment body count
            body_num += 1;

            // Continue to next loop as all elements are zero
            continue;
        }

        // get moment arm to locater point from C
        const Eigen::Matrix3d two_r_cross_mat = -2 * m_rbm_conn.block<3, 3>(6 * body_num + 3, 3 * i);

        // get E matrix representation of body quaternion from m_psi_conv_quat_ang
        const Eigen::Matrix<double, 3, 4> two_E_body =
            m_psi_conv_quat_ang.block<3, 4>(6 * body_num + 3, 7 * body_num + 3);

        // get change of variable matrix element D_{\alpha}
        const Eigen::Matrix<double, 7, 3> n_D_alpha = -m_conv_body_2_part_dof.block<7, 3>(7 * body_num, 3 * i);

        /* ANCHOR: tensor contractions to produce result */

        // convert Eigen::Matrix --> Eigen::Tensor
        const Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3>> tens_two_r_cross_mat = TensorCast(two_r_cross_mat);
        const Eigen::TensorFixedSize<double, Eigen::Sizes<3, 4>> tens_two_E_body      = TensorCast(two_E_body);
        const Eigen::TensorFixedSize<double, Eigen::Sizes<7, 3>> tens_n_D_alpha       = TensorCast(n_D_alpha);

        // 1) Compute result_{i m k} = levi_cevita{l i m} n_D_alpha_{k l}
        Eigen::array<Eigen::IndexPair<int>, 1> contract_lim_kl = {Eigen::IndexPair<int>(0, 1)}; // {i m k}
        Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3, 7>> d_rCrossMat_d_xi;
        d_rCrossMat_d_xi.device(device) = m_levi_cevita.contract(tens_n_D_alpha, contract_lim_kl);

        // 2) Compute result_{i j k} = d_rCrossMat_d_xi_{i m k} two_E_body{m j}
        Eigen::array<Eigen::IndexPair<int>, 1> contract_imk_mj = {Eigen::IndexPair<int>(1, 0)}; // {i k j}, must shuffle
        Eigen::array<int, 3>                   swap_last_two_indices({0, 2, 1});

        Eigen::TensorFixedSize<double, Eigen::Sizes<3, 7, 4>> d_rCrossMat_d_xi_times_two_E_body_preshuffle;
        d_rCrossMat_d_xi_times_two_E_body_preshuffle.device(device) =
            d_rCrossMat_d_xi.contract(tens_two_E_body, contract_imk_mj);

        Eigen::TensorFixedSize<double, Eigen::Sizes<3, 4, 7>> d_rCrossMat_d_xi_times_two_E_body;
        d_rCrossMat_d_xi_times_two_E_body.device(device) =
            d_rCrossMat_d_xi_times_two_E_body_preshuffle.shuffle(swap_last_two_indices);

        // 3) Compute result_{i j k} = two_r_cross_mat_{i m} Kappa_tilde_{m j k} (same indices as (2))
        Eigen::array<Eigen::IndexPair<int>, 1> contract_im_mjk = {Eigen::IndexPair<int>(1, 0)}; // {i j k}
        Eigen::TensorFixedSize<double, Eigen::Sizes<3, 4, 7>> rCrossMat_times_d_E_body_d_xi;
        rCrossMat_times_d_E_body_d_xi.device(device) = tens_two_r_cross_mat.contract(m_kappa_tilde, contract_im_mjk);

        // Compute and output matrix element
        Eigen::TensorFixedSize<double, Eigen::Sizes<3, 4, 7>> angular_gradient;
        angular_gradient.device(device) = -d_rCrossMat_d_xi_times_two_E_body - rCrossMat_times_d_E_body_d_xi;

        // output matrix indices
        Eigen::array<Eigen::Index, 3> offsets          = {3 * i, 7 * body_num + 3, 7 * body_num};
        Eigen::array<Eigen::Index, 3> extents          = {3, 4, 7};
        m_rbm_conn_T_quat_grad.slice(offsets, extents) = angular_gradient;
    }
    assert(body_num + 1 == m_num_bodies && "Not all bodies were indexed correctly");
}
