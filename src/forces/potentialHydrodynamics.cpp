//
// Created by Alec Glisman on 07/31/21
//

#include <potentialHydrodynamics.hpp>

// REVIEW[epic=Debug]: Uncomment line below to prevent all runtime checks from executing in debug
// #define NO_HYDRO_CHECK

potentialHydrodynamics::potentialHydrodynamics(std::shared_ptr<systemData> sys)
{
    // save classes
    m_system = sys;

    // Initialize logger
    m_logFile   = m_system->outputDir() + "/logs/" + m_logName + "-log.txt";
    auto logger = spdlog::basic_logger_mt(m_logName, m_logFile);
    spdlog::get(m_logName)->info("Initializing potential hydrodynamics");

    // Variables for for-loop
    m_num_pair_inter = m_system->numParticles() * (m_system->numParticles() - 1) / 2; // Number of interactions to count
    spdlog::get(m_logName)->info("Setting number of interactions to count: {0}", m_num_pair_inter);

    // tensor variables
    m_3N = 3 * m_system->numParticles();
    spdlog::get(m_logName)->info("Length of 3N tensor quantities: {0}", m_3N);
    m_7M = 7 * m_system->numBodies();
    spdlog::get(m_logName)->info("Length of 7M tensor quantities: {0}", m_7M);

    // set identity matrices
    m_I3N      = Eigen::MatrixXd::Identity(m_3N, m_3N);
    m_c1_2_I3N = m_c1_2 * m_I3N;

    // Initialize mass matrices
    spdlog::get(m_logName)->info("Initializing mass matrices");
    m_M_intrinsic = (m_system->particleDensity() * m_unitSphereVol) * m_I3N;
    m_M_added     = Eigen::MatrixXd::Zero(m_3N, m_3N);

    m_M_total.noalias() = m_M_added;
    m_M_total.noalias() += m_M_intrinsic;

    spdlog::get(m_logName)->info("Initializing mass tensors");
    m_grad_M_added = Eigen::Tensor<double, 3>(m_3N, m_3N, m_3N);
    m_grad_M_added.setZero();
    m_grad_M_added_body_coords = Eigen::Tensor<double, 3>(m_3N, m_3N, m_7M);
    m_grad_M_added_body_coords.setZero();
    m_tens_M_total = Eigen::Tensor<double, 2>(m_3N, m_3N);
    m_tens_M_total.setZero();

    spdlog::get(m_logName)->info("Initializing tensors used in hydrodynamic force calculations.");
    m_N1 = Eigen::Tensor<double, 3>(m_3N, m_3N, m_7M);
    m_N1.setZero();
    m_N2 = Eigen::Tensor<double, 3>(m_7M, m_3N, m_7M);
    m_N2.setZero();
    m_N3 = Eigen::Tensor<double, 3>(7 * m_system->numBodies(), 7 * m_system->numBodies(), 7 * m_system->numBodies());
    m_N3.setZero();

    m_M_tilde_tilde = Eigen::Tensor<double, 2>(7 * m_system->numBodies(), m_3N);
    m_M_tilde_tilde.setZero();
    m_M_tilde = Eigen::Tensor<double, 2>(m_7M, m_7M);
    m_M_tilde.setZero();
    m_mat_M_tilde = Eigen::MatrixXd::Zero(m_7M, m_7M);

    // Assign particle pair information
    spdlog::get(m_logName)->info("Initializing particle pair information vectors");
    m_alphaVec = Eigen::VectorXd::Zero(m_num_pair_inter);
    m_betaVec  = Eigen::VectorXd::Zero(m_num_pair_inter);
    m_r_mag_ab = Eigen::VectorXd::Zero(m_num_pair_inter);
    m_r_ab     = Eigen::MatrixXd::Zero(3, m_num_pair_inter);

    /* Fill the particle index vectors
     * Calculate ahead of time to save time during runtime
     */
    spdlog::get(m_logName)->info("Filling particle pair information tensors");
    for (int i = 0; i < m_num_pair_inter; i++)
    {
        /* alpha and beta convention to convert from linear coordinate to ordered pair
         * @Source:
         * https://stackoverflow.com/questions/33810187/openmp-and-c-private-variables/33836073#33836073
         */
        int alpha = i / m_system->numParticles();
        int beta  = i % m_system->numParticles();

        if (beta <= alpha)
        {
            alpha = m_system->numParticles() - alpha - 2;
            beta  = m_system->numParticles() - beta - 1;
        }

        assert(alpha >= 0 && "Calculated alpha must be non-negative");
        assert(beta >= 0 && "Calculated beta must be non-negative");

        m_alphaVec[i] = alpha;
        m_betaVec[i]  = beta;
    }

    // Compute all relevant quantities
    spdlog::get(m_logName)->critical("Setting up temporary single-thread eigen device");
    Eigen::ThreadPool       thread_pool = Eigen::ThreadPool(1);
    Eigen::ThreadPoolDevice single_core_device(&thread_pool, 1);

    spdlog::get(m_logName)->info("Calling update()");
    update(single_core_device);

    spdlog::get(m_logName)->info("Constructor complete");
    spdlog::get(m_logName)->flush();
}

potentialHydrodynamics::~potentialHydrodynamics()
{
    spdlog::get(m_logName)->info("Destructing potential hydrodynamics");
    spdlog::get(m_logName)->flush();
}

void
potentialHydrodynamics::update(Eigen::ThreadPoolDevice& device)
{
    calcParticleDistances();

    calcAddedMass();
    calcAddedMassGrad(device);

    calcTotalMass();

    calcBodyTensors(device);

    calcHydroForces(device);
}

void
potentialHydrodynamics::calcParticleDistances()
{
    /* NOTE: Fill Mass matrix elements one (3 x 3) block at a time (matrix elements between
     * particles \alpha and \beta) */
    for (int i = 0; i < m_num_pair_inter; i++)
    {
        m_r_ab.col(i).noalias() = m_system->positionsParticles()(Eigen::seqN(3 * m_alphaVec[i], 3));
        m_r_ab.col(i).noalias() -= m_system->positionsParticles()(Eigen::seqN(3 * m_betaVec[i], 3));

        m_r_mag_ab[i] = m_r_ab.col(i).norm(); //[1]; |r| between 2 particles

#if !defined(NDEBUG) && !defined(NO_HYDRO_CHECKS)
        spdlog::get(m_logName)->critical("Checking distance between pairs {0} & {1}", m_alphaVec[i], m_betaVec[i]);
        spdlog::get(m_logName)->critical("Interparticle distance is", m_r_mag_ab[i]);
        spdlog::get(m_logName)->flush();
#endif
    }
}

void
potentialHydrodynamics::calcAddedMass()
{
    // set matrices to zero
    m_M_added.setZero(m_3N, m_3N);

    /* Fill off-diagonal elements (without units ) */
    /* NOTE: Fill Mass matrix elements one (3 x 3) block at a time (matrix elements between
     * particles \alpha and \beta) */
    for (int k = 0; k < m_num_pair_inter; k++)
    {
        // Convert (\alpha, \beta) --> (i, j) by factor of 3
        int i_part = 3 * m_alphaVec[k];
        int j_part = 3 * m_betaVec[k];

        // Full distance between particles \alpha and \beta
        Eigen::Vector3d r_ij     = m_r_ab.col(k); // [1]
        double          r_mag_ij = m_r_mag_ab[k]; //[1]; |r| between 2 particles

        // M^{(1)} Matrix Element Constants:
        double M1_c1 = -m_c3_2 / std::pow(r_mag_ij, 5); // [1]
        double M1_c2 = m_c1_2 / std::pow(r_mag_ij, 3);  // [1]

        // Full matrix elements for M^{(1)}_{ij} (NOTE: missing factor of 1/2)
        Eigen::Matrix3d Mij = r_ij * r_ij.transpose(); //[1]; Outer product of \bm{r} \bm{r}
        Mij *= M1_c1;
        Mij.noalias() += M1_c2 * m_system->i3();

        // Output added mass element (symmetry of mass matrix)
        m_M_added.block<3, 3>(i_part, j_part).noalias() = Mij;
        m_M_added.block<3, 3>(j_part, i_part).noalias() = Mij;
    }

    /* Construct full added mass matrix
     * M = 1/2 I + M^{(1)}
     */
    m_M_added.noalias() += m_c1_2_I3N;                         // Add diagonal elements
    m_M_added *= (m_system->fluidDensity() * m_unitSphereVol); // mass units
}

void
potentialHydrodynamics::calcAddedMassGrad(Eigen::ThreadPoolDevice& device)
{
    // eigen tensor contraction variables
    Eigen::array<Eigen::IndexPair<int>, 1>  contract_ijl_lk  = {Eigen::IndexPair<int>(2, 0)};
    Eigen::array<Eigen::IndexPair<long>, 0> empty_index_list = {};
    Eigen::array<int, 3>                    shuffle_one_right({2, 0, 1});
    Eigen::array<int, 3>                    shuffle_one_left({1, 2, 1});

    // set matrices to zero
    m_grad_M_added.setZero();

    /* NOTE: Fill Mass matrix elements one (3 x 3 x 3) block at a time (matrix elements between
     * particles \alpha and \beta) */

    for (int k = 0; k < m_num_pair_inter; k++)
    {
        // Convert (\alpha, \beta) --> (i, j) by factor of 3
        int i_part = 3 * m_alphaVec[k];
        int j_part = 3 * m_betaVec[k];

        // Full distance between particles \alpha and \beta
        Eigen::Vector3d r_ij     = m_r_ab.col(k);           // [1]
        double          r_mag_ij = m_r_mag_ab[k];           //[1]; |r| between 2 particles
        Eigen::Matrix3d r_dyad_r = r_ij * r_ij.transpose(); // [1]; Outer product of \bm{r} \bm{r}
        Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3>> tens_rr = TensorCast(r_dyad_r);

        // Constants to use in Calculation
        double gradM1_c1 =
            -(m_system->fluidDensity() * m_unitSphereVol) * m_c3_2 * std::pow(r_mag_ij, -5); // mass units
        double gradM1_c2 =
            (m_system->fluidDensity() * m_unitSphereVol) * m_c15_2 * std::pow(r_mag_ij, -7); // mass units

        // outer products (I_{i j} r_{k}) and permutations
        Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3, 3>> delta_ij_r_k =
            gradM1_c1 * m_system->tensI3().contract(TensorCast(r_ij), empty_index_list);
        // shuffle all dimensions to the right by 1: (i, j, k) --> (k, i, j), (2, 0, 1)
        Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3, 3>> delta_jk_r_i = delta_ij_r_k.shuffle(shuffle_one_right);
        // shuffle all dimensions to the left by 1: (i, j, k) --> (k, i, j), (1, 2, 0)
        Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3, 3>> delta_ki_r_j = delta_ij_r_k.shuffle(shuffle_one_left);

        // full matrix element for M_{i j, i}
        Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3, 3>> Mij_i = delta_ij_r_k + delta_jk_r_i;
        Mij_i += delta_ki_r_j;
        Mij_i += gradM1_c2 * tens_rr.contract(TensorCast(r_ij), empty_index_list);

        // indices to start at
        Eigen::array<Eigen::Index, 3> offsets_ij_i = {i_part, j_part, i_part};
        Eigen::array<Eigen::Index, 3> offsets_ij_j = {i_part, j_part, j_part};
        Eigen::array<Eigen::Index, 3> offsets_ji_i = {j_part, i_part, i_part};
        Eigen::array<Eigen::Index, 3> offsets_ji_j = {j_part, i_part, j_part};
        // length of data to access
        Eigen::array<Eigen::Index, 3> extents = {3, 3, 3};

        // M_{ij, i}: Matrix Element (Anti-Symmetric upon exchange of derivative, Symmetric upon
        // exchange of first two indices)
        m_grad_M_added.slice(offsets_ij_i, extents) = Mij_i;

        // M_{ij, j} = - M_{ij, i}
        m_grad_M_added.slice(offsets_ij_j, extents) = -Mij_i;

        // M_{ji, j} = - M_{ij, i}
        m_grad_M_added.slice(offsets_ji_j, extents) = -Mij_i;

        // M_{ji, i} = M_{ij, i}
        m_grad_M_added.slice(offsets_ji_i, extents) = Mij_i;
    }

    m_grad_M_added_body_coords.device(device) =
        m_grad_M_added.contract(m_system->tensConvBody2PartDof(), contract_ijl_lk);
}

void
potentialHydrodynamics::calcTotalMass()
{
    m_M_total.noalias() = m_M_intrinsic;
    m_M_total.noalias() += m_M_added;

    m_tens_M_total = TensorCast(m_M_total);
}

void
potentialHydrodynamics::calcBodyTensors(Eigen::ThreadPoolDevice& device)
{
    // contraction indices
    Eigen::array<Eigen::IndexPair<int>, 1> contract_li_lj = {Eigen::IndexPair<int>(0, 0)}; // = A^T B
    Eigen::array<Eigen::IndexPair<int>, 1> contract_il_lj = {Eigen::IndexPair<int>(1, 0)}; // = A B

    // tensor index shuffling
    Eigen::array<int, 3> flip_first_two_indices({1, 0, 2}); // (i, j, k) --> (j, i, k)

    /* ANCHOR: Compute linear combinations of total mass matrix and rbm_conn_t_quat */
    m_M_tilde_tilde.device(device) = m_system->tensRbmConnTQuat().contract(m_tens_M_total, contract_li_lj);

    m_M_tilde.device(device) = m_M_tilde_tilde.contract(m_system->tensRbmConnTQuat(), contract_il_lj);
    m_mat_M_tilde            = MatrixCast(m_M_tilde_tilde, m_7M, m_7M);

    /* ANCHOR: Compute linear combinations of GRADIENTS of total mass matrix and rbm_conn_t_quat */
    // N^{(1)}
    m_N1 = m_grad_M_added_body_coords;

    // Get transpose (of first two indices) of \nabla C
    Eigen::Tensor<double, 3> rbm_conn_quat   = Eigen::Tensor<double, 3>(m_3N, m_3N, m_7M);
    Eigen::Tensor<double, 3> rbm_conn_quat_T = m_system->tensRbmConnTQuat();
    rbm_conn_quat.device(device)             = rbm_conn_quat_T.shuffle(flip_first_two_indices);

    // N^{(2)}
    Eigen::Tensor<double, 3> N2_term1 = Eigen::Tensor<double, 3>(m_7M, m_3N, m_7M);
    Eigen::Tensor<double, 3> N2_term2 = Eigen::Tensor<double, 3>(m_7M, m_3N, m_7M);

    N2_term1.device(device) = rbm_conn_quat.contract(m_tens_M_total, contract_il_lj);
    N2_term2.device(device) = m_system->tensRbmConnTQuat().contract(m_N1, contract_li_lj);
    m_N2.device(device)     = N2_term1 + N2_term2;

    // N^{(3)}
    Eigen::Tensor<double, 3> N3_term1 = Eigen::Tensor<double, 3>(m_7M, m_7M, m_7M);
    Eigen::Tensor<double, 3> N3_term2 = Eigen::Tensor<double, 3>(m_7M, m_7M, m_7M);
    Eigen::Tensor<double, 3> N3_term3 = Eigen::Tensor<double, 3>(m_7M, m_7M, m_7M);

    N3_term1.device(device) = N2_term1.contract(m_system->tensRbmConnTQuat(), contract_il_lj);
    N3_term2.device(device) = N2_term2.contract(m_system->tensRbmConnTQuat(), contract_il_lj);
    N3_term3.device(device) = m_M_tilde_tilde.contract(m_system->tensRbmConnTQuatGrad(), contract_il_lj);

    m_N3.device(device) = N3_term1 + N3_term2;
    m_N3.device(device) += N3_term3;
}

void
potentialHydrodynamics::calcHydroForces(Eigen::ThreadPoolDevice& device)
{
    // tensor contraction indices
    Eigen::array<Eigen::IndexPair<int>, 0> empty_index_list = {}; // outer product
    Eigen::array<Eigen::IndexPair<int>, 2> contract_jki_jk = {Eigen::IndexPair<int>(0, 0), Eigen::IndexPair<int>(1, 1)};
    Eigen::array<Eigen::IndexPair<int>, 2> contract_ijk_jk = {Eigen::IndexPair<int>(1, 0), Eigen::IndexPair<int>(2, 1)};
    Eigen::array<Eigen::IndexPair<int>, 1> contract_ij_j   = {Eigen::IndexPair<int>(1, 0)};

    // get kinematic tensors from systemData class
    Eigen::Tensor<double, 1> xi_dot  = TensorCast(m_system->velocitiesBodies());
    Eigen::Tensor<double, 1> xi_ddot = TensorCast(m_system->accelerationsBodies());

    Eigen::Tensor<double, 1> V     = TensorCast(m_system->velocitiesParticlesArticulation());
    Eigen::Tensor<double, 1> V_dot = TensorCast(m_system->accelerationsParticlesArticulation());

    // calculate 2nd order kinematic tensors
    Eigen::Tensor<double, 2> V_V = Eigen::Tensor<double, 2>(m_3N, m_3N);
    V_V.device(device)           = V.contract(V, empty_index_list);

    Eigen::Tensor<double, 2> xi_dot_xi_dot = Eigen::Tensor<double, 2>(m_7M, m_7M);
    xi_dot_xi_dot.device(device)           = xi_dot.contract(xi_dot, empty_index_list);

    Eigen::Tensor<double, 2> xi_dot_V = Eigen::Tensor<double, 2>(m_7M, m_3N);
    xi_dot_V.device(device)           = xi_dot.contract(V, empty_index_list);

    // hydrodynamic forces arising from locater point motion (3 terms; contains locater inertia term)
    Eigen::Tensor<double, 1> F_loc       = Eigen::Tensor<double, 1>(m_7M); // TODO
    Eigen::Tensor<double, 1> inertia_loc = Eigen::Tensor<double, 1>(m_7M); // TODO

    // hydrodynamic forces arising from coupling of locater and internal D.o.F. motion (2 terms)
    Eigen::Tensor<double, 1> F_loc_int = Eigen::Tensor<double, 1>(m_7M); // TODO

    // hydrodynamic forces arising from internal D.o.F. motion (2 terms; contains internal inertia term)
    Eigen::Tensor<double, 1> F_int = Eigen::Tensor<double, 1>(m_7M); // TODO

    // compute complete potential flow hydrodynamic force
    Eigen::Tensor<double, 1> F_hydro = Eigen::Tensor<double, 1>(m_7M);
    F_hydro.device(device)           = F_loc + F_loc_int;
    F_hydro.device(device) += F_int;
    m_F_hydro.noalias() = MatrixCast(F_hydro, m_7M, 1, device);

    // Compute non-inertial part of force (include internal D.o.F. inertia)
    Eigen::Tensor<double, 1> F_hydro_no_inertia = Eigen::Tensor<double, 1>(m_7M);
    F_hydro_no_inertia.device(device)           = F_hydro - inertia_loc;
    m_F_hydroNoInertia.noalias()                = MatrixCast(F_hydro_no_inertia, m_7M, 1, device);
}
