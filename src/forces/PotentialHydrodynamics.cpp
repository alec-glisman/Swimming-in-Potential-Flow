//
// Created by Alec Glisman on 07/31/21
//

#include <PotentialHydrodynamics.hpp>

PotentialHydrodynamics::PotentialHydrodynamics(std::shared_ptr<SystemData> sys)
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
    m_7N = 7 * m_system->numParticles();
    spdlog::get(m_logName)->info("Length of 7N tensor quantities: {0}", m_7N);
    m_7M = 7 * m_system->numBodies();
    spdlog::get(m_logName)->info("Length of 7M tensor quantities: {0}", m_7M);

    // set identity matrices
    m_I7N_linear  = Eigen::MatrixXd::Zero(m_7N, m_7N);
    m_I7N_angular = Eigen::MatrixXd::Zero(m_7N, m_7N);

    for (int particle_id = 0; particle_id < m_system->numParticles(); particle_id++)
    {
        const int particle_id_7{7 * particle_id};

        m_I7N_linear.block<3, 3>(particle_id_7, particle_id_7).noalias() = Eigen::MatrixXd::Identity(3, 3);

        /// @todo: (1, 1) entry is the added element to the mass matrix. Verify that changing this value does not affect
        /// the end result (Issue #6). Currently defaults to the same value as the intrinsic moment of inertia.
        m_I7N_angular.block<4, 4>(particle_id_7 + 3, particle_id_7 + 3).noalias() = Eigen::MatrixXd::Identity(4, 4);
    }

    m_c1_2_I7N_linear = m_c1_2 * m_I7N_linear;

    // Initialize mass matrices
    spdlog::get(m_logName)->info("Initializing mass matrices");

    m_M_intrinsic = (m_system->particleDensity() * m_unit_sphere_volume) * m_I7N_linear;

    m_J_intrinsic = (m_scalar_moment_inertia * m_system->particleDensity() * m_unit_sphere_volume) *
                    m_I7N_angular; // intrinsic moment of inertia for spheres

    m_M_added = Eigen::MatrixXd::Zero(m_7N, m_7N);

    m_M_total.noalias() = m_M_added;
    m_M_total.noalias() += m_M_intrinsic;
    m_M_total.noalias() += m_J_intrinsic;

    spdlog::get(m_logName)->info("Initializing mass tensors");
    m_grad_M_added = Eigen::Tensor<double, 3>(m_7N, m_7N, m_7N);
    m_grad_M_added.setZero();
    m_grad_M_added_body_coords = Eigen::Tensor<double, 3>(m_7N, m_7N, m_7M);
    m_grad_M_added_body_coords.setZero();
    m_tens_M_total = Eigen::Tensor<double, 2>(m_7N, m_7N);
    m_tens_M_total.setZero();

    spdlog::get(m_logName)->info("Initializing tensors used in hydrodynamic force calculations.");
    m_N1 = Eigen::Tensor<double, 3>(m_7N, m_7N, m_7M);
    m_N1.setZero();
    m_N2 = Eigen::Tensor<double, 3>(m_7M, m_7N, m_7M);
    m_N2.setZero();
    m_N3 = Eigen::Tensor<double, 3>(m_7M, m_7M, m_7M);
    m_N3.setZero();

    m_M2 = Eigen::Tensor<double, 2>(m_7M, m_7N);
    m_M2.setZero();
    m_M3 = Eigen::Tensor<double, 2>(m_7M, m_7M);
    m_M3.setZero();
    m_mat_M3 = Eigen::MatrixXd::Zero(m_7M, m_7M);

    // N^{(2)}
    N2_term1_preshuffle = Eigen::Tensor<double, 3>(m_7M, m_7M, m_7N);
    N2_term1_preshuffle.setZero();
    N2_term1 = Eigen::Tensor<double, 3>(m_7M, m_7N, m_7M);
    N2_term1.setZero();
    N2_term2 = Eigen::Tensor<double, 3>(m_7M, m_7N, m_7M);
    N2_term2.setZero();

    // N^{(3)}
    N3_term1_preshuffle = Eigen::Tensor<double, 3>(m_7M, m_7M, m_7M);
    N3_term3.setZero();
    N3_term2_preshuffle = Eigen::Tensor<double, 3>(m_7M, m_7M, m_7M);
    N3_term2_preshuffle.setZero();
    N3_term1 = Eigen::Tensor<double, 3>(m_7M, m_7M, m_7M);
    N3_term1.setZero();
    N3_term2 = Eigen::Tensor<double, 3>(m_7M, m_7M, m_7M);
    N3_term2.setZero();
    N3_term3 = Eigen::Tensor<double, 3>(m_7M, m_7M, m_7M);
    N3_term3.setZero();

    // Assign particle pair information
    spdlog::get(m_logName)->info("Initializing particle pair information vectors");
    m_alphaVec = Eigen::VectorXi::Zero(m_num_pair_inter);
    m_betaVec  = Eigen::VectorXi::Zero(m_num_pair_inter);
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

        m_alphaVec(i) = alpha;
        m_betaVec(i)  = beta;
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

PotentialHydrodynamics::~PotentialHydrodynamics()
{
    spdlog::get(m_logName)->info("Destructing potential hydrodynamics");
    spdlog::get(m_logName)->flush();
    spdlog::drop(m_logName);
}

void
PotentialHydrodynamics::update(Eigen::ThreadPoolDevice& device)
{
    calcParticleDistances();

    calcAddedMass();
    calcAddedMassGrad(device);

    calcTotalMass();

    calcBodyTensors(device);

    calcHydroForces(device);
}

void
PotentialHydrodynamics::calcParticleDistances()
{
    /* NOTE: Fill Mass matrix elements one (3 x 3) block at a time (matrix elements between
     * particles \alpha and \beta) */
    for (int i = 0; i < m_num_pair_inter; i++)
    {
        m_r_ab.col(i).noalias() = m_system->positionsParticles().segment<3>(3 * m_alphaVec(i));
        m_r_ab.col(i).noalias() -= m_system->positionsParticles().segment<3>(3 * m_betaVec(i));

        m_r_mag_ab(i) = m_r_ab.col(i).norm(); //(1); |r| between 2 particles

#if !defined(NDEBUG)
        spdlog::get(m_logName)->critical("Checking distance between particle pair [{0}, {1}]: {2:.3f}", m_alphaVec(i),
                                         m_betaVec(i), m_r_mag_ab(i));
        spdlog::get(m_logName)->flush();
#endif
    }
}

void
PotentialHydrodynamics::calcAddedMass()
{
    // set matrices to zero
    m_M_added.setZero();

    /* Fill off-diagonal elements (without units ) */
    /* NOTE: Fill Mass matrix elements one (3 x 3) block at a time (matrix elements between
     * particles \alpha and \beta) */
    for (int k = 0; k < m_num_pair_inter; k++)
    {
        // Convert (\alpha, \beta) --> (i, j) by factor of 3
        const int i_part{7 * m_alphaVec(k)};
        const int j_part{7 * m_betaVec(k)};

        // Full distance between particles \alpha and \beta
        const Eigen::Vector3d r_ij = m_r_ab.col(k);    // (1)
        const double          r_mag_ij{m_r_mag_ab(k)}; //(1); |r| between 2 particles

        // M^{(1)} Matrix Element Constants:
        const double M1_c1{-m_c3_2 / std::pow(r_mag_ij, 5)}; // (1)
        const double M1_c2{m_c1_2 / std::pow(r_mag_ij, 3)};  // (1)

        // Full matrix elements for M^{(1)}_{ij} (NOTE: missing factor of 1/2)
        Eigen::Matrix3d Mij = r_ij * r_ij.transpose(); //(1); Outer product of \bm{r} \bm{r}
        Mij *= M1_c1;
        Mij.noalias() += M1_c2 * m_system->i3();

        // Output added mass element (symmetry of mass matrix)
        m_M_added.block<3, 3>(i_part, j_part).noalias() = Mij;
        m_M_added.block<3, 3>(j_part, i_part).noalias() = Mij;
    }

    /* Construct full added mass matrix
     * M = 1/2 I + M^{(1)}
     */
    m_M_added.noalias() += m_c1_2_I7N_linear;                       // Add diagonal elements
    m_M_added *= (m_system->fluidDensity() * m_unit_sphere_volume); // mass units
}

void
PotentialHydrodynamics::calcAddedMassGrad(Eigen::ThreadPoolDevice& device)
{
    // `Eigen::Tensor` contraction indices
    const Eigen::array<Eigen::IndexPair<int>, 1>  contract_ijl_kl = {Eigen::IndexPair<int>(2, 1)};
    const Eigen::array<Eigen::IndexPair<long>, 0> outer_product   = {};

    // `Eigen::Tensor` permutation indices
    const Eigen::array<int, 3> permute_ijk_kij({2, 0, 1}); // {i, j, k} --> {k, i, j}
    const Eigen::array<int, 3> permute_ijk_jki({1, 2, 1}); // {i, j, k} --> {j, k, i}

    // set matrices to zero
    m_grad_M_added.setZero();

    /* NOTE: Fill Mass matrix elements one (3 x 3 x 3) block at a time (matrix elements between
     * particles \alpha and \beta) */

    for (int k = 0; k < m_num_pair_inter; k++)
    {
        // Convert (\alpha, \beta) --> (i, j) by factor of 3
        int i_part{7 * m_alphaVec(k)};
        int j_part{7 * m_betaVec(k)};

        // Full distance between particles \alpha and \beta
        const Eigen::Vector3d                                 r_ij   = m_r_ab.col(k); // (1)
        const Eigen::TensorFixedSize<double, Eigen::Sizes<3>> tens_r = TensorCast(r_ij, 3);

        const double          r_mag_ij{m_r_mag_ab(k)};            //(1); |r| between 2 particles
        const Eigen::Matrix3d r_dyad_r = r_ij * r_ij.transpose(); // (1); Outer product of \bm{r} \bm{r}
        const Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3>> tens_rr = TensorCast(r_dyad_r, 3, 3);

        // Constants to use in Calculation
        double gradM1_c1{-(m_system->fluidDensity() * m_unit_sphere_volume) * m_c3_2 *
                         std::pow(r_mag_ij, -5)}; // mass units
        double gradM1_c2{(m_system->fluidDensity() * m_unit_sphere_volume) * m_c15_2 *
                         std::pow(r_mag_ij, -7)}; // mass units

        // outer products (I_{i j} r_{k}) and permutations
        Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3, 3>> delta_ij_r_k;
        delta_ij_r_k.device(device) = gradM1_c1 * m_system->tensI3().contract(tens_r, outer_product);
        // shuffle all dimensions to the right by 1: (i, j, k) --> (k, i, j), (2, 0, 1)
        Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3, 3>> delta_jk_r_i;
        delta_jk_r_i.device(device) = delta_ij_r_k.shuffle(permute_ijk_kij);
        // shuffle all dimensions to the left by 1: (i, j, k) --> (k, i, j), (1, 2, 0)
        Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3, 3>> delta_ki_r_j;
        delta_ki_r_j.device(device) = delta_ij_r_k.shuffle(permute_ijk_jki);

        // full matrix element for M_{i j, i}
        Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3, 3>> Mij_i;
        Mij_i.device(device) = delta_ij_r_k + delta_jk_r_i;
        Mij_i.device(device) += delta_ki_r_j;
        Mij_i.device(device) += gradM1_c2 * tens_rr.contract(tens_r, outer_product);

        // indices to start at
        const Eigen::array<Eigen::Index, 3> offsets_ij_i = {i_part, j_part, i_part};
        const Eigen::array<Eigen::Index, 3> offsets_ij_j = {i_part, j_part, j_part};
        const Eigen::array<Eigen::Index, 3> offsets_ji_i = {j_part, i_part, i_part};
        const Eigen::array<Eigen::Index, 3> offsets_ji_j = {j_part, i_part, j_part};
        // length of data to access
        const Eigen::array<Eigen::Index, 3> extents = {3, 3, 3};

        // M_{ij, i}: Matrix Element (Anti-Symmetric upon exchange of derivative, Symmetric upon
        // exchange of first two indices)
        m_grad_M_added.slice(offsets_ij_i, extents).device(device) = Mij_i;

        // M_{ij, j} = - M_{ij, i}
        m_grad_M_added.slice(offsets_ij_j, extents).device(device) = -Mij_i;

        // M_{ji, j} = - M_{ij, i}
        m_grad_M_added.slice(offsets_ji_j, extents).device(device) = -Mij_i;

        // M_{ji, i} = M_{ij, i}
        m_grad_M_added.slice(offsets_ji_i, extents).device(device) = Mij_i;
    }

    m_grad_M_added_body_coords.device(device) = m_grad_M_added.contract(m_system->tensChi(), contract_ijl_kl);
}

void
PotentialHydrodynamics::calcTotalMass()
{
    m_M_total.noalias() = m_M_intrinsic;
    m_M_total.noalias() += m_J_intrinsic;
    m_M_total.noalias() += m_M_added;

    m_tens_M_total = TensorCast(m_M_total, m_7N, m_7N);
}

void
PotentialHydrodynamics::calcBodyTensors(Eigen::ThreadPoolDevice& device)
{
    /* ANCHOR: Indices */
    // `Eigen::Tensor` contraction indices
    const Eigen::array<Eigen::IndexPair<int>, 1> contract_il_lj = {Eigen::IndexPair<int>(1, 0)}; // = A B
    const Eigen::array<Eigen::IndexPair<int>, 1> contract_il_jl = {Eigen::IndexPair<int>(1, 1)}; // = A B^T

    // `Eigen::Tensor` permutation indices
    const Eigen::array<int, 3> permute_ikj_ijk({0, 2, 1}); // (i, j, k) --> (i, j, k)

    /* ANCHOR: Compute linear combinations of total mass matrix and zeta */
    m_M2.device(device) = m_system->tensZeta().contract(m_tens_M_total, contract_il_lj);

    m_M3.device(device) = m_M2.contract(m_system->tensZeta(), contract_il_jl);
    m_mat_M3            = MatrixCast(m_M3, m_7M, m_7M, device);

    /* ANCHOR: Compute linear combinations of GRADIENTS of total mass matrix and zeta */
    // N^{(1)}
    m_N1.device(device) = m_grad_M_added_body_coords;

    // N^{(2)}
    N2_term1_preshuffle.device(device) = m_system->tensGradZeta().contract(m_tens_M_total, contract_il_lj);

    m_N2.device(device) = N2_term1_preshuffle.shuffle(permute_ikj_ijk);
    m_N2.device(device) += m_system->tensZeta().contract(m_N1, contract_il_lj);

    // N^{(3)}
    N3_term1_preshuffle.device(device) = N2_term1.contract(m_system->tensZeta(), contract_il_jl);
    N3_term1_preshuffle.device(device) += N2_term2.contract(m_system->tensZeta(), contract_il_jl);

    m_N3.device(device) = N3_term1_preshuffle.shuffle(permute_ikj_ijk);
    m_N3.device(device) += m_M2.contract(m_system->tensGradZeta(), contract_il_jl);
}

void
PotentialHydrodynamics::calcHydroForces(Eigen::ThreadPoolDevice& device)
{
    // `Eigen::Tensor` contraction indices
    // 0th order contraction
    const Eigen::array<Eigen::IndexPair<int>, 0> outer_product = {};
    // 1st order contraction
    const Eigen::array<Eigen::IndexPair<int>, 1> contract_ij_j = {Eigen::IndexPair<int>(1, 0)};
    // 2nd order contraction
    const Eigen::array<Eigen::IndexPair<int>, 2> contract_jki_jk = {Eigen::IndexPair<int>(0, 0),
                                                                    Eigen::IndexPair<int>(1, 1)};
    const Eigen::array<Eigen::IndexPair<int>, 2> contract_ijk_jk = {Eigen::IndexPair<int>(1, 0),
                                                                    Eigen::IndexPair<int>(2, 1)};

    // `Eigen::Tensor` permutation indices
    const Eigen::array<int, 3> permute_ij_ji({1, 0});
    const Eigen::array<int, 3> permute_ijk_jik({1, 0, 2});

    // get kinematic tensors from SystemData class
    Eigen::Tensor<double, 1> xi_dot  = TensorCast(m_system->velocitiesBodies(), m_7M);    // (7M x 1)
    Eigen::Tensor<double, 1> xi_ddot = TensorCast(m_system->accelerationsBodies(), m_7M); // (7M x 1)

    Eigen::Tensor<double, 1> V     = TensorCast(m_system->velocitiesParticlesArticulation(), m_7N);    // (7N x 1)
    Eigen::Tensor<double, 1> V_dot = TensorCast(m_system->accelerationsParticlesArticulation(), m_7N); // (7N x 1)

    // calculate 2nd order kinematic tensors
    Eigen::Tensor<double, 2> V_V = Eigen::Tensor<double, 2>(m_7N, m_7N); // (7N x 7N)
    V_V.device(device)           = V.contract(V, outer_product);

    Eigen::Tensor<double, 2> xi_dot_xi_dot = Eigen::Tensor<double, 2>(m_7M, m_7M); // (7M x 7M)
    xi_dot_xi_dot.device(device)           = xi_dot.contract(xi_dot, outer_product);

    Eigen::Tensor<double, 2> xi_dot_V = Eigen::Tensor<double, 2>(m_7M, m_7N); // (7M x 7N)
    xi_dot_V.device(device)           = xi_dot.contract(V, outer_product);

    Eigen::Tensor<double, 2> V_xi_dot = Eigen::Tensor<double, 2>(m_7N, m_7M); // (7N x 7M)
    V_xi_dot.device(device)           = xi_dot_V.shuffle(permute_ij_ji);

    // hydrodynamic forces arising from locater point motion (3 terms; contains locater inertia term)
    Eigen::Tensor<double, 1> inertia_loc = Eigen::Tensor<double, 1>(m_7M); // (7M x 1)
    inertia_loc.device(device)           = -m_M3.contract(xi_ddot, contract_ij_j);

    Eigen::Tensor<double, 1> F_loc = Eigen::Tensor<double, 1>(m_7M); // (7M x 1)
    F_loc.device(device)           = inertia_loc;
    F_loc.device(device) += 0.50 * m_N3.contract(xi_dot_xi_dot, contract_jki_jk);
    F_loc.device(device) -= m_N3.contract(xi_dot_xi_dot, contract_ijk_jk);

    // hydrodynamic forces arising from coupling of locater and internal D.o.F. motion (2 terms)
    Eigen::Tensor<double, 1> F_loc_int = Eigen::Tensor<double, 1>(m_7M); // (7M x 1)
    F_loc_int.device(device)           = m_N2.contract(xi_dot_V, contract_jki_jk);
    F_loc_int.device(device) -= m_N2.contract(V_xi_dot, contract_ijk_jk);

    // hydrodynamic forces arising from internal D.o.F. motion (2 terms; contains internal inertia term)
    Eigen::Tensor<double, 1> F_int = Eigen::Tensor<double, 1>(m_7M); // (7M x 1)
    F_int.device(device)           = 0.50 * m_N1.contract(V_V, contract_jki_jk);
    F_int.device(device) -= m_M2.contract(V_dot, contract_ij_j);

    // compute complete potential flow hydrodynamic force
    Eigen::Tensor<double, 1> F_hydro = Eigen::Tensor<double, 1>(m_7M); // (7M x 1)
    F_hydro.device(device)           = F_loc + F_loc_int;
    F_hydro.device(device) += F_int;
    m_F_hydro.noalias() = MatrixCast(F_hydro, m_7M, 1, device);

    // Compute non-inertial part of force (include internal D.o.F. inertia)
    Eigen::Tensor<double, 1> F_hydro_no_inertia = Eigen::Tensor<double, 1>(m_7M); // (7M x 1)
    F_hydro_no_inertia.device(device)           = F_hydro - inertia_loc;
    m_F_hydroNoInertia.noalias()                = MatrixCast(F_hydro_no_inertia, m_7M, 1, device);
}
