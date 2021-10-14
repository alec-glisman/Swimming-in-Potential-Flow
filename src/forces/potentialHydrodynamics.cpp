//
// Created by Alec Glisman on 07/31/21
//

#include <potentialHydrodynamics.hpp>

// REVIEW[epic=Debug]: Uncomment line below to prevent all runtime checks from executing in debug
#define NO_HYDRO_CHECK

potentialHydrodynamics::potentialHydrodynamics(std::shared_ptr<systemData> sys)
{
    // save classes
    m_system = sys;

    // Initialize logger
    m_logFile   = m_system->outputDir() + "/logs/" + m_logName + "-log.txt";
    auto logger = spdlog::basic_logger_mt(m_logName, m_logFile);
    spdlog::get(m_logName)->info("Initializing potential hydrodynamics");

    // Variables for for-loop
    m_num_inter = m_system->numParticles() * (m_system->numParticles() - 1) / 2; // Number of interactions to count
    spdlog::get(m_logName)->info("Setting number of interactions to count: {0}", m_num_inter);

    // tensor variables
    m_3N = 3 * m_system->numParticles(); // length of tensor quantities
    spdlog::get(m_logName)->info("Length of tensor quantities: {0}", m_3N);

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

    // Initialize force vectors
    spdlog::get(m_logName)->info("Initializing hydrodynamic force vectors");
    m_F_hydro          = Eigen::VectorXd::Zero(m_3N);
    m_F_hydroNoInertia = Eigen::VectorXd::Zero(m_3N);
    m_t1_Inertia       = Eigen::VectorXd::Zero(m_3N);

    spdlog::get(m_logName)->info("Initializing hydrodynamic force tensors");
    m_t2_VelGrad = Eigen::Tensor<double, 1>(m_3N);
    m_t2_VelGrad.setZero();
    m_t3_PosGrad = Eigen::Tensor<double, 1>(m_3N);
    m_t3_PosGrad.setZero();

    // Assign particle pair information
    spdlog::get(m_logName)->info("Initializing particle pair information vectors");
    m_alphaVec = Eigen::VectorXd::Zero(m_num_inter);
    m_betaVec  = Eigen::VectorXd::Zero(m_num_inter);
    m_r_mag_ab = Eigen::VectorXd::Zero(m_num_inter);
    m_r_ab     = Eigen::MatrixXd::Zero(3, m_num_inter);

    /* Fill the particle index vectors
     * Calculate ahead of time to save time during runtime
     */
    spdlog::get(m_logName)->info("Filling particle pair information tensors");
    for (int i = 0; i < m_num_inter; i++)
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
    spdlog::get(m_logName)->info("Calling update()");
    update();

    spdlog::get(m_logName)->info("Constructor complete");
    spdlog::get(m_logName)->flush();
}

potentialHydrodynamics::~potentialHydrodynamics()
{
    spdlog::get(m_logName)->info("Destructing potential hydrodynamics");
    spdlog::get(m_logName)->flush();
}

void
potentialHydrodynamics::update()
{
    calcParticleDistances();

    calcAddedMass();
    calcAddedMassGrad();
    calcHydroTensors();

    calcHydroForces();
}

void
potentialHydrodynamics::updateForcesOnly()
{
    calcHydroForces();
}

void
potentialHydrodynamics::calcParticleDistances()
{
    /* NOTE: Fill Mass matrix elements one (3 x 3) block at a time (matrix elements between
     * particles \alpha and \beta) */
    for (int i = 0; i < m_num_inter; i++)
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
    for (int k = 0; k < m_num_inter; k++)
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
        Mij.noalias() += M1_c2 * m_I3;

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
potentialHydrodynamics::calcAddedMassGrad()
{
    // eigen tensor contraction variables
    Eigen::array<Eigen::IndexPair<long>, 0> empty_index_list = {};
    Eigen::array<int, 3>                    shuffle_one_right({2, 0, 1});
    Eigen::array<int, 3>                    shuffle_one_left({1, 2, 1});

    // set matrices to zero
    m_grad_M_added.setZero();

    /* NOTE: Fill Mass matrix elements one (3 x 3 x 3) block at a time (matrix elements between
     * particles \alpha and \beta) */

    for (int k = 0; k < m_num_inter; k++)
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
            gradM1_c1 * m_system->tensI().contract(TensorCast(r_ij), empty_index_list);
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
}

void
potentialHydrodynamics::calcHydroTensors()
{
    // total mass
    m_M_total.noalias() = m_M_intrinsic;
    m_M_total.noalias() += m_M_added;
}

void
potentialHydrodynamics::calcHydroForces()
{
    // calculate requisite configurational tensors
    Eigen::MatrixXd          U_dyad_U = m_system->velocitiesParticles() * m_system->velocitiesParticles().transpose();
    Eigen::Tensor<double, 2> tens_UU  = TensorCast(U_dyad_U);

    // Term 1: - M_kj * U_dot_j
    m_t1_Inertia.noalias() = -m_M_total * m_system->accelerationsParticles(); // dim = (3N) x 1

    // Term 2: - d(M_kj)/d(R_l) * U_l * U_j;
    Eigen::array<Eigen::IndexPair<int>, 2> contract_kjl_jl = {Eigen::IndexPair<int>(1, 0), Eigen::IndexPair<int>(2, 1)};
    m_t2_VelGrad                                           = -m_grad_M_added.contract(tens_UU, contract_kjl_jl);

    // Term 3: (1/2) U_i * d(M_ij)/d(R_k) * U_j
    Eigen::array<Eigen::IndexPair<int>, 2> contract_ijk_ij = {Eigen::IndexPair<int>(0, 0), Eigen::IndexPair<int>(1, 1)};
    m_t3_PosGrad                                           = 0.50 * m_grad_M_added.contract(tens_UU, contract_ijk_ij);

    // Compute non-inertial part of force
    m_F_hydroNoInertia.noalias() = MatrixCast(m_t2_VelGrad, m_3N, 1);
    m_F_hydroNoInertia.noalias() += MatrixCast(m_t3_PosGrad, m_3N, 1);

    // Compute complete potential pressure force
    m_F_hydro.noalias() = m_t1_Inertia; // Full hydrodynamic force, dim = (3N) x 1
    m_F_hydro.noalias() += m_F_hydroNoInertia;
}
