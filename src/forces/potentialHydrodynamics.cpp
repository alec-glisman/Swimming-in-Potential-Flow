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
    m_num_inter = m_system->numParticles() * (m_system->numParticles() - 1) /
                  2; // Number of interactions to count
    spdlog::get(m_logName)->info("Setting number of interactions to count: {0}", m_num_inter);

    // tensor variables
    m_len_tensor = 3 * m_system->numParticles(); // length of tensor quantities
    spdlog::get(m_logName)->info("Length of tensor quantities: {0}", m_len_tensor);

    // set identity matrices
    m_I3N      = Eigen::MatrixXd::Identity(m_len_tensor, m_len_tensor);
    m_c1_2_I3N = m_I3N;
    m_c1_2_I3N *= m_c1_2;

    // Initialize mass matrices
    spdlog::get(m_logName)->info("Initializing mass tensors");
    m_M_intrinsic = m_system->particleDensity() * m_unitSphereVol * m_I3N;
    m_M_added     = Eigen::MatrixXd::Zero(m_len_tensor, m_len_tensor);

    m_M_total.noalias() = m_M_added;
    m_M_total.noalias() += m_M_intrinsic;

    m_grad_M_added = Eigen::MatrixXd::Zero(m_len_tensor, m_len_tensor * m_len_tensor);

    // Initialize force vectors
    spdlog::get(m_logName)->info("Initializing hydrodynamic force tensors");
    m_F_hydro          = Eigen::VectorXd::Zero(m_len_tensor);
    m_F_hydroNoInertia = Eigen::VectorXd::Zero(m_len_tensor);
    m_t1_Inertia       = Eigen::VectorXd::Zero(m_len_tensor);
    m_t2_VelGrad       = Eigen::VectorXd::Zero(m_len_tensor);
    m_t3_PosGrad       = Eigen::VectorXd::Zero(m_len_tensor);

    // Assign particle pair information
    spdlog::get(m_logName)->info("Initializing particle pair information tensors");
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
}

potentialHydrodynamics::~potentialHydrodynamics()
{
    spdlog::get(m_logName)->info("Destructing potential hydrodynamics");
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
        m_r_ab.col(i).noalias() = m_system->positions(Eigen::seqN(3 * m_alphaVec[i], 3));
        m_r_ab.col(i).noalias() -= m_system->positions(Eigen::seqN(3 * m_betaVec[i], 3));

        m_r_mag_ab[i] = m_r_ab.col(i).norm(); //! [1]; |r| between 2 particles

#if !defined(NDEBUG) && !defined(NO_HYDRO_CHECKS)
        spdlog::get(m_logName)->debug("Checking distance between pairs {0} & {1}", m_alphaVec[i],
                                      m_betaVec[i]);
        spdlog::get(m_logName)->debug("Interparticle distance is", m_r_mag_ab[i]);

        assert(m_r_mag_ab[i] >= 2.0 && "Particle overlap found");
#endif
    }
}

void
potentialHydrodynamics::calcAddedMass()
{
    // set matrices to zero
    m_M_added.setZero(m_len_tensor, m_len_tensor);

    /* Fill off-diagonal elements (without units ) */
    /* NOTE: Fill Mass matrix elements one (3 x 3) block at a time (matrix elements between
     * particles \alpha and \beta) */
    for (int k = 0; k < m_num_inter; k++)
    {
        //! Convert (\alpha, \beta) --> (i, j) by factor of 3
        int i_part = 3 * m_alphaVec[k];
        int j_part = 3 * m_betaVec[k];

        //! Full distance between particles \alpha and \beta
        Eigen::Vector3d r_ij     = m_r_ab.col(k); // [1]
        double          r_mag_ij = m_r_mag_ab[k]; //! [1]; |r| between 2 particles

        //! M^{(1)} Matrix Element Constants:
        double M1_c1 = -m_c3_2 / std::pow(r_mag_ij, 5); // [1]
        double M1_c2 = m_c1_2 / std::pow(r_mag_ij, 3);  // [1]

        //! Full matrix elements for M^{(1)}_{ij} (NOTE: missing factor of 1/2)
        Eigen::Matrix3d Mij = r_ij * r_ij.transpose(); //! [1]; Outer product of \bm{r} \bm{r}
        Mij *= M1_c1;
        Mij.noalias() += M1_c2 * m_I3;

        //! Output added mass element (symmetry of mass matrix)
        m_M_added.block<3, 3>(i_part, j_part).noalias() = Mij;
        m_M_added.block<3, 3>(j_part, i_part).noalias() = Mij;
    }

    /* Construct full added mass matrix
     * M = 1/2 I + M^{(1)}
     */
    m_M_added.noalias() += m_c1_2_I3N;                       // Add diagonal elements
    m_M_added *= m_system->fluidDensity() * m_unitSphereVol; // mass units
}

void
potentialHydrodynamics::calcAddedMassGrad()
{
    // set matrices to zero
    m_grad_M_added.setZero(m_len_tensor, m_len_tensor * m_len_tensor);

    /* NOTE: Fill Mass matrix elements one (3 x 3) block at a time (matrix elements between
     * particles \alpha and \beta) */
    for (int k = 0; k < m_num_inter; k++)
    {
        //! Convert (\alpha, \beta) --> (i, j) by factor of 3
        int i_part = 3 * m_alphaVec[k];
        int j_part = 3 * m_betaVec[k];

        //! Full distance between particles \alpha and \beta
        Eigen::Vector3d r_ij     = m_r_ab.col(k); // [1]
        double          r_mag_ij = m_r_mag_ab[k]; //! [1]; |r| between 2 particles

        //! Matrices to use in Calculation
        Eigen::Matrix3d delta_Ri_x, delta_Ri_y, delta_Ri_z; // [1]
        delta_Ri_x.noalias() = Eigen::Matrix3d::Zero(); //! (\delta_{j x} r_i) + (\delta_{i x} r_j)
        delta_Ri_y.noalias() = Eigen::Matrix3d::Zero(); //! (\delta_{j y} r_i) + (\delta_{i y} r_j)
        delta_Ri_z.noalias() = Eigen::Matrix3d::Zero(); //! (\delta_{j z} r_i) + (\delta_{i z} r_j)

        /* (\delta_{j x} r_i) Adds all components of \bm{r} as a column vector to the first column
         * (j=x) The indexing comes from a M_{i j, k} index scheme, where k=x in this case */
        delta_Ri_x.col(0).noalias() += r_ij;
        delta_Ri_x.row(0).noalias() += r_ij.transpose();
        delta_Ri_y.col(1).noalias() += r_ij;
        delta_Ri_y.row(1).noalias() += r_ij.transpose();
        delta_Ri_z.col(2).noalias() += r_ij;
        delta_Ri_z.row(2).noalias() += r_ij.transpose();

        //! Constants to use in Calculation
        double gradM1_c1 = -m_c3_2 / std::pow(r_mag_ij, 5); // [1]
        double gradM1_c2 = m_c15_2 / std::pow(r_mag_ij, 7); // [1]

        //! Full matrix elements for M_{ij, i}
        Eigen::Matrix3d Mij_ix, Mij_iy, Mij_iz;             // [1]
        Eigen::Matrix3d r_dyad_r = r_ij * r_ij.transpose(); // [1]; Outer product of \bm{r} \bm{r}

        // Mij_ix
        Mij_ix.noalias() = delta_Ri_x;
        Mij_ix.noalias() += r_ij[0] * m_I3;
        Mij_ix *= gradM1_c1;                                  // sub-term 1 complete
        Mij_ix.noalias() += (gradM1_c2 * r_ij[0]) * r_dyad_r; // sub-term 2 complete

        // Mij_iy
        Mij_iy.noalias() = delta_Ri_y;
        Mij_iy.noalias() += r_ij[1] * m_I3;
        Mij_iy *= gradM1_c1;                                  // sub-term 1 complete
        Mij_iy.noalias() += (gradM1_c2 * r_ij[1]) * r_dyad_r; // sub-term 2 complete

        // Mij_iz
        Mij_iz.noalias() = delta_Ri_z;
        Mij_iz.noalias() += r_ij[2] * m_I3;
        Mij_iz *= gradM1_c1;                                  // sub-term 1 complete
        Mij_iz.noalias() += (gradM1_c2 * r_ij[2]) * r_dyad_r; // sub-term 2 complete

        //! Rows to start data access at
        int row_Ri = i_part;
        int row_Rj = j_part;

        //! Columns for each block to start at: M_{ij, i}
        int col_j_dRi_x = flattenedCol(j_part, i_part, 0, m_len_tensor);
        int col_j_dRi_y = flattenedCol(j_part, i_part, 1, m_len_tensor);
        int col_j_dRi_z = flattenedCol(j_part, i_part, 2, m_len_tensor);

        //! M_{ij, j} elements
        int col_j_dRj_x = flattenedCol(j_part, j_part, 0, m_len_tensor);
        int col_j_dRj_y = flattenedCol(j_part, j_part, 1, m_len_tensor);
        int col_j_dRj_z = flattenedCol(j_part, j_part, 2, m_len_tensor);

        //! M_{ji, j} elements
        int col_i_dRj_x = flattenedCol(i_part, j_part, 0, m_len_tensor);
        int col_i_dRj_y = flattenedCol(i_part, j_part, 1, m_len_tensor);
        int col_i_dRj_z = flattenedCol(i_part, j_part, 2, m_len_tensor);

        //! M_{ji, i} elements
        int col_i_dRi_x = flattenedCol(i_part, i_part, 0, m_len_tensor);
        int col_i_dRi_y = flattenedCol(i_part, i_part, 1, m_len_tensor);
        int col_i_dRi_z = flattenedCol(i_part, i_part, 2, m_len_tensor);

        //! M_{ij, i}: Matrix Element (Anti-Symmetric upon exchange of derivative, Symmetric upon
        //! exchange of first two indices)
        m_grad_M_added.block<3, 3>(row_Ri, col_j_dRi_x).noalias() = Mij_ix; //! M_{ij, i_x}
        m_grad_M_added.block<3, 3>(row_Ri, col_j_dRi_y).noalias() = Mij_iy; //! M_{ij, i_y}
        m_grad_M_added.block<3, 3>(row_Ri, col_j_dRi_z).noalias() = Mij_iz; //! M_{ij, i_z}

        //! M_{ij, j} = - M_{ij, i}
        m_grad_M_added.block<3, 3>(row_Ri, col_j_dRj_x).noalias() = -Mij_ix;
        m_grad_M_added.block<3, 3>(row_Ri, col_j_dRj_y).noalias() = -Mij_iy;
        m_grad_M_added.block<3, 3>(row_Ri, col_j_dRj_z).noalias() = -Mij_iz;

        //! M_{ji, j} = - M_{ij, i}
        m_grad_M_added.block<3, 3>(row_Rj, col_i_dRj_x).noalias() = -Mij_ix;
        m_grad_M_added.block<3, 3>(row_Rj, col_i_dRj_y).noalias() = -Mij_iy;
        m_grad_M_added.block<3, 3>(row_Rj, col_i_dRj_z).noalias() = -Mij_iz;

        //! M_{ji, i} = M_{ij, i}
        m_grad_M_added.block<3, 3>(row_Rj, col_i_dRi_x).noalias() = Mij_ix;
        m_grad_M_added.block<3, 3>(row_Rj, col_i_dRi_y).noalias() = Mij_iy;
        m_grad_M_added.block<3, 3>(row_Rj, col_i_dRi_z).noalias() = Mij_iz;
    }

    /* Construct full added mass matrix
     * \nabla M = \nabla M^{(1)}
     */
    m_grad_M_added *= m_system->fluidDensity() * m_unitSphereVol; // mass units
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
    Eigen::MatrixXd U_dyad_U = m_system->velocities * m_system->velocities.transpose();

    //! Term 1: - M_kj * U_dot_j
    m_t1_Inertia.noalias() = -m_M_total * m_system->accelerations; //! dim = (3N) x 1

    //! Term 2: - d(M_kj)/d(R_l) * U_l * U_j;
    Eigen::MatrixXd hat_dMkj = Eigen::MatrixXd::Zero(m_len_tensor, m_len_tensor);
    // Reduce along derivative vars (index: l) via matrix-vector product for each derivative
    // variable
    for (int l = 0; l < m_len_tensor; l += 1)
    {
        hat_dMkj.noalias() += m_system->velocities[l] *
                              m_grad_M_added.block(0, l * m_len_tensor, m_len_tensor, m_len_tensor);
    }
    // Reduce along gradient (index: l) via matrix-vector product
    m_t2_VelGrad.noalias() = -hat_dMkj * m_system->velocities;

    //! Term 3: (1/2) U_i * d(M_ij)/d(R_k) * U_j
    m_t3_PosGrad.noalias() = Eigen::VectorXd::Zero(m_len_tensor);

    for (int k = 0; k < m_len_tensor; k += 3)
    { //! k: derivative variable, 3N loops
        // Unroll loop slightly (do 3 entries at a time)
        m_t3_PosGrad[k] = m_grad_M_added.block(0, k * m_len_tensor, m_len_tensor, m_len_tensor)
                              .cwiseProduct(U_dyad_U)
                              .sum();
        m_t3_PosGrad[k + 1] =
            m_grad_M_added.block(0, (k + 1) * m_len_tensor, m_len_tensor, m_len_tensor)
                .cwiseProduct(U_dyad_U)
                .sum();
        m_t3_PosGrad[k + 2] =
            m_grad_M_added.block(0, (k + 2) * m_len_tensor, m_len_tensor, m_len_tensor)
                .cwiseProduct(U_dyad_U)
                .sum();
    }
    /* Coefficient-wise operations for the sub-terms */
    m_t3_PosGrad *= m_c1_2;

    /* Compute non-inertial part of force */
    m_F_hydroNoInertia.noalias() = m_t2_VelGrad;
    m_F_hydroNoInertia.noalias() += m_t3_PosGrad;

    /* Compute complete potential pressure force */
    m_F_hydro.noalias() = m_t1_Inertia; // Full hydrodynamic force, dim = (3N) x 1
    m_F_hydro.noalias() += m_F_hydroNoInertia;
}

int
potentialHydrodynamics::flattenedCol(int realColParticle, int depthParticle, int spatialDim,
                                     int dim3N)
{
    /*
     * Converts index a in (a, b, c) to index b' in (a, b')
     *
     * Function converts between three dimensional matrix (rowParticle, realColumnParticle,
     * depthParticle) into a flattened two dimensional representation. "Layers" are concatenated
     * together horizontally to make one short and very wide two dimensional matrix.
     *
     * int realColParticle: 3 * the particle number of the col in 3D matrix
     * int depthParticle: 3 * the particle number of the 3rd dim in 3D matrix
     * int spatialDim: spatial dimension of derivative var (x, y, z) denoted as (0, 1, or 2),
     * respectively. int dim3N: nDim
     */
    return (realColParticle) + (depthParticle + spatialDim) * dim3N;
}