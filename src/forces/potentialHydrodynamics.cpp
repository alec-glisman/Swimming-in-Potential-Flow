//
// Created by Alec Glisman on 07/31/21
//

#include <potentialHydrodynamics.hpp>

// REVIEW[epic=Debug]: Uncomment line below to prevent all runtime checks from executing in debug
#define NO_HYDRO_CHECK

potentialHydrodynamics::potentialHydrodynamics(systemData& sys)
{
    // save classes
    system = &sys;

    // Initialize logger
    m_logFile   = system->outputDir() + "/logs/" + m_logName + "-log.txt";
    auto logger = spdlog::basic_logger_mt(m_logName, m_logFile);
    spdlog::get(m_logName)->info("Initializing potential hydrodynamics");

    // Variables for for-loop
    num_inter = system->numParticles() * (system->numParticles() - 1) /
                2; // Number of interactions to count
    spdlog::get(m_logName)->info("Setting number of interactions to count: {0}", num_inter);

    // tensor variables
    len_tensor = 3 * system->numParticles(); // length of tensor quantities
    spdlog::get(m_logName)->info("Length of tensor quantities: {0}", len_tensor);

    // set identity matrices
    I3N = Eigen::MatrixXd::Identity(len_tensor, len_tensor);

    // Initialize mass matrices
    spdlog::get(m_logName)->info("Initializing mass tensors");
    M_intrinsic = system->particleDensity() * unitSphereVol * I3N;
    M_added     = Eigen::MatrixXd::Zero(len_tensor, len_tensor);

    M_total.noalias() = M_added;
    M_total.noalias() += M_intrinsic;

    grad_M_added = Eigen::MatrixXd::Zero(len_tensor, len_tensor * len_tensor);

    // Assign particle pair information
    spdlog::get(m_logName)->info("Initializing particle pair information tensors");
    alphaVec = Eigen::VectorXd::Zero(num_inter);
    betaVec  = Eigen::VectorXd::Zero(num_inter);
    r_mag_ab = Eigen::VectorXd::Zero(num_inter);
    r_ab     = Eigen::MatrixXd::Zero(3, num_inter);

    /* Fill the particle index vectors
     * Calculate ahead of time to save time during runtime
     */
    spdlog::get(m_logName)->info("Filling particle pair information tensors");
    for (int i = 0; i < num_inter; i++)
    {
        /* alpha and beta convention to convert from linear coordinate to ordered pair
         * @Source:
         * https://stackoverflow.com/questions/33810187/openmp-and-c-private-variables/33836073#33836073
         */
        int alpha = i / system->numParticles();
        int beta  = i % system->numParticles();

        if (beta <= alpha)
        {
            alpha = system->numParticles() - alpha - 2;
            beta  = system->numParticles() - beta - 1;
        }

        assert(alpha >= 0 && "Calculated alpha must be non-negative");
        assert(beta >= 0 && "Calculated beta must be non-negative");

        alphaVec[i] = alpha;
        betaVec[i]  = beta;
    }

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
potentialHydrodynamics::calcParticleDistances()
{
    /* NOTE: Fill Mass matrix elements one (3 x 3) block at a time (matrix elements between
     * particles \alpha and \beta) */
    for (int i = 0; i < num_inter; i++)
    {
        r_ab.col(i).noalias() = system->positions(Eigen::seqN(3 * alphaVec[i], 3));
        r_ab.col(i).noalias() -= system->positions(Eigen::seqN(3 * betaVec[i], 3));

        r_mag_ab[i] = r_ab.col(i).norm(); //! [1]; |r| between 2 particles

#if !defined(NDEBUG) && !defined(NO_HYDRO_CHECKS)
        spdlog::get(m_logName)->debug("Checking distance between pairs {0} & {1}", alphaVec[i],
                                      betaVec[i]);
        spdlog::get(m_logName)->debug("Interparticle distance is", r_mag_ab[i]);

        assert(r_mag_ab[i] >= 2.0 && "Particle overlap found");
#endif
    }
}

void
potentialHydrodynamics::calcAddedMass()
{
    // set matrices to zero
    M_added.setZero(len_tensor, len_tensor);

    /* Fill off-diagonal elements (without units ) */
    /* NOTE: Fill Mass matrix elements one (3 x 3) block at a time (matrix elements between
     * particles \alpha and \beta) */
    for (int k = 0; k < num_inter; k++)
    {
        //! Convert (\alpha, \beta) --> (i, j) by factor of 3
        int i_part = 3 * alphaVec[k];
        int j_part = 3 * betaVec[k];

        //! Full distance between particles \alpha and \beta
        Eigen::Vector3d r_ij     = r_ab.col(k); // [1]
        double          r_mag_ij = r_mag_ab[k]; //! [1]; |r| between 2 particles

        //! M^{(1)} Matrix Element Constants:
        double M1_c1 = -c3_2 / std::pow(r_mag_ij, 5); // [1]
        double M1_c2 = c1_2 / std::pow(r_mag_ij, 3);  // [1]

        //! Full matrix elements for M^{(1)}_{ij} (NOTE: missing factor of 1/2)
        Eigen::Matrix3d Mij = r_ij * r_ij.transpose(); //! [1]; Outer product of \bm{r} \bm{r}
        Mij *= M1_c1;
        Mij.noalias() += M1_c2 * I3;

        //! Output added mass element (symmetry of mass matrix)
        M_added.block<3, 3>(i_part, j_part).noalias() = Mij;
        M_added.block<3, 3>(j_part, i_part).noalias() = Mij;
    }

    /* Construct full added mass matrix
     * M = 1/2 I + M^{(1)}
     */
    M_added.noalias() += I3N;                                 // Add diagonal elements
    M_added *= c1_2 * system->fluidDensity() * unitSphereVol; // factor of 1/2 and mass units
}

void
potentialHydrodynamics::calcAddedMassGrad()
{
    // set matrices to zero
    grad_M_added.setZero(len_tensor, len_tensor * len_tensor);

    /* NOTE: Fill Mass matrix elements one (3 x 3) block at a time (matrix elements between
     * particles \alpha and \beta) */
    for (int k = 0; k < num_inter; k++)
    {
        //! Convert (\alpha, \beta) --> (i, j) by factor of 3
        int i_part = 3 * alphaVec[k];
        int j_part = 3 * betaVec[k];

        //! Full distance between particles \alpha and \beta
        Eigen::Vector3d r_ij     = r_ab.col(k); // [1]
        double          r_mag_ij = r_mag_ab[k]; //! [1]; |r| between 2 particles

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
        double gradM1_c1 = -c3_2 / std::pow(r_mag_ij, 5); // [1]
        double gradM1_c2 = c15_2 / std::pow(r_mag_ij, 7); // [1]

        //! Full matrix elements for M_{ij, i}
        Eigen::Matrix3d Mij_ix, Mij_iy, Mij_iz;             // [1]
        Eigen::Matrix3d r_dyad_r = r_ij * r_ij.transpose(); // [1]; Outer product of \bm{r} \bm{r}

        // Mij_ix
        Mij_ix.noalias() = delta_Ri_x;
        Mij_ix.noalias() += r_ij[0] * I3;
        Mij_ix *= gradM1_c1;                                  // sub-term 1 complete
        Mij_ix.noalias() += (gradM1_c2 * r_ij[0]) * r_dyad_r; // sub-term 2 complete

        // Mij_iy
        Mij_iy.noalias() = delta_Ri_y;
        Mij_iy.noalias() += r_ij[1] * I3;
        Mij_iy *= gradM1_c1;                                  // sub-term 1 complete
        Mij_iy.noalias() += (gradM1_c2 * r_ij[1]) * r_dyad_r; // sub-term 2 complete

        // Mij_iz
        Mij_iz.noalias() = delta_Ri_z;
        Mij_iz.noalias() += r_ij[2] * I3;
        Mij_iz *= gradM1_c1;                                  // sub-term 1 complete
        Mij_iz.noalias() += (gradM1_c2 * r_ij[2]) * r_dyad_r; // sub-term 2 complete

        //! Rows to start data access at
        int row_Ri = i_part;
        int row_Rj = j_part;

        //! Columns for each block to start at: M_{ij, i}
        int col_j_dRi_x = flattenedCol(j_part, i_part, 0, len_tensor);
        int col_j_dRi_y = flattenedCol(j_part, i_part, 1, len_tensor);
        int col_j_dRi_z = flattenedCol(j_part, i_part, 2, len_tensor);

        //! M_{ij, j} = - M_{ij, i} -- Anti-symmetry under exchange of derivative variable, due to j
        //! particle
        int col_j_dRj_x = flattenedCol(j_part, j_part, 0, len_tensor);
        int col_j_dRj_y = flattenedCol(j_part, j_part, 1, len_tensor);
        int col_j_dRj_z = flattenedCol(j_part, j_part, 2, len_tensor);

        //! M_{ji, j} = - M_{ji, i} -- Anti-Symmetry upon exchange of derivative variable
        int col_i_dRj_x = flattenedCol(i_part, j_part, 0, len_tensor);
        int col_i_dRj_y = flattenedCol(i_part, j_part, 1, len_tensor);
        int col_i_dRj_z = flattenedCol(i_part, j_part, 2, len_tensor);

        //! M_{ji, i} =   M_{ij, i} -- Symmetry under exchange of mass matrix elements
        int col_i_dRi_x = flattenedCol(i_part, i_part, 0, len_tensor);
        int col_i_dRi_y = flattenedCol(i_part, i_part, 1, len_tensor);
        int col_i_dRi_z = flattenedCol(i_part, i_part, 2, len_tensor);

        //! M_{ij, i}: Matrix Element (Anti-Symmetric upon exchange of derivative, Symmetric upon
        //! exchange of first two indices)
        grad_M_added.block<3, 3>(row_Ri, col_j_dRi_x).noalias() = Mij_ix; //! M_{ij, i_x}
        grad_M_added.block<3, 3>(row_Ri, col_j_dRi_y).noalias() = Mij_iy; //! M_{ij, i_y}
        grad_M_added.block<3, 3>(row_Ri, col_j_dRi_z).noalias() = Mij_iz; //! M_{ij, i_z}

        //! M_{ij, j} = - M_{ij, i}
        grad_M_added.block<3, 3>(row_Ri, col_j_dRj_x).noalias() = -Mij_ix;
        grad_M_added.block<3, 3>(row_Ri, col_j_dRj_y).noalias() = -Mij_iy;
        grad_M_added.block<3, 3>(row_Ri, col_j_dRj_z).noalias() = -Mij_iz;

        //! M_{ji, j} = - M_{ij, i}
        grad_M_added.block<3, 3>(row_Rj, col_i_dRj_x).noalias() = -Mij_ix;
        grad_M_added.block<3, 3>(row_Rj, col_i_dRj_y).noalias() = -Mij_iy;
        grad_M_added.block<3, 3>(row_Rj, col_i_dRj_z).noalias() = -Mij_iz;

        //! M_{ji, i} = M_{ij, i}
        grad_M_added.block<3, 3>(row_Rj, col_i_dRi_x).noalias() = Mij_ix;
        grad_M_added.block<3, 3>(row_Rj, col_i_dRi_y).noalias() = Mij_iy;
        grad_M_added.block<3, 3>(row_Rj, col_i_dRi_z).noalias() = Mij_iz;
    }

    /* Construct full added mass matrix
     * \nabla M = \nabla M^{(1)}
     */
    grad_M_added *= c1_2 * system->fluidDensity() * unitSphereVol; // factor of 1/2 and mass units
}

void
potentialHydrodynamics::calcHydroTensors()
{
    // total mass
    M_total.noalias() = M_intrinsic;
    M_total.noalias() += M_added;
}

void
potentialHydrodynamics::calcHydroForces()
{
    // calculate requisite configurational tensors
    Eigen::MatrixXd U_dyad_U = system->velocities * system->velocities.transpose();

    //! Term 1: - M_kj * U_dot_j
    Eigen::VectorXd t1_Intertia = -M_total * system->accelerations; //! dim = (3N) x 1

    //! Term 2: - d(M_kj)/d(R_l) * U_l * U_j;
    Eigen::MatrixXd hat_dMkj = Eigen::MatrixXd::Zero(len_tensor, len_tensor);
    // Reduce along derivative vars (index: l) via matrix-vector product for each derivative
    // variable
    for (int l = 0; l < len_tensor; l += 1)
    {
        hat_dMkj.noalias() +=
            system->velocities[l] * grad_M_added.block(0, l * len_tensor, len_tensor, len_tensor);
    }
    // Reduce along gradient (index: l) via matrix-vector product
    Eigen::VectorXd t2_VelGrad = -hat_dMkj * system->velocities;

    //! Term 3: (1/2) U_i * d(M_ij)/d(R_k) * U_j
    Eigen::VectorXd t3_PosGrad = Eigen::VectorXd::Zero(len_tensor);

    for (int k = 0; k < len_tensor; k += 3)
    { //! k: derivative variable, 3N loops
        // Unroll loop slightly (do 3 entries at a time)
        t3_PosGrad[k] = grad_M_added.block(0, k * len_tensor, len_tensor, len_tensor)
                            .cwiseProduct(U_dyad_U)
                            .sum();
        t3_PosGrad[k + 1] = grad_M_added.block(0, (k + 1) * len_tensor, len_tensor, len_tensor)
                                .cwiseProduct(U_dyad_U)
                                .sum();
        t3_PosGrad[k + 2] = grad_M_added.block(0, (k + 2) * len_tensor, len_tensor, len_tensor)
                                .cwiseProduct(U_dyad_U)
                                .sum();
    }
    /* Coefficient-wise operations for the sub-terms */
    t3_PosGrad *= c1_2;

    /* Compute non-inertial part of force */
    F_hydroNoInertia.noalias() = t2_VelGrad;
    F_hydroNoInertia.noalias() += t3_PosGrad;

    /* Compute complete potential pressure force */
    F_hydro.noalias() = t1_Intertia; // Full hydrodynamic force, dim = (3N) x 1
    F_hydro.noalias() += F_hydroNoInertia;
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