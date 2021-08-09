//
// Created by Alec Glisman on 07/31/21
//

#include <potentialHydrodynamics.hpp>

// REVIEW[epic=Debug]: Uncomment line below to prevent all runtime checks from executing in debug
// #define NO_HYDRO_CHECK

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
        spdlog::get(m_logName)->info("Checking distance between pairs {0} & {1}", alphaVec[i],
                                     betaVec[i]);
        spdlog::get(m_logName)->info("Interparticle distance is", r_mag_ab[i]);

        assert(r_mag_ab[i] >= 2.0 && "Particle overlap found");
#endif
    }
}

void
potentialHydrodynamics::calcAddedMass()
{ // FIXME: Check the sign of this
    // set matrices to zero
    M_added.setZero(len_tensor, len_tensor);

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

        //! Full matrix elements for M^{(1)}_{ij}
        Eigen::Matrix3d Mij = r_ij * r_ij.transpose(); //! [1]; Outer product of \bm{r} \bm{r}
        Mij *= M1_c1;
        Mij.noalias() += M1_c2 * I3;

        //! Output added mass element
        M_added.block<3, 3>(i_part, j_part).noalias() = Mij;
    }

    M_added += M_added.transpose(); // Symmetry to get other mass elements
    // TODO: Change sign of M_added from M^{(1)}
    // TODO: Add diagonal elements
    // TODO: Multiply by factor of 1/2
    M_added *= system->fluidDensity() * unitSphereVol; // units
}

void
potentialHydrodynamics::calcAddedMassGrad()
{
    // TODO
}

void
potentialHydrodynamics::calcHydroTensors()
{
    // total mass
    M_total.noalias() = M_intrinsic;
    M_total.noalias() += M_added;

    // TODO
}

void
potentialHydrodynamics::calcHydroForces()
{
    // TODO
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