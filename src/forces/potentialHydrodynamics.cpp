//
// Created by Alec Glisman on 07/31/21
//

#include <potentialHydrodynamics.hpp>

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

    // Set identity matrices
    I3N = Eigen::MatrixXd::Identity(system->numParticles(), system->numParticles());

    // Initialize mass matrices
    spdlog::get(m_logName)->info("Initializing mass tensors");
    M_intrinsic = (system->particleDensity() * (4.0 / 3.0) * M_PI) * I3N;
    M_added     = Eigen::MatrixXd::Zero(system->numParticles(), system->numParticles());

    M_total.noalias() = M_added;
    M_total.noalias() += M_intrinsic;

    grad_M_added = Eigen::MatrixXd::Zero(system->numParticles(),
                                         system->numParticles() * system->numParticles());

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
    calcAddedMass();
    calcAddedMassGrad();

    calcHydroTensors();
    calcHydroForces();
}

void
potentialHydrodynamics::calcAddedMass()
{
    // TODO
}

void
potentialHydrodynamics::calcAddedMassGrad()
{
    // TODO
}

void
potentialHydrodynamics::calcHydroTensors()
{
    // TODO
}

void
potentialHydrodynamics::calcHydroForces()
{
    // TODO
}