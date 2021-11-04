//
// Created by Alec Glisman on 07/30/21
// Building off of GSDReader class in HOOMD-Blue
//

#ifndef BODIES_IN_POTENTIAL_FLOW_POTENTIAL_HYDRODYNAMICS_H
#define BODIES_IN_POTENTIAL_FLOW_POTENTIAL_HYDRODYNAMICS_H

#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif
/* Include all internal project dependencies */
#include <systemData.hpp>

/* Include all external project dependencies */
// Intel MKL
#if __has_include("mkl.h")
#define EIGEN_USE_MKL_ALL
#else
#pragma message(" !! COMPILING WITHOUT INTEL MKL OPTIMIZATIONS !! ")
#endif
// eigen3(Linear algebra)
#define EIGEN_NO_AUTOMATIC_RESIZING
#define EIGEN_USE_THREADS
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include <eigen3/unsupported/Eigen/CXX11/ThreadPool>
// eigen3 conversion between Eigen::Tensor (unsupported) and Eigen::Matrix
#include <helper_eigenTensorConversion.hpp>
// Logging
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/spdlog.h>

/* Forward declarations */
class systemData;

/**
 * @class potentialHydrodynamics
 *
 * @brief Computes hydrodynamic tensors (mass, forces) on a finite number of spheres in potential flow.
 *
 */
class potentialHydrodynamics
{
  public:
    /**
     * @brief Construct a new potential hydrodynamics object
     *
     * @param sys systemData class to gather data from
     */
    explicit potentialHydrodynamics(std::shared_ptr<systemData> sys);

    /**
     * @brief Destroy the potential Hydrodynamics object
     *
     */
    ~potentialHydrodynamics();

    /**
     * @brief Updates all hydrodynamic quantities at current configuration (data from `systemData`)
     *
     * @details Functions must be called in a certain grouping.
     * Functions within a group can be called in any order.
     * Groups must be called in ascending order.
     * (1) `calcParticleDistances()`.
     * (2) `calcAddedMass()`, `calcAddedMassGrad()`.
     * (3) `calcTotalMass()`.
     * (4) `calcBodyTensors()`.
     * (5) `calcHydroForces()`.
     *
     * @param device device (CPU thread-pool or GPU) used to speed up tensor calculations
     *
     */
    void
    update(Eigen::ThreadPoolDevice& device);

  private:
    /**
     * @brief Calculates `m_r_mag_ab` and `m_r_ab` at current configuration
     *
     * @details Configuration specified in `systemData` class
     *
     */
    void
    calcParticleDistances();

    /**
     * @brief Calculates `m_M_added`
     *
     * @details Must call `calcParticleDistances()` before.
     * Configuration specified in `systemData` class.
     *
     */
    void
    calcAddedMass();

    /**
     * @brief Calculates `m_M_added`
     *
     * @details Must call `calcAddedMass()` before to ensure added mass is up-to-date.
     *
     */
    void
    calcTotalMass();

    /**
     * @brief Calculates `m_M_added`
     *
     * @details Must call `calcParticleDistances()` before.
     * Configuration specified in `systemData` class
     *
     * @param device device (CPU thread-pool or GPU) used to speed up tensor calculations
     *
     */
    void
    calcAddedMassGrad(Eigen::ThreadPoolDevice& device);

    /**
     * @brief Calculates \{`m_N1`, `m_N2`, `m_N3`, `m_M3`, and `m_M2`\}
     *
     * @details Must call `calcTotalMass()` and assumes `systemData` rigid body motion tensors are up to date.
     *
     * @param device device (CPU thread-pool or GPU) used to speed up tensor calculations
     *
     */
    void
    calcBodyTensors(Eigen::ThreadPoolDevice& device);

    /**
     * @brief Calculates `m_F_hydro` and `m_F_hydroNoInertia`
     *
     * @details Must call `calcBodyTensors()` before.
     * Kinematics specified in `systemData` class.
     *
     * @param device device (CPU thread-pool or GPU) used to speed up tensor calculations
     *
     */
    void
    calcHydroForces(Eigen::ThreadPoolDevice& device);

    // classes
    /// shared pointer reference to systemData class
    std::shared_ptr<systemData> m_system;

    // logging
    /// path of logfile for spdlog to write to
    std::string m_logFile;
    /// filename of logfile for spdlog to write to
    const std::string m_logName{"potentialHydrodynamics"};

    // For-loop variables
    /// = s. Number of pairwise interactions to count: @f$s = 1/2 \, N \, (N - 1) @f$
    int m_num_pair_inter{-1};

    // tensor variables
    /// 6N length of tensor quantities
    int m_6N{-1}; // FIXME: Change to 7N and make sure all indexing reflects this
    /// 7M length of tensor quantities
    int m_7M{-1};

    // ANCHOR: distance between particle pairs
    /// (s x 1) first particle number in pairwise interactions
    Eigen::VectorXi m_alphaVec;
    /// (s x 1) second particle number in pairwise interactions
    Eigen::VectorXi m_betaVec;
    /// (s x 1) relative displacement between particle pairs
    Eigen::VectorXd m_r_mag_ab;
    /// (3 x s) relative displacement between particle pairs
    Eigen::MatrixXd m_r_ab;

    // ANCHOR: identity Matrices
    /// (6N x 6N) identity matrix for (3N x 3N) subset of linear elements
    Eigen::MatrixXd m_I6N_linear;
    /// (6N x 6N) identity matrix for (3N x 3N) subset of angular elements
    Eigen::MatrixXd m_I6N_angular;
    /// (6N x 6N) 1/2 identity matrix for (3N x 3N) subset of linear elements
    Eigen::MatrixXd m_c1_2_I6N_linear;

    // ANCHOR: mass matrices
    /// (6N x 6N) added mass matrix
    Eigen::MatrixXd m_M_added;
    /// (6N x 6N) intrinsic mass matrix
    Eigen::MatrixXd m_M_intrinsic;
    /// (6N x 6N) intrinsic moment of inertia matrix (for spheres)
    Eigen::MatrixXd m_J_intrinsic;
    /// (6N x 6N) total mass matrix
    Eigen::MatrixXd m_M_total;

    /// (6N x 6N) tensor version of total mass matrix
    Eigen::Tensor<double, 2> m_tens_M_total;
    /// (6N x 6N x 6N) gradient of total mass matrix (only added mass components) in particle coordinates
    Eigen::Tensor<double, 3> m_grad_M_added;

    /// (6N x 6N x 7M) gradient of total mass matrix (only added mass components) in body coordinates
    Eigen::Tensor<double, 3> m_grad_M_added_body_coords;

    // ANCHOR: Hydrodynamic forces
    /// (7M x 1) total hydrodynamic force
    Eigen::VectorXd m_F_hydro;
    /// (7M x 1) hydrodynamic force in absence of inertial term
    Eigen::VectorXd m_F_hydroNoInertia;

    // ANCHOR: linear combinations of gradient of rbm and total mass matrix
    /// (6N x 6N x 7M) @f$ \nabla_{\xi} \, \boldsymbol{M} @f$
    Eigen::Tensor<double, 3> m_N1;
    /// (7M x 6N x 7M) @f$ \nabla_{\xi} \, \boldsymbol{\zeta} \, \boldsymbol{M} @f$
    Eigen::Tensor<double, 3> m_N2;
    /// (7M x 7M x 7M) @f$ \nabla_{\xi} \, \boldsymbol{\zeta} \, \boldsymbol{M} \, \boldsymbol{\zeta}^{\mathrm{T}} @f$
    Eigen::Tensor<double, 3> m_N3;

    /// (7M x 6N) @f$ \boldsymbol{\zeta} \, \boldsymbol{M} \, \boldsymbol{\zeta}^{\mathrm{T}} @f$
    Eigen::Tensor<double, 2> m_M3;
    /// (7M x 7M) @f$ \boldsymbol{\zeta} \, \boldsymbol{M} @f$
    Eigen::Tensor<double, 2> m_M2;

    /// (7M x 7M) `Eigen::Matrix` form of `m_M3`
    Eigen::MatrixXd m_mat_M3;

    // ANCHOR: constants
    /// volume of a unit sphere
    const double m_unitSphereVol{4.0 / 3.0 * M_PI};
    /// = 1/2
    const double m_c1_2{0.50};
    /// = 3/2
    const double m_c3_2{1.50};
    /// = 15/2
    const double m_c15_2{7.50};

  public:
    const Eigen::VectorXd&
    fHydro() const
    {
        return m_F_hydro;
    }

    const Eigen::VectorXd&
    fHydroNoInertia() const
    {
        return m_F_hydroNoInertia;
    }

    const Eigen::MatrixXd&
    mTotal() const
    {
        return m_M_total;
    }

    const Eigen::MatrixXd&
    mTotalBodyCoords() const
    {
        return m_mat_M3;
    }
};

#endif // BODIES_IN_POTENTIAL_FLOW_POTENTIAL_HYDRODYNAMICS_H