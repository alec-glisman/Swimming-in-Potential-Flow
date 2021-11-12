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
#include <SystemData.hpp>

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

// FIXME: Remove
#include <iostream>

/* Forward declarations */
class SystemData;

/**
 * @class PotentialHydrodynamics
 *
 * @brief Computes hydrodynamic tensors (mass, forces) on a finite number of spheres in potential flow.
 *
 */
class PotentialHydrodynamics
{
  public:
    /**
     * @brief Construct a new potential hydrodynamics object
     *
     * @param sys SystemData class to gather data from
     */
    explicit PotentialHydrodynamics(std::shared_ptr<SystemData> sys);

    /**
     * @brief Destroy the potential Hydrodynamics object
     *
     */
    ~PotentialHydrodynamics();

    /**
     * @brief Updates all hydrodynamic quantities at current configuration (data from `SystemData`)
     *
     * @details Functions must be called in a certain grouping.
     * Functions within a group can be called in any order.
     * Groups must be called in ascending order.
     * (1) `calcParticleDistances()`.
     * (2) `calcAddedMass()`, `calcAddedMassGrad()`.
     * (3) `calcTotalMass()`.
     * (4) `calcBodyTensors()`.
     * (5) `calcHydroForces()`.
     * (6) `calcHydroEnergy()`
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
     * @details Configuration specified in `SystemData` class
     *
     */
    void
    calcParticleDistances();

    /**
     * @brief Calculates `m_M_added`
     *
     * @details Must call `calcParticleDistances()` before.
     * Configuration specified in `SystemData` class.
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
     * Configuration specified in `SystemData` class
     *
     * @param device device (CPU thread-pool or GPU) used to speed up tensor calculations
     *
     */
    void
    calcAddedMassGrad(Eigen::ThreadPoolDevice& device);

    /**
     * @brief Calculates \{`m_N1`, `m_N2`, `m_N3`, `m_M3`, and `m_M2`\}
     *
     * @details Must call `calcTotalMass()` and assumes `SystemData` rigid body motion tensors are up to date.
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
     * Kinematics specified in `SystemData` class.
     *
     * @param device device (CPU thread-pool or GPU) used to speed up tensor calculations
     *
     */
    void
    calcHydroForces(Eigen::ThreadPoolDevice& device);

    /**
     * @brief Calculates and sets energetic components in `SystemData` class
     *
     * @details Must call `calcBodyTensors()` before.
     *
     * @param device device (CPU thread-pool or GPU) used to speed up tensor calculations
     *
     */
    void
    calcHydroEnergy(Eigen::ThreadPoolDevice& device);

    // classes
    /// shared pointer reference to SystemData class
    std::shared_ptr<SystemData> m_system;

    // logging
    /// path of logfile for spdlog to write to
    std::string m_logFile;
    /// filename of logfile for spdlog to write to
    const std::string m_logName{"PotentialHydrodynamics"};

    // For-loop variables
    /// = s. Number of pairwise interactions to count: @f$s = 1/2 \, N \, (N - 1) @f$
    int m_num_pair_inter{-1};

    // tensor variables
    /// 7N length of tensor quantities
    int m_7N{-1};
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
    /// (7N x 7N) identity matrix for (3N x 3N) subset of linear elements
    Eigen::MatrixXd m_I7N_linear;
    /// (7N x 7N) identity matrix for (3N x 3N) subset of angular elements
    Eigen::MatrixXd m_I7N_angular;
    /// (7N x 7N) 1/2 identity matrix for (3N x 3N) subset of linear elements
    Eigen::MatrixXd m_c1_2_I7N_linear;

    // ANCHOR: mass matrices
    /// (7N x 7N) added mass matrix
    Eigen::MatrixXd m_M_added;
    /// (7N x 7N) intrinsic mass matrix
    Eigen::MatrixXd m_M_intrinsic;
    /// (7N x 7N) intrinsic moment of inertia matrix (for spheres)
    Eigen::MatrixXd m_J_intrinsic;
    /// (7N x 7N) total mass matrix
    Eigen::MatrixXd m_M_total;

    /// (7N x 7N) tensor version of total mass matrix
    Eigen::Tensor<double, 2> m_tens_M_total;
    /// (7N x 7N x 3N) gradient of total mass matrix (only added mass components) in particle coordinates
    Eigen::Tensor<double, 3> m_grad_M_added;

    /// (7N x 7N x 7M) gradient of total mass matrix (only added mass components) in body coordinates
    Eigen::Tensor<double, 3> m_grad_M_added_body_coords;

    // ANCHOR: Hydrodynamic forces
    /// (7M x 1) total hydrodynamic force
    Eigen::VectorXd m_F_hydro;
    /// (7M x 1) hydrodynamic force in absence of inertial term
    Eigen::VectorXd m_F_hydroNoInertia;

    // ANCHOR: linear combinations of gradient of rbm and total mass matrix
    /// (7N x 7N x 7M) @f$ \nabla_{\xi} \, \boldsymbol{M} @f$
    Eigen::Tensor<double, 3> m_N1;
    /// (7M x 7N x 7M) @f$ \nabla_{\xi} \, \boldsymbol{\zeta} \, \boldsymbol{M} @f$
    Eigen::Tensor<double, 3> m_N2;
    /// (7M x 7M x 7M) @f$ \nabla_{\xi} \, \boldsymbol{\zeta} \, \boldsymbol{M} \, \boldsymbol{\zeta}^{\mathrm{T}} @f$
    Eigen::Tensor<double, 3> m_N3;

    /// (7M x 7M) @f$ \boldsymbol{\zeta} \, \boldsymbol{M} \, \boldsymbol{\zeta}^{\mathrm{T}} @f$
    Eigen::Tensor<double, 2> m_M3;
    /// (7M x 7N) @f$ \boldsymbol{\zeta} \, \boldsymbol{M} @f$
    Eigen::Tensor<double, 2> m_M2;

    /// (7M x 7M) `Eigen::Matrix` form of `m_M3`
    Eigen::MatrixXd m_mat_M3;

    // ANCHOR: N sub-terms
    Eigen::Tensor<double, 3> N2_term1_preshuffle;
    Eigen::Tensor<double, 3> N2_term1;
    Eigen::Tensor<double, 3> N2_term2;

    Eigen::Tensor<double, 3> N3_term1_preshuffle;
    Eigen::Tensor<double, 3> N3_term2_preshuffle;
    Eigen::Tensor<double, 3> N3_term1;
    Eigen::Tensor<double, 3> N3_term2;
    Eigen::Tensor<double, 3> N3_term3;

    // ANCHOR: constants
    /// volume of a unit sphere
    const double m_unit_sphere_volume{4.0 / 3.0 * M_PI};
    /// moment of inertia for a unit sphere
    const double m_scalar_moment_inertia{0.40};
    /// REVIEW: arbitrary, positive scalar to use in extending tensors from 6N to 7N (detailed in Udwadia, Schutte. J.
    /// Appl. Mech. Jul 2010)
    /// @todo: (1, 1) entry is the added element to the mass matrix. Verify that changing this value does not affect
    /// the end result (Issue #6). Currently defaults to 1.
    const double m_scalar_w{1.00};
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

    const Eigen::Tensor<double, 3>&
    gradMAdded() const
    {
        return m_grad_M_added;
    }
};

#endif // BODIES_IN_POTENTIAL_FLOW_POTENTIAL_HYDRODYNAMICS_H
