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
// eigen3(Linear algebra)
#define EIGEN_NO_AUTOMATIC_RESIZING
#define EIGEN_USE_MKL_ALL
#define EIGEN_USE_THREADS
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include <eigen3/unsupported/Eigen/CXX11/ThreadPool>
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
     */
    void
    update();

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
     */
    void
    calcAddedMassGrad();

    /**
     * @brief Calculates \{`m_N1`, `m_N2`, `m_N3`, `m_M_tilde`, and `m_M_tilde_tilde`\}
     *
     * @details Must call `calcTotalMass()` and assumes `systemData` rigid body motion tensors are up to date.
     *
     */
    void
    calcBodyTensors(); // TODO

    /**
     * @brief Calculates `m_F_hydro` and `m_F_hydroNoInertia`
     *
     * @details Must call `calcBodyTensors()` before.
     * Kinematics specified in `systemData` class.
     *
     */
    void
    calcHydroForces();

    // classes
    std::shared_ptr<systemData> m_system; ///< shared pointer reference to systemData class

    // logging
    std::string       m_logFile;                           ///< path of logfile for spdlog to write to
    const std::string m_logName{"potentialHydrodynamics"}; ///< filename of logfile for spdlog to write to

    // For-loop variables
    int m_num_pair_inter{-1}; ///< = s. Number of pairwise interactions to count: @f$s = 1/2 \, N \, (N - 1) @f$

    // tensor variables
    int m_3N{-1}; ///< length of tensor quantities @f$ 3 \, N @f$

    // Distance between particle pairs
    Eigen::VectorXd m_alphaVec; ///< \[s x 1\] first particle number in pairwise interactions
    Eigen::VectorXd m_betaVec;  ///< \[s x 1\] second particle number in pairwise interactions
    Eigen::VectorXd m_r_mag_ab; ///< \[s x 1\] relative displacement between particle pairs
    Eigen::MatrixXd m_r_ab;     ///< \[3 x s\] relative displacement between particle pairs

    // Identity Matrices
    Eigen::MatrixXd m_I3N;      ///< \[3N x 3N\] identity matrix
    Eigen::MatrixXd m_c1_2_I3N; ///< \[3N x 3N\] 1/2 identity matrix

    // Mass matrices
    Eigen::MatrixXd          m_M_added;      ///< \[3N x 3N\] added mass matrix
    Eigen::MatrixXd          m_M_intrinsic;  ///< \[3N x 3N\] intrinsic mass matrix
    Eigen::MatrixXd          m_M_total;      ///< \[3N x 3N\] total mass matrix
    Eigen::Tensor<double, 2> m_tens_M_total; ///< \[3N x 3N\] tensor version of total mass matrix
    Eigen::Tensor<double, 3> m_grad_M_added; ///< \[3N x 3N\] gradient of total mass matrix (only added mass components)

    // Hydrodynamic forces
    Eigen::VectorXd m_F_hydro;          ///< \[3N x 1\] total hydrodynamic force
    Eigen::VectorXd m_F_hydroNoInertia; ///< \[3N x 1\] hydrodynamic force in absence of inertial term

    Eigen::VectorXd          m_t1_Inertia; ///< \[3N x 1\] inertial term in hydrodynamic force
    Eigen::Tensor<double, 1> m_t2_VelGrad; ///< \[3N x 1\] velocity gradient term in hydrodynamic force
    Eigen::Tensor<double, 1> m_t3_PosGrad; ///< \[3N x 1\] position gradient term in hydrodynamic force

    // linear combinations of gradient of rbm and total mass matrix
    Eigen::Tensor<double, 3> m_N1; ///< \[3N x 3N x 7M\] @f$ \nabla_{\xi} \, \boldsymbol{M} @f$ TODO
    Eigen::Tensor<double, 3>
        m_N2; ///< \[7M x 3N x 7M\] @f$ \nabla_{\xi} \, \boldsymbol{A}^{\mathrm{T}} \, \boldsymbol{M} @f$ TODO
    Eigen::Tensor<double, 3> m_N3; ///< \[7M x 7M x 7M\] @f$ \nabla_{\xi} \, \boldsymbol{A}^{\mathrm{T}} \,
                                   /// \boldsymbol{M} \, \boldsymbol{A} @f$ TODO

    Eigen::Tensor<double, 2>
        m_M_tilde; ///< \[7M x 3N\] @f$ \boldsymbol{A}^{\mathrm{T}} \, \boldsymbol{M} \, \boldsymbol{A} @f$ TODO
    Eigen::Tensor<double, 2>
                    m_M_tilde_tilde;        ///< \[7M x 7M\] @f$ \boldsymbol{A}^{\mathrm{T}} \, \boldsymbol{M} @f$ TODO
    Eigen::MatrixXd m_M_tilde_tilde_matrix; ///< \[7M x 7M\] @f$ `Eigen::Matrix` form of `m_M_tilde_tilde`
    // constants
    const double m_unitSphereVol{4.0 / 3.0 * M_PI}; ///< volume of a unit sphere
    const double m_c1_2{0.50};                      ///< = 1/2
    const double m_c2_3{2.0 / 3.0};                 ///< = 2/3
    const double m_c3_2{1.50};                      ///< = 3/2
    const double m_c15_2{7.50};                     ///< = 15/2

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
    mTildeTilde()
    {
        m_M_tilde_tilde_matrix = MatrixCast(m_M_tilde_tilde, 7 * m_system->numBodies(), 7 * m_system->numBodies());

        return m_M_tilde_tilde_matrix;
    }
};

#endif // BODIES_IN_POTENTIAL_FLOW_POTENTIAL_HYDRODYNAMICS_H