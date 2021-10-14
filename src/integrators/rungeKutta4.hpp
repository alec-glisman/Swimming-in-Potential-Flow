//
// Created by Alec Glisman on 07/30/21
// Building off of GSDReader class in HOOMD-Blue
//

#ifndef BODIES_IN_POTENTIAL_FLOW_RUNGE_KUTTA_4_H
#define BODIES_IN_POTENTIAL_FLOW_RUNGE_KUTTA_4_H

#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif
/* Include all internal project dependencies */
#include <potentialHydrodynamics.hpp>
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
// Debugging
#include <iostream>

/* Forward declarations */
class systemData;

/**
 * @class rungeKutta4
 *
 * @brief Standard Runge-Kutta 4th order integration.
 *
 * @details Integrates both velocity and positional degrees of freedom simultaneously from known acceleration
 * components.
 * Orientational components are handled with quaterions and the Udwadia-Kalaba constrained Lagrangian Dynamics formalism
 * is used to ensure unitary norms of each quaternion. The true D.o.F. integrated are those of the body components and
 * individual particle components are calculated via the rigid body motion connectivity tensors.
 *
 */
class rungeKutta4
{
  public:
    /**
     * @brief Construct a new runge Kutta4 object
     *
     * @param sys systemData class to gather data from
     * @param hydro potentialHydrodynamics class to get hydrodynamic force data from
     */
    explicit rungeKutta4(std::shared_ptr<systemData> sys, std::shared_ptr<potentialHydrodynamics> hydro);

    /**
     * @brief Destroy the runge Kutta4 object
     *
     */
    ~rungeKutta4();

    void
    integrate(Eigen::ThreadPoolDevice& device);

  private:
    void
    integrateSecondOrder(Eigen::ThreadPoolDevice& device);

    void
    accelerationUpdate(Eigen::VectorXd& acc, Eigen::ThreadPoolDevice& device);

    void
    udwadiaKalaba(Eigen::VectorXd& acc);

    // classes
    /// shared pointer reference to systemData class
    std::shared_ptr<systemData> m_system;
    /// shared pointer reference to potentialHydrodynamics class
    std::shared_ptr<potentialHydrodynamics> m_potHydro;

    // logging
    /// path of logfile for spdlog to write to
    std::string m_logFile;
    /// filename of logfile for spdlog to write to
    const std::string m_logName{"rungeKutta4"};

    // time step variables
    /// (dimensional) integrator finite time step
    double m_dt{-1.0};
    /// (dimensional) 1/2 integrator finite time step
    double m_c1_2_dt{-1.0};
    /// (dimensional) 1/6 integrator finite time step
    double m_c1_6_dt{-1.0};

    // constants
    /// = 1/2
    const double m_c1_2{0.50};
    /// = 1/6
    const double m_c1_6{1.0 / 6.0};
};

#endif // BODIES_IN_POTENTIAL_FLOW_RUNGE_KUTTA_4_H