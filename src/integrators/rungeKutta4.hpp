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

    /**
     * @brief Public function to call internal numerical Runge Kutta integration functions.
     *
     * @param device `Eigen::ThreadPoolDevice` to use for `Eigen::Tensor` computations
     */
    void
    integrate(Eigen::ThreadPoolDevice& device);

  private:
    /**
     * @brief 2nd order Runge-Kutta 4th order integration. Integrates known acceleration body components into velocity
     * and then positional components.
     *
     * @details Solve system of form: @f$ y'(t) = f( y(t),  t ) @f$
     *
     * @see **Reference:**
     * https://www.physicsforums.com/threads/using-runge-kutta-method-for-position-calc.553663/post-3634957
     *
     * @param device `Eigen::ThreadPoolDevice` to use for `Eigen::Tensor` computations
     */
    void
    integrateSecondOrder(Eigen::ThreadPoolDevice& device);

    /**
     * @brief Helper function for `integrateSecondOrder()`.
     *
     * @details Function will update simulation time, then system data, then hydrodynamic data, and finally call
     * `udwadiaKalaba()` to calculate the body acceleration components at a given system time.
     *
     * @param t current integration time
     * @param pos body position vector that will be overwritten (if image system)
     * @param vel body velocity vector that will be overwritten (if image system)
     * @param acc (output) body acceleration vector that will be overwritten
     * @param device `Eigen::ThreadPoolDevice` to use for `Eigen::Tensor` computations
     */
    void
    accelerationUpdate(const double t, Eigen::VectorXd& pos, Eigen::VectorXd& vel, Eigen::VectorXd& acc,
                       Eigen::ThreadPoolDevice& device);

    /**
     * @brief Replaces 2nd 1/2 of body position and velocity D.o.F. with image of 1st 1/2 assuming the reflection plane
     * is @f$ z = 0 @f$.
     *
     * @see For discussion of how to invert z-axis of quaternion: https://stackoverflow.com/a/33999726
     *
     * @param pos (output) body position vector that will be overwritten
     * @param vel (output) body velocity vector that will be overwritten
     */
    void
    imageBodyPosVel(Eigen::VectorXd& pos, Eigen::VectorXd& vel);

    /**
     * @brief Replaces 2nd 1/2 of body acceleration D.o.F. with image of 1st 1/2 assuming the reflection plane
     * is @f$ z = 0 @f$.
     *
     * @see For discussion of how to invert z-axis of quaternion: https://stackoverflow.com/a/33999726
     *
     * @param acc (output) body acceleration vector that will be overwritten
     */
    void
    imageBodyAcc(Eigen::VectorXd& acc);

    /**
     * @brief Calculates the body acceleration components given the constraints established by the Udwadia linear
     * system (calculated in the `systemData` class).
     *
     * @see **Original paper on constrained Lagrangian dynamics:** Udwadia, Firdaus E., and Robert E. Kalaba. "A new
     * perspective on constrained motion." Proceedings of the Royal Society of London. Series A: Mathematical and
     * Physical Sciences 439.1906 (1992): 407-410.
     *
     * @see **Application of algorithm to unit quaternion vectors:** Udwadia, Firdaus E., and Aaron D. Schutte. "An
     * alternative derivation of the quaternion equations of motion for rigid-body rotational dynamics."
     * (2010): 044505.
     *
     * @param acc (output) body acceleration vector that will be overwritten
     */
    void
    udwadiaKalaba(Eigen::VectorXd& acc);

    // classes
    /// shared pointer reference to `systemData` class
    std::shared_ptr<systemData> m_system;
    /// shared pointer reference to `potentialHydrodynamics` class
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

    // tensor length variables
    /// = 7M
    int m_7M{-1};

    // constants
    /// = 1/2
    const double m_c1_2{0.50};
    /// = 1/6
    const double m_c1_6{1.0 / 6.0};
};

#endif // BODIES_IN_POTENTIAL_FLOW_RUNGE_KUTTA_4_H