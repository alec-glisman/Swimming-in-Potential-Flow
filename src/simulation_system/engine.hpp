//
// Created by Alec Glisman on 07/31/21
//

#ifndef BODIES_IN_POTENTIAL_FLOW_ENGINE_H
#define BODIES_IN_POTENTIAL_FLOW_ENGINE_H

#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif

/* Include all internal project dependencies */
#include <potentialHydrodynamics.hpp>
#include <progressBar.hpp>
#include <rungeKutta4.hpp>
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
// STL
#include <math.h>    // isinf, sqr
#include <memory>    // for std::unique_ptr and std::shared_ptr
#include <stdexcept> // std::errors
#include <string>    // std::string

/* Forward declarations */
class systemData;

/**
 * @class engine
 *
 * @brief Manages the time integration of the simulation system as well as when to output data to GSD.
 *
 */
class engine
{
  public:
    /**
     * @brief Construct a new engine object
     *
     * @param sys systemData class to collect data from. Must be fully initialized and have GSD data loaded.
     */
    explicit engine(std::shared_ptr<systemData> sys);

    /**
     * @brief Destroy the engine object
     *
     */
    ~engine();

    /**
     * @brief Runs the simulation from @f$ t_0 @f$ to @f$ t_0f @f$.
     *
     * @details Method creates an `Eigen::ThreadPool` and `Eigen::ThreadPoolDevice` that is passed
     * to the `integrate()` method to speed up `Eigen::Tensor` computations.
     * Method also calculates the total number of integration steps required and manages the output
     * of the `progressBar` class.
     *
     */
    void
    run();

  private:
    /**
     * @brief Updates the simulation framework 1 time step.
     *
     * @param device device (CPU thread-pool or GPU) used to speed up tensor calculations
     */
    void
    integrate(Eigen::ThreadPoolDevice& device);

    // classes
    std::shared_ptr<systemData>             m_system;
    std::shared_ptr<potentialHydrodynamics> m_potHydro;
    std::shared_ptr<rungeKutta4>            m_rk4Integrator;
    std::shared_ptr<ProgressBar>            m_progressBar;

    // logging
    std::string       m_logFile;
    const std::string m_logName{"engine"};

    // eigen parallelization parameters
    /// number of physical CPU cores to use in tensor calculations
    const int m_num_physical_cores = std::thread::hardware_concurrency();

    // ProgressBar output
    /// Percentage of simulation progress at which to output a GSD frame
    const double m_outputPercentile{0.01};
};

#endif // BODIES_IN_POTENTIAL_FLOW_ENGINE_H