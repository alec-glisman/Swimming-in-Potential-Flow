//
// Created by Alec Glisman on 07/31/21
//

#ifndef BODIES_IN_POTENTIAL_FLOW_ENGINE_H
#define BODIES_IN_POTENTIAL_FLOW_ENGINE_H

#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif

/* Include all internal project dependencies */
#include <rungeKutta4.hpp>
#include <systemData.hpp>

/* Include all external project dependencies */
// eigen3 (Linear algebra)
#include <Eigen/Core>
#include <Eigen/Eigen>
#define EIGEN_USE_MKL_ALL
#include <Eigen/src/Core/util/MKL_support.h>
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
class rungeKutta;

class engine
{
  public:
    engine(systemData& sys);

    ~engine();

    void
    run();

  private:
    void
    acceleration_update();

    void
    integrate();

    // classes
    systemData*                  system;
    std::shared_ptr<rungeKutta4> rk4Integrator;

    // logging
    std::string m_logFile;
    std::string m_logName{"engine"};
};

#endif // BODIES_IN_POTENTIAL_FLOW_ENGINE_H