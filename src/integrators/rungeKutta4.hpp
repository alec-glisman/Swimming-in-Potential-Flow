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
// eigen3 (Linear algebra)
#include <Eigen/Core>
#include <Eigen/Eigen>
#define EIGEN_USE_MKL_ALL
#include <Eigen/src/Core/util/MKL_support.h>
// Logging
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/spdlog.h>

/* Forward declarations */
class systemData;

class rungeKutta4
{
  public:
    rungeKutta4(systemData& sys);

    rungeKutta4(systemData& sys, potentialHydrodynamics& hydro);

    ~rungeKutta4();

    void
    integrate();

  private:
    void
    acceleration_update();

    // classes
    systemData*             system;
    potentialHydrodynamics* potHydro;

    // logging
    std::string m_logFile;
    std::string m_logName{"rungeKutta4"};
};

#endif // BODIES_IN_POTENTIAL_FLOW_RUNGE_KUTTA_4_H