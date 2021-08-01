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

class potentialHydrodynamics
{
  public:
    potentialHydrodynamics(systemData& sys);

    ~potentialHydrodynamics();

  private:
    // classes
    systemData* system;

    // logging
    std::string m_logFile;
    std::string m_logName{"potentialHydrodynamics"};
};

#endif // BODIES_IN_POTENTIAL_FLOW_POTENTIAL_HYDRODYNAMICS_H