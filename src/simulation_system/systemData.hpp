//
// Created by Alec Glisman on 07/30/21
//

#ifndef BODIES_IN_POTENTIAL_FLOW_SYSTEM_DATA_H
#define BODIES_IN_POTENTIAL_FLOW_SYSTEM_DATA_H

#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif

/* Include all internal project dependencies */
#include <gsd.h> // GSD File

/* Include all external project dependencies */
// Logging
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/spdlog.h>
// eigen3(Linear algebra)
#include <Eigen/Core>
#include <Eigen/Eigen>
#define EIGEN_USE_MKL_ALL
#include <Eigen/src/Core/util/MKL_support.h>
// STL
#include <memory> // for std::unique_ptr and std::shared_ptr
#include <string> // std::string

class systemData
{
  public:
    systemData(std::string inputGSDFile, std::string outputDir);

    ~systemData();

  private:
    void
    check_gsd_return();

    // constructor
    std::string m_inputGSDFile;
    std::string m_outputDir;

    // logging
    std::string                     m_logFile;
    std::string                     m_logName{"systemData"};
    std::shared_ptr<spdlog::logger> m_logger;

    // GSD
    std::shared_ptr<gsd_handle> m_handle{new gsd_handle};
    int                         m_return_val{-1};

    // degrees of freedom
    int num_dim{-1};
    int num_particles{-1};

    // integrator
    double dt{-1};
    double tf{-1};
    int    num_steps_output{-1};

    // material parameters
    double fluid_density{-1};
    double particle_density{-1};

    // potential parameters
    double wca_epsilon{-1};
    double wca_sigma{-1};

    // kinematics
    Eigen::VectorXd positions;
    Eigen::VectorXd velocities;
    Eigen::VectorXd accelerations;
};

#endif
