//
// Created by Alec Glisman on 07/30/21
//

#ifndef BODIES_IN_POTENTIAL_FLOW_SYSTEM_DATA_H
#define BODIES_IN_POTENTIAL_FLOW_SYSTEM_DATA_H

#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif

/* Include all internal project dependencies */
#include <GSDReader.hpp> // GSD parser
#include <gsd.h>         // GSD File

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
#include <memory>    // for std::unique_ptr and std::shared_ptr
#include <stdexcept> // std::errors
#include <string>    // std::string

/* Forward declarations */
class GSDReader;

class systemData
{
  public:
    systemData(std::string inputGSDFile, std::string outputDir);

    ~systemData();

    void
    check_gsd_return();

    std::string
    outputDir() const
    {
        return m_outputDir;
    }

    void
    setReturnVal(int return_val)
    {
        m_return_val = return_val;
    }

    std::shared_ptr<gsd_handle>
    handle() const
    {
        return m_handle;
    }

    int
    timestep() const
    {
        return m_timestep;
    }
    void
    setTimestep(int timestep)
    {
        m_timestep = timestep;
    }

    void
    setReturnBool(bool return_bool)
    {
        m_return_bool = return_bool;
    }

    std::string
    inputGSDFile() const
    {
        return m_inputGSDFile;
    }

  private:
    // classes
    std::shared_ptr<GSDReader> gsdReader;

    // constructor
    std::string m_inputGSDFile;
    std::string m_outputDir;

    // logging
    std::string m_logFile;
    std::string m_logName{"systemData"};

    // GSD
    std::shared_ptr<gsd_handle> m_handle{new gsd_handle};
    int                         m_return_val{0};
    bool                        m_return_bool{true};

    // degrees of freedom
    int m_num_dim{-1};
    int m_num_particles{-1};

    // integrator
    double m_dt{-1};
    double m_tf{-1};
    double m_t{0.0};
    int    m_timestep{-1};
    int    m_num_steps_output{-1};

    // material parameters
    double m_fluid_density{-1};
    double m_particle_density{-1};

    // potential parameters
    double m_wca_epsilon{-1};
    double m_wca_sigma{-1};

    // kinematics
    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_velocities;
    Eigen::VectorXd m_accelerations;
};

#endif
