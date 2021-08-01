//
// Created by Alec Glisman on 07/30/21
//

#ifndef BODIES_IN_POTENTIAL_FLOW_SYSTEM_DATA_H
#define BODIES_IN_POTENTIAL_FLOW_SYSTEM_DATA_H

#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif

/* Include all internal project dependencies */
#include <GSDUtil.hpp> // GSD parser
#include <gsd.h>       // GSD File

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
class GSDUtil;

class systemData
{
  public:
    systemData(std::string inputGSDFile, std::string outputDir);

    ~systemData();

    void
    check_gsd_return();

    void
    resizeTensors();

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

    int
    numDim()
    {
        return m_num_dim;
    }
    void
    setNumDim(int num_dim)
    {
        m_num_dim = num_dim;
    }

    int
    numParticles() const
    {
        return m_num_particles;
    }
    void
    setNumParticles(int num_particles)
    {
        m_num_particles = num_particles;
    }

    double
    dt() const
    {
        return m_dt;
    }
    void
    setDt(double dt)
    {
        m_dt = dt;
    }

    std::shared_ptr<Eigen::VectorXd>
    positions() const
    {
        return m_positions;
    }
    void
    setPositions(const std::shared_ptr<Eigen::VectorXd>& positions)
    {
        m_positions = positions;
    }

    std::shared_ptr<Eigen::VectorXd>
    velocities() const
    {
        return m_velocities;
    }
    void
    setVelocities(const std::shared_ptr<Eigen::VectorXd>& velocities)
    {
        m_velocities = velocities;
    }

    std::shared_ptr<Eigen::VectorXd>
    accelerations() const
    {
        return m_accelerations;
    }
    void
    setAccelerations(const std::shared_ptr<Eigen::VectorXd>& accelerations)
    {
        m_accelerations = accelerations;
    }

    double
    t() const
    {
        return m_t;
    }
    void
    setT(double t)
    {
        m_t = t;
    }

    double
    tf() const
    {
        return m_tf;
    }
    void
    setTf(double tf)
    {
        m_tf = tf;
    }

    int
    numStepsOutput() const
    {
        return m_num_steps_output;
    }
    void
    setNumStepsOutput(int num_steps_output)
    {
        m_num_steps_output = num_steps_output;
    }

    double
    fluidDensity() const
    {
        return m_fluid_density;
    }
    void
    setFluidDensity(double fluid_density)
    {
        m_fluid_density = fluid_density;
    }

    double
    particleDensity() const
    {
        return m_particle_density;
    }
    void
    setParticleDensity(double particle_density)
    {
        m_particle_density = particle_density;
    }

    double
    wcaSigma() const
    {
        return m_wca_sigma;
    }
    void
    setWcaSigma(double wca_sigma)
    {
        m_wca_sigma = wca_sigma;
    }

    double
    wcaEpsilon() const
    {
        return m_wca_epsilon;
    }
    void
    setWcaEpsilon(double wca_epsilon)
    {
        m_wca_epsilon = wca_epsilon;
    }

  private:
    // classes
    std::shared_ptr<GSDUtil> m_gsdUtil;

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
    std::shared_ptr<Eigen::VectorXd> m_positions;
    std::shared_ptr<Eigen::VectorXd> m_velocities;
    std::shared_ptr<Eigen::VectorXd> m_accelerations;
};

#endif
