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
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>
#define EIGEN_USE_MKL_ALL
#include <eigen3/Eigen/src/Core/util/MKL_support.h>
// STL
#include <memory>    // for std::unique_ptr and std::shared_ptr
#include <stdexcept> // std::errors
#include <string>    // std::string

/* Forward declarations */
class GSDUtil;

class systemData : public std::enable_shared_from_this<systemData>
{
  public:
    systemData(std::string inputGSDFile, std::string outputDir);

    ~systemData();

    void
    parseGSD();

    void
    check_gsd_return();

    void
    resizeTensors();

  private:
    void
    checkInput();

    // data i/o
    std::string m_inputGSDFile;
    std::string m_outputDir;

    // logging
    std::string       m_logFile;
    const std::string m_logName{"systemData"};

    // GSD data
    std::shared_ptr<GSDUtil>    m_gsdUtil;
    std::shared_ptr<gsd_handle> m_handle{new gsd_handle};
    int                         m_return_val{0};
    bool                        m_return_bool{true};
    bool                        m_GSD_parsed{false};

    // kinematics
    Eigen::VectorXd m_orientations; // unit quaternions
    Eigen::VectorXd m_positions;
    Eigen::VectorXd m_velocities;
    Eigen::VectorXd m_accelerations;

    // particle parameters
    Eigen::VectorXi m_particle_type_id; // REVIEW[epic=assumptions] {0: locater particle, 1:
                                        // constrained particle}

    // degrees of freedom
    int m_num_DoF{-1};
    int m_num_spatial_dim{-1};
    int m_num_particles{-1};
    int m_num_bodies{-1};

    // integrator
    double m_dt{-1.0};
    double m_tf{-1.0};
    double m_t{0.0};
    double m_tau{-1.0};
    int    m_timestep{-1};
    int    m_num_steps_output{-1};

    // material parameters
    double m_fluid_density{-1};
    double m_particle_density{-1};

    // potential parameters
    double m_wca_epsilon{-1};
    double m_wca_sigma{-1};

    // setters/getters
  public:
    // data i/o
    std::string
    outputDir() const
    {
        return m_outputDir;
    }

    std::string
    inputGSDFile() const
    {
        return m_inputGSDFile;
    }

    // GSD data
    std::shared_ptr<gsd_handle>
    handle() const
    {
        return m_handle;
    }

    std::shared_ptr<GSDUtil>
    gsdUtil() const
    {
        return m_gsdUtil;
    }

    bool
    gSDParsed() const
    {
        return m_GSD_parsed;
    }

    void
    setReturnVal(int return_val)
    {
        m_return_val = return_val;
    }

    void
    setReturnBool(bool return_bool)
    {
        m_return_bool = return_bool;
    }

    // kinematics
    const Eigen::VectorXd&
    orientations() const
    {
        return m_orientations;
    }
    void
    setOrientations(const Eigen::VectorXd& orientations)
    {
        m_orientations = orientations;
    }

    const Eigen::VectorXd&
    positions() const
    {
        return m_positions;
    }
    void
    setPositions(const Eigen::VectorXd& positions)
    {
        m_positions = positions;
    }

    const Eigen::VectorXd&
    velocities() const
    {
        return m_velocities;
    }
    void
    setVelocities(const Eigen::VectorXd& velocities)
    {
        m_velocities = velocities;
    }

    const Eigen::VectorXd&
    accelerations() const
    {
        return m_accelerations;
    }
    void
    setAccelerations(const Eigen::VectorXd& accelerations)
    {
        m_accelerations = accelerations;
    }

    // particle parameters
    const Eigen::VectorXi&
    particleTypeId() const
    {
        return m_particle_type_id;
    }
    void
    setParticleTypeId(const Eigen::VectorXi& particle_type_id)
    {
        m_particle_type_id = particle_type_id;
    }

    // degrees of freedom
    int
    numDoF()
    {
        return m_num_DoF;
    }

    int
    numSpatialDim() const
    {
        return m_num_spatial_dim;
    }
    void
    setNumSpatialDim(int num_spatial_dim)
    {
        m_num_spatial_dim = num_spatial_dim;
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

    int
    numBodies() const
    {
        return m_num_bodies;
    }
    void
    setNumBodies(int num_bodies)
    {
        m_num_bodies = num_bodies;
    }

    // integrator
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

    double
    tau() const
    {
        return m_tau;
    }
    void
    setTau(double tau)
    {
        m_tau = tau;
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

    // material parameters
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

    // potential parameters
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
};

#endif
