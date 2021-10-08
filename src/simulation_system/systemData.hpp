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
#define EIGEN_NO_AUTOMATIC_RESIZING
#define EIGEN_USE_MKL_ALL
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
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
    initializeData();

    void
    updateConstraints(double time);

    /**
     * Function takes in vector in vector cross-product expression: c = a \times b
     * vec must be 'a' in above equation
     * Output is the matrix representation of 'a \times' operator */
    static void
    crossProdMat(const Eigen::Vector3d& vec, Eigen::Matrix3d& mat)
    {
        mat << 0, -vec(2), vec(1), vec(2), 0, -vec(0), -vec(1), vec(0), 0;
    };

    /**
     * Function converts index b in (a, b, c) -> index b' in (a, b').
     *
     * Conversion between 3D matrix (row_idx_3d, column_idx_3d, layer_idx_3d)
     * into a flattened 2D representation. Layers are concatenated together
     * horizontally along 2D column axis to make one short and very wide two
     * dimensional matrix.
     *
     * int column_idx_3d:       column index of 3D matrix
     * int layer_idx_3d:        layer index of 3D matrix
     * int layer_deriv_dim:     spatial dimension of derivative var (x, y, z) denoted as (0, 1, 2)
     * int layer_column_width:  number of columns in a layer of 3D matrix
     */
    static int
    convert3dIdxTo2d(int column_idx_3d, int layer_idx_3d, int layer_deriv_dim,
                     int layer_column_width)
    {
        return column_idx_3d + (layer_idx_3d + layer_deriv_dim) * layer_column_width;
    };

  private:
    void
    parseGSD();

    void
    checkInput();

    void
    velocitiesArticulation(double time);

    void
    accelerationsArticulation(double time);

    void
    locaterPointLocations();

    void
    udwadiaLinearSystem(double time);

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

    /* REVIEW[epic=Change,order=0]: System specifi data, change parameters stored for different
     * systems */
    // Swimming kinematic constraints
    double m_sys_spec_U0{-1.0};
    double m_sys_spec_omega{-1.0};
    double m_sys_spec_phase_shift{-1.0};
    double m_sys_spec_R_avg{-1.0};

    // Udwadia constraint linear system
    Eigen::MatrixXd m_Udwadia_A;
    Eigen::VectorXd m_Udwadia_b;

    // "identity" tensors
    const Eigen::Matrix3d m_I = Eigen::Matrix3d::Identity(3, 3);
    Eigen::Matrix3d       m_I_tilde;
    Eigen::Matrix3d       m_I_tilde_tilde;

    // general-use tensors
    Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3, 3>> levi_cevita;
    Eigen::TensorFixedSize<double, Eigen::Sizes<3, 4, 7>> kappa_tilde;

    // kinematics
    Eigen::VectorXd m_positions_bodies;     // [7M x 1] (both linear and angular D.o.F.)
    Eigen::VectorXd m_velocities_bodies;    // [7M x 1] (both linear and angular D.o.F.)
    Eigen::VectorXd m_accelerations_bodies; // [7M x 1] (both linear and angular D.o.F.)

    Eigen::VectorXd m_orientations_particles;  // [4N x 1]
    Eigen::VectorXd m_positions_particles;     // [3N x 1]
    Eigen::VectorXd m_velocities_particles;    // [3N x 1]
    Eigen::VectorXd m_accelerations_particles; // [3N x 1]

    Eigen::VectorXd m_positions_locater_particles;          // [3M x 1]
    Eigen::VectorXd m_velocities_particles_articulation;    // [3N x 1]
    Eigen::VectorXd m_accelerations_particles_articulation; // [3N x 1]

    // particle parameters
    Eigen::VectorXi m_particle_type_id; // REVIEW[epic=assumptions] {0: locater particle,
                                        // 1: constrained particle}

    // degrees of freedom
    int m_num_spatial_dim{-1}; // = 3
    int m_num_particles{-1};   // = N
    int m_num_bodies{-1};      // = M
    int m_num_DoF{-1};         // = 6M
    int m_num_constraints{-1}; // = M

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
    int
    returnVal()
    {
        return m_return_val;
    }

    void
    setReturnBool(bool return_bool)
    {
        m_return_bool = return_bool;
    }
    bool
    returnBool()
    {
        return m_return_bool;
    }

    // kinematics
    const Eigen::VectorXd&
    orientationsParticles() const
    {
        return m_orientations_particles;
    }
    void
    setOrientationsParticles(const Eigen::VectorXd& orientations_particles)
    {
        m_orientations_particles = orientations_particles;
    }

    const Eigen::VectorXd&
    positionsParticles() const
    {
        return m_positions_particles;
    }
    void
    setPositionsParticles(const Eigen::VectorXd& positions_particles)
    {
        m_positions_particles = positions_particles;
    }

    const Eigen::VectorXd&
    velocitiesParticles() const
    {
        return m_velocities_particles;
    }
    void
    setVelocitiesParticles(const Eigen::VectorXd& velocities_particles)
    {
        m_velocities_particles = velocities_particles;
    }

    const Eigen::VectorXd&
    accelerationsParticles() const
    {
        return m_accelerations_particles;
    }
    void
    setAccelerationsParticles(const Eigen::VectorXd& accelerations_particles)
    {
        m_accelerations_particles = accelerations_particles;
    }

    Eigen::VectorXd
    positionsBodies() const
    {
        return m_positions_bodies;
    }
    void
    setPositionsBodies(const Eigen::VectorXd& positions_bodies)
    {
        m_positions_bodies = positions_bodies;
    }

    Eigen::VectorXd
    velocitiesBodies() const
    {
        return m_velocities_bodies;
    }
    void
    setVelocitiesBodies(const Eigen::VectorXd& velocities_bodies)
    {
        m_velocities_bodies = velocities_bodies;
    }

    Eigen::VectorXd
    accelerationsBodies() const
    {
        return m_accelerations_bodies;
    }
    void
    setAccelerationsBodies(const Eigen::VectorXd& accelerations_bodies)
    {
        m_accelerations_bodies = accelerations_bodies;
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

    const Eigen::VectorXd&
    velocitiesParticlesArticulation() const
    {
        return m_velocities_particles_articulation;
    }

    const Eigen::VectorXd&
    accelerationsParticlesArticulation() const
    {
        return m_accelerations_particles_articulation;
    }

    double
    sysSpecRAvg() const
    {
        return m_sys_spec_R_avg;
    }
    void
    setSysSpecRAvg(double sys_spec_RAvg)
    {
        m_sys_spec_R_avg = sys_spec_RAvg;
    }

    double
    sysSpecPhaseShift() const
    {
        return m_sys_spec_phase_shift;
    }
    void
    setSysSpecPhaseShift(double sys_spec_phaseShift)
    {
        m_sys_spec_phase_shift = sys_spec_phaseShift;
    }

    double
    sysSpecOmega() const
    {
        return m_sys_spec_omega;
    }
    void
    setSysSpecOmega(double sys_spec_omega)
    {
        m_sys_spec_omega = sys_spec_omega;
    }

    double
    sysSpecU0() const
    {
        return m_sys_spec_U0;
    }
    void
    setSysSpecU0(double sys_spec_U0)
    {
        m_sys_spec_U0 = sys_spec_U0;
    }

    const Eigen::MatrixXd&
    udwadiaA() const
    {
        return m_Udwadia_A;
    }

    const Eigen::VectorXd&
    udwadiaB() const
    {
        return m_Udwadia_b;
    }
};

#endif
