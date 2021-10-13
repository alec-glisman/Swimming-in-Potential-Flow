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
#define EIGEN_USE_THREADS
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include <eigen3/unsupported/Eigen/CXX11/ThreadPool>
#include <helper_eigenTensorConversion.hpp>
// STL
#include <memory>    // for std::unique_ptr and std::shared_ptr
#include <stdexcept> // std::errors
#include <string>    // std::string
#include <thread>    // std::thread::hardware_concurrency(); number of physical cores

/* Forward declarations */
class GSDUtil;

/**
 * @class systemData
 *
 * @brief Class contains all data structures relevant to the simulation framework.
 * Note that after construction, `initializeData()` must be called to properly load
 * data from a GSD file.
 *
 */
class systemData : public std::enable_shared_from_this<systemData>
{
  public:
    systemData(std::string inputGSDFile, std::string outputDir);

    ~systemData();

    /**
     * @brief Function loads data from GSD file using `GSDUtil` class.
     *
     * @details Must be called after class construction before data can properly
     * be integrated using the `engine` class.
     *
     */
    void
    initializeData();

    /**
     * @brief Updates all relevant rigid body motion tensors, respective gradients, and kinematic/Udwadia constraints
     *
     * @details Many functions are called and there is a dependency chain between them. All function calls can be run in
     * any order besides gradientChangeOfVariableTensors(), which depends on rigidBodyMotionTensors() being run first.
     * This is assuming the Udwadia linear constraint system (@f$ \mathbf{A} \, \ddot{\boldsymbol{\xi}} = \mathbf{b}
     * @f$ ) is independent of configuration.
     *
     * @param time (dimensionless) simulation time to update parameters.
     */
    void
    updateConstraints(double time);

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
    rigidBodyMotionTensors();

    void
    gradientChangeOfVariableTensors();

    void
    nMatrices();

    void
    udwadiaLinearSystem(double time);

    /* SECTION: Static functions */
    /**
     * @brief Computes the E matrix of quaternion (4-vector) input
     *
     * @param vec Input 4-vector
     * @param mat Output matrix representation
     */
    static void
    eMatrix(const Eigen::Vector4d& vec, Eigen::Matrix<double, 3, 4>& mat)
    {
        mat << -vec(1), vec(0), -vec(3), vec(2), -vec(2), vec(3), vec(0), -vec(1), -vec(3), -vec(2), vec(1), vec(0);
    };

    /**
     * @brief Function takes in vector in vector cross-product expression: $c = a \times b$
     *
     * @param vec Input 3-vector must be $a$ in above equation.
     * @param mat Matrix representation of $a \times$ operator.
     */
    static void
    crossProdMat(const Eigen::Vector3d& vec, Eigen::Matrix3d& mat)
    {
        mat << 0, -vec(2), vec(1), vec(2), 0, -vec(0), -vec(1), vec(0), 0;
    };
    /* !SECTION */

    /* ANCHOR: general attributes */
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

    /* ANCHOR: eigen parallelization parameters */
    int m_num_physical_cores = std::thread::hardware_concurrency();
    ///< number of threads in pool for eigen calculations (set to number of physical cores on machine)
    Eigen::ThreadPool m_thread_pool = Eigen::ThreadPool(m_num_physical_cores);

    /* REVIEW[epic=Change,order=0]: System specific data, change parameters stored for different
     * systems */
    // Swimming kinematic constraints
    double m_sys_spec_U0{-1.0};          ///< velocity amplitude of kinematic constraint between particle pairs
    double m_sys_spec_omega{-1.0};       ///< oscillation frequency
    double m_sys_spec_phase_shift{-1.0}; ///< phase shift (in radians) between oscillators
    double m_sys_spec_R_avg{-1.0};       ///< Time-average spatial separation between a particle pair during oscillation

    /* ANCHOR: Udwadia constraint linear system */
    /// \[M x N\] linear operator defining relationship between constraints on @f$ \ddot{\boldsymbol{\xi}} @f$.
    Eigen::MatrixXd m_Udwadia_A;
    /// \[M x 1\] Result of @f$ \mathbf{A} \, \ddot{\boldsymbol{\xi}} @f$
    Eigen::VectorXd m_Udwadia_b;

    /* ANCHOR: Tensors set in constructor */
    // "identity" tensors
    const Eigen::Matrix3d m_I = Eigen::Matrix3d::Identity(3, 3); ///< \[3 x 3\] 2nd order identity tensor
    Eigen::Matrix3d       m_I_tilde;                             ///< \[3 x 3\] reflection about z-axis tensor

    // general-use tensors
    /// \[3 x 3 x 3\] (skew-symmetric) 3rd order identity tensor
    Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3, 3>> levi_cevita;
    /// \[3 x 4 x 7\] @f$ \nabla_{\xi_{\alpha}} \boldsymbol{E}{(\boldsymbol{\theta})} @f$
    Eigen::TensorFixedSize<double, Eigen::Sizes<3, 4, 7>> kappa_tilde;

    /* ANCHOR: rigid body motion tensors */
    /// \[6M x 3N\] @f$ \boldsymbol{\Sigma} @f$ rigid body motion connectivity tensor
    Eigen::MatrixXd m_rbm_conn;
    /// \[6M x 7M\] @f$ \boldsymbol{\Psi} @f$ converts linear/quaternion body velocity D.o.F. to linear/angular
    /// velocity D.o.F.
    Eigen::MatrixXd m_psi_conv_quat_ang;
    /// \[3N x 7M\] @f$ \boldsymbol{C} @f$ converts linear/quaternion body velocity D.o.F. to linear particle
    /// velocities (NOTE: this was \f$ \boldsymbol{A} \f$ in written work)
    Eigen::MatrixXd m_C_conv_quat_part;
    /// \[3N x 7M x 7M\] @f$ \nabla_{\xi} \boldsymbol{C} @f$
    Eigen::Tensor<double, 3> m_C_conv_quat_part_grad;

    /* ANCHOR: gradient tensors in E.o.M. */

    // change of gradient variable tensors
    /// \[7M x 3N\] converts particle position D.o.F. to body position/quaternion D.o.F. (NOTE: this was
    /// @f$ \boldsymbol{\beta} @f$ in written work)
    Eigen::MatrixXd m_D_conv_quat_part;

    // linear combinations of gradient of rbm and C
    Eigen::Tensor<double, 3> m_N1; ///< \[3N x 3N x 7M\] TODO
    Eigen::Tensor<double, 3> m_N2; ///< \[7M x 3N x 7M\] TODO
    Eigen::Tensor<double, 3> m_N3; ///< \[7M x 7M x 7M\] TODO

    /* ANCHOR: kinematic vectors */
    Eigen::VectorXd m_positions_bodies;     ///< \[7M x 1\] both linear and angular D.o.F.
    Eigen::VectorXd m_velocities_bodies;    ///< \[7M x 1\] both linear and angular D.o.F.
    Eigen::VectorXd m_accelerations_bodies; ///< \[7M x 1\] both linear and angular D.o.F.

    Eigen::VectorXd m_orientations_particles;  ///< \[4N x 1\]
    Eigen::VectorXd m_positions_particles;     ///< \[3N x 1\]
    Eigen::VectorXd m_velocities_particles;    ///< \[3N x 1\]
    Eigen::VectorXd m_accelerations_particles; ///< \[3N x 1\]

    Eigen::VectorXd m_positions_locater_particles;          ///< \[3M x 1\]
    Eigen::VectorXd m_velocities_particles_articulation;    ///< \[3N x 1\]
    Eigen::VectorXd m_accelerations_particles_articulation; ///< \[3N x 1\]

    /* ANCHOR: particle parameters */
    ///\[N x 1\] REVIEW[epic=assumptions] 1) {0: constrained particle, 1: locater particle}.
    /// 2) locater particle listed first in order and denotes when to switch body index to next body.
    /// ex: (1, 0, 0, 1, 0, 0)
    Eigen::VectorXi m_particle_type_id;

    /* ANCHOR: degrees of freedom */
    int m_num_spatial_dim{-1}; ///< = 3
    int m_num_particles{-1};   ///< = N
    int m_num_bodies{-1};      ///< = M
    int m_num_DoF{-1};         ///< = 6M
    int m_num_constraints{-1}; ///< = M

    /* ANCHOR: integrator parameters */
    double m_dt{-1.0};             ///< (dimensionless) integration @f$ \Delta t @f$
    double m_tf{-1.0};             ///< (dimensionless) final simulation time
    double m_t{0.0};               ///< (dimensionless) current simulation time
    double m_tau{-1.0};            ///< simulation system characteristic timescale
    int    m_timestep{-1};         ///< current simulation timestep (iteration)
    int    m_num_steps_output{-1}; ///< how many simulation steps to output to GSD

    /* ANCHOR: material parameters */
    double m_fluid_density{-1};    ///< mass density of fluid
    double m_particle_density{-1}; ///< mass density of solid spheres

    /* ANCHOR: potential parameters */
    double m_wca_epsilon{-1}; ///< @f$ \epsilon_{\mathrm{WCA}} @f$
    double m_wca_sigma{-1};   ///< @f$ \sigma_{\mathrm{WCA}} @f$

    /* SECTION: Setters and getters */
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

    /* !SECTION */
};

#endif
