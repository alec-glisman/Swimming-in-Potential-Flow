//
// Created by Alec Glisman on 07/30/21
//

#ifndef BODIES_IN_POTENTIAL_FLOW_SYSTEM_DATA_H
#define BODIES_IN_POTENTIAL_FLOW_SYSTEM_DATA_H

/* SECTION: Header */
#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif

/* Include all internal project dependencies */
#include <GSDUtil.hpp> // GSD parser
#include <gsd.h>       // GSD File

/* Include all external project dependencies */
// Intel MKL
#if __has_include("mkl.h")
#define EIGEN_USE_MKL_ALL
#else
#pragma message(" !! COMPILING WITHOUT INTEL MKL OPTIMIZATIONS !! ")
#endif
// eigen3(Linear algebra)
#define EIGEN_NO_AUTOMATIC_RESIZING
#define EIGEN_USE_THREADS
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include <eigen3/unsupported/Eigen/CXX11/ThreadPool>
// eigen3 conversion between Eigen::Tensor (unsupported) and Eigen::Matrix
#include <helper_eigenTensorConversion.hpp>
// Logging
#include <spdlog/fmt/ostr.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/spdlog.h>
// STL
#include <iostream>  // std::cout; std::endl;
#include <memory>    // for std::unique_ptr; std::shared_ptr
#include <stdexcept> // std::errors
#include <string>    // std::string
#include <thread>    // std::thread::hardware_concurrency(); number of physical cores

/* Forward declarations */
class GSDUtil;
/* !SECTION (Header) */

/**
 * @class SystemData
 *
 * @brief Class contains all data structures relevant to the simulation framework.
 * Note that after construction, `initializeData()` must be called to properly load
 * data from a GSD file.
 *
 * @see For a visualization of unit quaternion rotation: https://quaternions.online/
 *
 */
class SystemData : public std::enable_shared_from_this<SystemData>
{
    /* SECTION: Public methods */
  public:
    /**
     * @brief Construct a new system Data object. Default constructor.
     *
     */
    SystemData() = default;

    /**
     * @brief Construct a new system Data object
     *
     * @param inputGSDFile string path to (already created and set-up) GSD frame
     * @param outputDir string path to output directory for I/O
     */
    SystemData(std::string inputGSDFile, std::string outputDir);

    /**
     * @brief Destroy the system Data object
     *
     */
    ~SystemData();

    /**
     * @brief Function loads data from GSD file using `GSDUtil` class.
     *
     * @details Must be called after class construction before data can properly
     * be integrated using the `Engine` class.
     */
    void
    initializeData();

    /**
     * @brief Logs private attributes to logfile
     *
     */
    void
    logData();

    /**
     * @brief Updates all relevant rigid body motion tensors, respective gradients, and kinematic/Udwadia constraints.
     * Assumes `m_t` is current simulation time to update variables at.
     *
     * @details Many functions are called and there is a dependency chain between them.
     * The functions are grouped in the implementation and clearly denote the relative ordering that must occur.
     *
     * @param device `Eigen::ThreadPoolDevice` to use for `Eigen::Tensor` computations
     */
    void
    update(Eigen::ThreadPoolDevice& device);

    /**
     * @brief Normalizes all body coordinate quaternions.
     *
     * @details Modifies `m_positions_bodies`
     *
     */
    void
    normalizeQuaternions();
    /* !SECTION (Public methods) */

    /* SECTION: Private methods */
  private:
    /**
     * @brief Passes `this` (`SystemData` shared pointer instance) into `GSDUtil` constructor to load data from input
     * GSD file into `this`.
     *
     * @details Function called by `initializeData()` and is the reason that data can only be initialized after
     * constructor has completed. This stems from an issue in the following line of code: `m_gsdUtil    =
     * std::make_shared<GSDUtil>(shared_from_this());`. From the C++ standard, the `shared_from_this()` function can
     * only execute once the object referenced by `this` has been fully constructed. If I could get around this issue,
     * then only the class constructor would need to be called. This is safer and nicer. The code currently modifies the
     * attribute `m_GSD_parsed` to true to give some check that the system has been properly initialized. The `Engine`
     * class also checks that this boolean has been set.
     *
     * @todo (future) Find a better way to handle passing of `this` to `GSDUtil` so that this can be done in
     * `SystemData` constructor.
     */
    void
    parseGSD();

    /**
     * @brief Various assertions to check data loaded from GSD is at least physical.
     *
     * @details This is by no means exhaustive and should not be taken to mean data has been loaded correctly or as
     * expected.
     */
    void
    checkInput();

    /* SECTION: Constraints */
    /**
     * @brief Computes the articulation (linear) positions of the particles relative to their respective locater
     * points.
     *
     * @review_swimmer Change assignment of `m_positions_particles_articulation` for
     */
    void
    positionsArticulation();

    /**
     * @brief Computes the articulation (linear) velocities of the particles relative to their respective locater
     * points.
     *
     * @review_swimmer Change assignment of `m_velocities_particles_articulation` for
     */
    void
    velocitiesArticulation();

    /**
     * @brief Computes the articulation (linear) accelerations of the particles relative to their respective locater
     * points.
     *
     * @review_swimmer Change assignment of `m_accelerations_particles_articulation` for
     * different systems
     */
    void
    accelerationsArticulation();

    /**
     * @brief Computes the Udwadia constraint matrix and vector to uphold the quaternion unitary norm for each body.
     *
     * @review_swimmer Change assignment of `m_Udwadia_A` and `m_Udwadia_B` for different systems.
     *
     * @see **Original paper on constrained Lagrangian dynamics:** Udwadia, Firdaus E., and Robert E. Kalaba. "A new
     * perspective on constrained motion." Proceedings of the Royal Society of London. Series A: Mathematical and
     * Physical Sciences 439.1906 (1992): 407-410.
     *
     * @see **Application of algorithm to unit quaternion vectors:** Udwadia, Firdaus E., and Aaron D. Schutte. "An
     * alternative derivation of the quaternion equations of motion for rigid-body rotational dynamics."
     * (2010): 044505.
     */
    void
    udwadiaLinearSystem();
    /* !SECTION (Constraints) */

    /* SECTION: Rigid body motion */
    /**
     * @brief Computes \Sigma_{i \alpha} matrix for given particle number and body number.
     * Sigma matrix represents the transformation of coordinates from (linear/angular) body
     * coordinates to (linear/angular) particle relative configuration position coordinates.
     *
     * @details Modifies `m_rbm_conn`
     *
     * @param particle_id Particle number (alpha, i is body number)
     */
    void
    rbmMatrixElement(const int particle_id);

    /**
     * @brief Computes \Psi matrix for given body number.
     * \Psi matrix represents the transformation of coordinates from (linear/angular) body
     * coordinates to (linear/quaternion) body coordinates.
     *
     * @details Modifies `m_psi_conv_quat_ang`
     *
     * @param body_id Body number
     */
    void
    psiMatrixElement(const int body_id);

    /**
     * @brief Computes \chi_{i \alpha} matrix for given particle number and body number.
     * chi matrix represents the transformation of gradient coordinates from (linear/quaternion) body
     * coordinates to (linear) particle relative configuration position coordinates.
     *
     * @details Modifies `m_chi`
     *
     * @see For Wikipedia typeset version of rotated position w.r.t. body quaternion:
     * https://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation#Differentiation_with_respect_to_the_rotation_quaternion
     *
     * @see For detailed derivation: Lee, Byung-Uk (1991), "Differentiation with Quaternions, Appendix B",
     * Ph. D. Thesis, Stereo Matching of Skull Landmarks, Stanford University: 57–58
     * http://home.ewha.ac.kr/~bulee/quaternion.pdf
     *
     * @param particle_id Particle number (alpha, i is body number)
     */
    void
    chiMatrixElement(const int particle_id);

    /**
     * @brief Computes \nabla_{\xi} \Sigma_{i \alpha} tensor for given particle number and body number.
     *
     * @details Modifies `m_tens_grad_rbm_conn`
     * Some variables must be pre-computed:
     *     `m_rbm_conn`
     *     `m_psi_conv_quat_ang`
     *     `m_chi`
     *
     * @param particle_id Particle number (alpha, i is body number)
     */
    void
    gradRbmConnTensorElement(const int particle_id, Eigen::ThreadPoolDevice& device);

    /**
     * @brief Computes the rigid body motion connectivity tensors.
     *
     * @details This function computes the tensors required to convert body locater kinematics into particle
     * kinematics. It relies on `m_particle_type_id` being set with the correct convention.
     *
     * @see **Analogous paper for swimming rbm tensors** (NOTE: they also included torques associated with spherical
     * rotation about internal axes. We do not worry about that here due to no-flux boundary conditions on surfaces
     * rather than no-slip boundary conditions): Swan, James W., et al. "Modeling hydrodynamic self-propulsion with
     * Stokesian Dynamics. Or teaching Stokesian Dynamics to swim." Physics of Fluids 23.7 (2011): 071901.
     *
     * @param device `Eigen::ThreadPoolDevice` to use for `Eigen::Tensor` computations
     */
    void
    rigidBodyMotionTensors(Eigen::ThreadPoolDevice& device);

    /**
     * @brief Computes the gradients rigid body motion connectivity tensors and tensors associated with change
     * of variable degrees of freedom.
     *
     * @details This function computes the tensors required to convert gradient body locater kinematics into
     * particle kinematics. It relies on `m_particle_type_id` being set with the correct convention.
     *
     * @param device `Eigen::ThreadPoolDevice` to use for `Eigen::Tensor` computations
     */
    void
    gradientChangeOfVariableTensors(Eigen::ThreadPoolDevice& device);
    /* !SECTION (Rigid body motion) */

    /* SECTION: Convert between body and particle degrees of freedom */
    /**
     * @brief Computes the orientation (unit vector) of a particle with respect to its locater point
     *
     */
    void
    convertBody2ParticleOrient();

    /**
     * @brief Computes the positions of all particles from given locater positions and body orientations
     *
     */
    void
    convertBody2ParticlePos();

    /**
     * @brief Computes the particle D.o.F. from the body D.o.F.
     *
     * @details Currently only computes the velocity and acceleration D.o.F. components
     *
     * @param device `Eigen::ThreadPoolDevice` to use for `Eigen::Tensor` computations
     */
    void
    convertBody2ParticleVelAcc(Eigen::ThreadPoolDevice& device);
    /* !SECTION (Convert between body and particle degrees of freedom) */

    /* SECTION: Static methods */
    /**
     * @brief Function takes in vector in vector cross-product expression: @f$ c = a \times b @f$
     *
     * @param vec Input 3-vector must be @f$ a @f$ in above equation.
     * @param mat Matrix representation of @f$ a \times @f$ operator.
     */
    static inline void
    crossProdMat(const Eigen::Vector3d& vec, Eigen::Matrix3d& mat)
    {
        // clang-format off
        mat << 
            0, -vec(2), vec(1), 
            vec(2), 0, -vec(0), 
            -vec(1), vec(0), 0;
        // clang-format on
    };

    /**
     * @brief Computes the E matrix of quaternion (4-vector) input
     *
     * @see Eq. (2) of
     * Udwadia, Firdaus E., and Aaron D. Schutte. "An alternative derivation of the
     * quaternion equations of motion for rigid-body rotational dynamics." (2010)
     *
     * @param theta Input 4-vector (unit quaternion)
     * @param E_body Output matrix representation
     */
    static inline void
    eMatrix(const Eigen::Vector4d& theta, Eigen::Matrix4d& E_body)
    {
        // clang-format off
        E_body << 
            theta(0), theta(1), theta(2), theta(3), 
            -theta(1), theta(0), theta(3), -theta(2), 
            -theta(2), -theta(3), theta(0), theta(1), 
            -theta(3), theta(2), -theta(1), theta(0);
        // clang-format on
    };
    /* !SECTION (Static methods) */
    /* !SECTION (Private methods) */

    /* SECTION: Friend classes */
    friend class TestSystemData;
    /* !SECTION (Friend classes) */

    /* SECTION: Attributes */
    /* ANCHOR: Simulation hyperparameters */
    /// If the simulation system is constrained to be that of an image (neglect second 1/2 of DoF)
    /// @review_swimmer change if not using image-system constraints
    bool m_image_system{false};

    /* ANCHOR: general attributes */
    // data i/o
    std::string m_inputGSDFile;
    std::string m_outputDir;

    // logging
    /// path of logfile for spdlog to write to
    std::string m_logFile;
    /// filename of logfile for spdlog to write to
    const std::string m_logName{"SystemData"};

    // GSD data
    /// shared pointer reference to `GSDUtil` class
    std::shared_ptr<GSDUtil>    m_gsdUtil;
    std::shared_ptr<gsd_handle> m_handle{new gsd_handle};
    /// defaults to GSD_SUCCESS return value
    int m_return_val{0};
    /// defaults to successful parse GSD flag
    bool m_return_bool{true};
    // defaults to no GSD being parsed. Must be changed using `parseGSD()`
    bool m_GSD_parsed{false};

    /* ANCHOR: System specific data, change parameters stored for different systems */
    /**
     * @brief Swimming kinematic constraint velocity amplitude of kinematic constraint between particle pairs
     *
     * @review_swimmer Change parameters stored for different systems
     *
     */
    double m_sys_spec_U0{-1.0};
    /**
     * @brief Swimming kinematic constraint oscillation frequency
     *
     * @review_swimmer Change parameters stored for different systems
     *
     */
    double m_sys_spec_omega{-1.0};
    /**
     * @brief Swimming kinematic constraint phase shift (in radians between oscillators)
     *
     * @review_swimmer Change parameters stored for different systems
     *
     */
    double m_sys_spec_phase_shift{-1.0};
    /**
     * @brief Swimming kinematic constraint time-average spatial separation between a particle pair during oscillation
     *
     * @review_swimmer Change parameters stored for different systems
     *
     */
    double m_sys_spec_R_avg{-1.0};

    /* ANCHOR: particle parameters */
    /// (N x 1) REVIEW[epic=assumptions] 1) {0: constrained particle, 1: locater particle}.
    /// 2) locater particle listed first in order and denotes when to switch body index to next body.
    /// ex: (1, 0, 0, 1, 0, 0)
    Eigen::VectorXi m_particle_type_id;
    /// (N x 1) Group assembly (swimmer) number that each particle belongs to
    Eigen::VectorXi m_particle_group_id;

    /* ANCHOR: material parameters */
    /// mass density of fluid
    double m_fluid_density{-1};
    /// mass density of solid spheres
    double m_particle_density{-1};

    /* ANCHOR: potential parameters */
    /// WCA energy scale @f$ \epsilon_{\mathrm{WCA}} @f$
    double m_wca_epsilon{-1};
    /// WCA length scale @f$ \sigma_{\mathrm{WCA}} @f$
    double m_wca_sigma{-1};

    /* ANCHOR: degrees of freedom */
    /// = 3
    int m_num_spatial_dim{-1};
    /// = N
    int m_num_particles{-1};
    /// = M
    int m_num_bodies{-1};
    /// = 6M
    int m_num_DoF{-1};
    /// = M
    int m_num_constraints{-1};

    /* ANCHOR: integrator parameters */
    /// (dimensionless) integration @f$ \Delta t @f$
    double m_dt{-1.0};
    /// (dimensionless) final simulation time
    double m_tf{-1.0};
    /// (dimensionless) initial simulation time
    double m_t0{0.0};
    /// (dimensionless) current simulation time
    double m_t{0.0};
    /// simulation system characteristic timescale
    double m_tau{-1.0};
    /// current simulation timestep (iteration)
    int m_timestep{-1};
    /// how many simulation steps to output to GSD
    int m_num_steps_output{-1};

    /* ANCHOR: energetic parameters */
    double m_E_hydro_loc{0.0};
    double m_E_hydro_loc_int{0.0};
    double m_E_hydro_int{0.0};
    double m_E_hydro_simple{0.0};

    /* ANCHOR: Tensors set in constructor */
    // "identity" tensors
    /// (3 x 3) 2nd order identity tensor
    const Eigen::Matrix3d m_I3 = Eigen::Matrix3d::Identity(3, 3);
    /// (3 x 3) tensor version of `m_I3`
    const Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3>> m_tens_I3 = TensorCast(m_I3, 3, 3);

    // general-use tensors
    /// (3 x 3 x 3) (skew-symmetric) 3rd order identity tensor
    Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3, 3>> m_levi_cevita;
    /// (4 x 4 x 7) @f$ \nabla_{\xi_{\alpha}} \boldsymbol{E}^{T}{(\boldsymbol{\theta})} @f$
    Eigen::TensorFixedSize<double, Eigen::Sizes<4, 4, 7>> m_kappa;

    /* ANCHOR: kinematic vectors */
    /// (7M x 1) both linear and quaternion D.o.F. of body locaters
    Eigen::VectorXd m_positions_bodies;
    /// (7M x 1) both linear and quaternion D.o.F. of body locaters
    Eigen::VectorXd m_velocities_bodies;
    /// (7M x 1) both linear and quaternion D.o.F. of body locaters
    Eigen::VectorXd m_accelerations_bodies;

    /// (3N x 1) (linear) positions of all particles
    Eigen::VectorXd m_positions_particles;
    /// (7N x 1) (linear/angular) velocities of all particles
    Eigen::VectorXd m_velocities_particles;
    /// (7N x 1) (linear/angular) accelerations of all particles
    Eigen::VectorXd m_accelerations_particles;

    /// (3N x 1) orientations of all particles
    Eigen::VectorXd m_orientations_particles;
    /// (4N x 1) quaternions of all particles
    Eigen::VectorXd m_quaternions_particles;

    /// (3N x 1) (linear) *initial* (normalized) articulation positions of all particles
    Eigen::VectorXd m_positions_particles_articulation_init_norm;

    /// (3N x 1) (linear) articulation positions of all particles
    Eigen::VectorXd m_positions_particles_articulation;
    /// (7N x 1) (linear/angular) articulation velocities of all particles
    Eigen::VectorXd m_velocities_particles_articulation;
    /// (7N x 1) (linear/angular) articulation accelerations of all particles
    Eigen::VectorXd m_accelerations_particles_articulation;

    /* ANCHOR: rigid body motion tensors */
    /// (7M x 7N) @f$ \boldsymbol{\Sigma} @f$ rigid body motion connectivity tensor
    Eigen::MatrixXd m_rbm_conn;

    /// (7M x 7N) tensor version of `m_rbm_conn`
    Eigen::Tensor<double, 2> m_tens_rbm_conn;

    /// (7M x 7N x 7M) @f$ \nabla_{\xi} \boldsymbol{\sigma} @f$
    Eigen::Tensor<double, 3> m_tens_grad_rbm_conn;

    /// (7M x 3N) @f$ \boldsymbol{\chi} @f$ converts particle position D.o.F. to body position/quaternion D.o.F.
    Eigen::MatrixXd m_chi;

    /// (7M x 3N) tensor version of `m_chi`
    Eigen::Tensor<double, 2> m_tens_chi;

    /* ANCHOR: Udwadia constraint linear system */
    /// (number_constraints x 7M) linear operator defining relationship between constraints on @f$
    /// \ddot{\boldsymbol{\xi}} @f$.
    Eigen::MatrixXd m_Udwadia_A;

    /// (number of constraints x 1) Result of @f$ \mathbf{A} \, \ddot{\boldsymbol{\xi}} @f$
    Eigen::VectorXd m_Udwadia_b;
    /* !SECTION (Attributes) */

    /* SECTION: Setters and getters */
  public:
    /* ANCHOR: general attributes */
    bool
    imageSystem() const
    {
        return m_image_system;
    }
    void
    setImageSystem(bool image_system)
    {
        m_image_system = image_system;
    }

    // data i/o
    std::string
    inputGSDFile() const
    {
        return m_inputGSDFile;
    }

    std::string
    outputDir() const
    {
        return m_outputDir;
    }
    void
    setOutputDir(const std::string& outputDir)
    {
        m_outputDir = outputDir;
    }

    // GSD data
    std::shared_ptr<GSDUtil>
    gsdUtil() const
    {
        return m_gsdUtil;
    }

    std::shared_ptr<gsd_handle>
    handle() const
    {
        return m_handle;
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

    bool
    gSDParsed() const
    {
        return m_GSD_parsed;
    }

    /* ANCHOR: System specific data, change parameters stored for different systems */
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
    sysSpecRAvg() const
    {
        return m_sys_spec_R_avg;
    }
    void
    setSysSpecRAvg(double sys_spec_RAvg)
    {
        m_sys_spec_R_avg = sys_spec_RAvg;
    }

    /* ANCHOR: particle parameters */
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

    /* ANCHOR: material parameters */
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

    /* ANCHOR: potential parameters */
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

    /* ANCHOR: degrees of freedom */
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

    int
    numDoF()
    {
        return m_num_DoF;
    }

    /* ANCHOR: Tensors set in constructor */
    const Eigen::TensorFixedSize<double, Eigen::Sizes<3, 3>>&
    tensI3() const
    {
        return m_tens_I3;
    }

    const Eigen::Matrix3d&
    i3() const
    {
        return m_I3;
    }

    /* ANCHOR: integrator parameters */
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

    /* ANCHOR: energetic parameters */
    double
    eHydroLoc() const
    {
        return m_E_hydro_loc;
    }
    void
    setEHydroLoc(double E_hydro_loc)
    {
        m_E_hydro_loc = E_hydro_loc;
    }

    double
    eHydroLocInt() const
    {
        return m_E_hydro_loc_int;
    }
    void
    setEHydroLocInt(double E_hydro_loc_int)
    {
        m_E_hydro_loc_int = E_hydro_loc_int;
    }

    double
    eHydroInt() const
    {
        return m_E_hydro_int;
    }
    void
    setEHydroInt(double E_hydro_int)
    {
        m_E_hydro_int = E_hydro_int;
    }

    double
    eHydroSimple() const
    {
        return m_E_hydro_simple;
    }
    void
    setEHydroSimple(double E_hydro_simple)
    {
        m_E_hydro_simple = E_hydro_simple;
    }

    /* ANCHOR: kinematic vectors */
    const Eigen::MatrixXd&
    rbmConn() const
    {
        return m_rbm_conn;
    }
    // bodies
    const Eigen::VectorXd&
    positionsBodies() const
    {
        return m_positions_bodies;
    }
    void
    setPositionsBodies(const Eigen::VectorXd& positions_bodies)
    {
        m_positions_bodies = positions_bodies;
    }

    const Eigen::VectorXd&
    velocitiesBodies() const
    {
        return m_velocities_bodies;
    }
    void
    setVelocitiesBodies(const Eigen::VectorXd& velocities_bodies)
    {
        m_velocities_bodies = velocities_bodies;
    }

    const Eigen::VectorXd&
    accelerationsBodies() const
    {
        return m_accelerations_bodies;
    }
    void
    setAccelerationsBodies(const Eigen::VectorXd& accelerations_bodies)
    {
        m_accelerations_bodies = accelerations_bodies;
    }

    // particles
    const Eigen::VectorXd&
    quaternionsParticles() const
    {
        return m_quaternions_particles;
    }
    void
    setQuaternionsParticles(const Eigen::VectorXd& quaternions_particles)
    {
        m_quaternions_particles = quaternions_particles;
    }

    const Eigen::VectorXd&
    orientationsParticles() const
    {
        return m_quaternions_particles;
    }
    void
    setOrientationsParticles(const Eigen::VectorXd& orientations_particles)
    {
        m_quaternions_particles = orientations_particles;
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

    // particle articulations
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

    /* ANCHOR: rigid body motion tensors */
    const Eigen::Tensor<double, 2>&
    tensRbmConn() const
    {
        return m_tens_rbm_conn;
    }

    const Eigen::Tensor<double, 3>&
    tensGradRbmConn() const
    {
        return m_tens_grad_rbm_conn;
    }

    const Eigen::Tensor<double, 2>&
    tensChi() const
    {
        return m_tens_chi;
    }

    /* ANCHOR: Udwadia constraint linear system */
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
    /* !SECTION (Setters and getters) */
};

#endif
