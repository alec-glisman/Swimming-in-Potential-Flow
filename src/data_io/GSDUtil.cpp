//
// Created by Alec Glisman on 07/30/21
// Building off of GSDReader class in HOOMD-Blue
//

#include <GSDUtil.hpp>

GSDUtil::GSDUtil(std::shared_ptr<systemData> sys)
{
    // save classes
    m_system = sys;

    // Set member variables
    m_frame = 0;

    // Initialize logger
    m_logFile   = m_system->outputDir() + "/logs/" + m_logName + "-log.txt";
    auto logger = spdlog::basic_logger_mt(m_logName, m_logFile);
    spdlog::get(m_logName)->info("Initializing GSD reader");

    // Load GSD frame
    spdlog::get(m_logName)->info("Loading input GSD file");
    spdlog::get(m_logName)->info("GSD file path: {0}", m_system->inputGSDFile());
    auto return_val = gsd_open(m_system->handle().get(), m_system->inputGSDFile().c_str(), GSD_OPEN_READWRITE);

    m_system->setReturnVal(return_val);
    checkGSDReturn();

    // validate number of frames
    uint64_t nframes = gsd_get_nframes(m_system->handle().get());
    spdlog::get(m_logName)->info("{0} has {1} frames", m_system->inputGSDFile(),
                                 gsd_get_nframes(m_system->handle().get()));

    if (m_frame >= nframes)
    {
        spdlog::get(m_logName)->error("data.gsd_snapshot: Cannot read frame {0} {1} only has {2} frames", m_frame,
                                      m_system->inputGSDFile(), gsd_get_nframes(m_system->handle().get()));
        throw std::runtime_error("Error opening GSD file");
    }

    readHeader();
    readParameters();
    readParticles();

    readSystemSpecifics();

    spdlog::get(m_logName)->info("Constructor complete");
    spdlog::get(m_logName)->flush();
}

GSDUtil::GSDUtil(std::shared_ptr<systemData> sys, uint64_t frame) : GSDUtil(sys)
{
    m_frame = frame;

    // validate number of frames
    uint64_t nframes = gsd_get_nframes(m_system->handle().get());
    spdlog::get(m_logName)->info("{0} has {1} frames", m_system->inputGSDFile(),
                                 gsd_get_nframes(m_system->handle().get()));

    if (m_frame >= nframes)
    {
        spdlog::get(m_logName)->error("data.gsd_snapshot: Cannot read frame {0} {1} only has {2} frames", m_frame,
                                      m_system->inputGSDFile(), gsd_get_nframes(m_system->handle().get()));
        throw std::runtime_error("Error opening GSD file");
    }
}

GSDUtil::~GSDUtil()
{
    spdlog::get(m_logName)->info("GSDUtil destructor called");
    spdlog::get(m_logName)->flush();
}

void
GSDUtil::truncateGSD()
{
    spdlog::get(m_logName)->critical("truncating GSD file: {0}", m_system->inputGSDFile());
    auto return_val = gsd_truncate(m_system->handle().get());
    m_system->setReturnVal(return_val);
    checkGSDReturn();
}

void
GSDUtil::checkGSDReturn()
{
    if (m_system->returnVal() != 0)
    {
        spdlog::get(m_logName)->error("return_val = {0}", m_system->returnVal());
        spdlog::get(m_logName)->flush();
        throw std::runtime_error("Error parsing GSD file");
    }
    if (m_system->returnBool() == false)
    {
        spdlog::get(m_logName)->error("return_bool = {0}", m_system->returnBool());
        spdlog::get(m_logName)->flush();
        throw std::runtime_error("Error parsing GSD file");
    }
}

/* NOTE:
 */
bool
GSDUtil::readChunk(void* data, uint64_t frame, const char* name, size_t expected_size, unsigned int cur_n)
{
    const struct gsd_index_entry* entry = gsd_find_chunk(m_system->handle().get(), frame, name);

    if (entry == NULL && frame != 0)
    {
        entry = gsd_find_chunk(m_system->handle().get(), 0, name);
    }

    if (entry == NULL || (cur_n != 0 && entry->N != cur_n))
    {
        spdlog::get(m_logName)->warn("data.gsd_snapshot: chunk not found ");
        return false;
    }
    else
    {
        spdlog::get(m_logName)->info("data.gsd_snapshot: reading chunk {0}", name);
        size_t actual_size = entry->N * entry->M * gsd_sizeof_type((enum gsd_type)entry->type);

        if (actual_size != expected_size)
        {
            spdlog::get(m_logName)->error("data.gsd_snapshot: Expecting {0} bytes in {1} but found {2}", expected_size,
                                          name, actual_size);
        }

        auto return_val = gsd_read_chunk(m_system->handle().get(), data, entry);
        m_system->setReturnVal(return_val);
        checkGSDReturn();
    }

    return true;
}

void
GSDUtil::readHeader()
{
    spdlog::get(m_logName)->info("GSD parsing timestep");
    uint64_t timestep    = 0;
    auto     return_bool = readChunk(&timestep, m_frame, "log/configuration/step", 8);
    m_system->setReturnBool(return_bool);
    checkGSDReturn();
    m_system->setTimestep(int(timestep));
    spdlog::get(m_logName)->info("time step: {0}", timestep);
    assert(int(timestep) == m_system->timestep() && "time step not properly set");

    spdlog::get(m_logName)->info("GSD parsing dimensions");
    uint8_t dim = 0;
    return_bool = readChunk(&dim, m_frame, "log/configuration/dimensions", 1);
    m_system->setReturnBool(return_bool);
    checkGSDReturn();
    m_system->setNumSpatialDim(int(dim));
    spdlog::get(m_logName)->info("dim : {0}", dim);
    assert(int(dim) == m_system->numSpatialDim() && "number of dimensions not properly set");
    assert(int(dim) == 3 && "number of spatial dimensions must be 3 for  potential hydrodynamics");

    spdlog::get(m_logName)->info("GSD parsing number of particles");
    uint32_t N  = 0;
    return_bool = readChunk(&N, m_frame, "particles/N", 4);
    m_system->setReturnBool(return_bool);
    checkGSDReturn();
    m_system->setNumParticles(int(N));
    spdlog::get(m_logName)->info("Number of particles : {0}", N);
    if (N == 0)
    {
        spdlog::get(m_logFile)->error("data.gsd_snapshot: cannot read a file with 0 particles");
        throw std::runtime_error("Error reading GSD file");
    }
    assert(int(N) == m_system->numParticles() && "number of particles not properly set");
}

void
GSDUtil::readParameters()
{
    spdlog::get(m_logName)->info("GSD parsing dt");
    double dt{-1.0};
    auto   return_bool = readChunk(&dt, m_frame, "log/integrator/dt", 8);
    m_system->setReturnBool(return_bool);
    checkGSDReturn();
    m_system->setDt(dt);
    spdlog::get(m_logName)->info("dt : {0}", dt);
    assert(dt == m_system->dt() && "dt not properly set");

    spdlog::get(m_logName)->info("GSD parsing t");
    double t{-1.0};
    return_bool = readChunk(&t, m_frame, "log/integrator/t", 8);
    m_system->setReturnBool(return_bool);
    checkGSDReturn();
    m_system->setT(t);
    spdlog::get(m_logName)->info("t : {0}", t);
    assert(t == m_system->t() && "t not properly set");

    spdlog::get(m_logName)->info("GSD parsing tf");
    double tf{-1.0};
    return_bool = readChunk(&tf, m_frame, "log/integrator/tf", 8);
    m_system->setReturnBool(return_bool);
    checkGSDReturn();
    m_system->setTf(tf);
    spdlog::get(m_logName)->info("tf : {0}", tf);
    assert(tf == m_system->tf() && "tf not properly set");

    spdlog::get(m_logName)->info("GSD parsing tau");
    double tau{-1.0};
    return_bool = readChunk(&tau, m_frame, "log/integrator/tau", 8);
    m_system->setReturnBool(return_bool);
    checkGSDReturn();
    m_system->setTau(tau);
    spdlog::get(m_logName)->info("tau : {0}", tau);
    assert(tau == m_system->tau() && "tau not properly set");

    spdlog::get(m_logName)->info("GSD parsing num_steps_output");
    uint64_t num_steps_output{0};
    return_bool = readChunk(&num_steps_output, m_frame, "log/integrator/num_steps_output", 8);
    m_system->setReturnBool(return_bool);
    checkGSDReturn();
    m_system->setNumStepsOutput(int(num_steps_output));
    spdlog::get(m_logName)->info("num_steps_output : {0}", num_steps_output);
    assert(int(num_steps_output) == m_system->numStepsOutput() && "num_steps_output not properly set");

    spdlog::get(m_logName)->info("GSD parsing fluid_density");
    double fluid_density{-1.0};
    return_bool = readChunk(&fluid_density, m_frame, "log/material_parameters/fluid_density", 8);
    m_system->setReturnBool(return_bool);
    checkGSDReturn();
    m_system->setFluidDensity(fluid_density);
    spdlog::get(m_logName)->info("fluid_density : {0}", fluid_density);
    assert(fluid_density == m_system->fluidDensity() && "fluid_density not properly set");

    spdlog::get(m_logName)->info("GSD parsing particle_density");
    double particle_density{-1.0};
    return_bool = readChunk(&particle_density, m_frame, "log/material_parameters/particle_density", 8);
    m_system->setReturnBool(return_bool);
    checkGSDReturn();
    m_system->setParticleDensity(particle_density);
    spdlog::get(m_logName)->info("particle_density : {0}", particle_density);
    assert(particle_density == m_system->particleDensity() && "particle_density not properly set");

    spdlog::get(m_logName)->info("GSD parsing wca_epsilon");
    double wca_epsilon{-1.0};
    return_bool = readChunk(&wca_epsilon, m_frame, "log/wca/epsilon", 8);
    m_system->setReturnBool(return_bool);
    checkGSDReturn();
    m_system->setWcaEpsilon(wca_epsilon);
    spdlog::get(m_logName)->info("wca_epsilon : {0}", wca_epsilon);
    assert(wca_epsilon == m_system->wcaEpsilon() && "wca_epsilon not properly set");

    spdlog::get(m_logName)->info("GSD parsing wca_sigma");
    double wca_sigma{-1.0};
    return_bool = readChunk(&wca_sigma, m_frame, "log/wca/sigma", 8);
    m_system->setReturnBool(return_bool);
    checkGSDReturn();
    m_system->setWcaSigma(wca_sigma);
    spdlog::get(m_logName)->info("wca_sigma : {0}", wca_sigma);
    assert(wca_sigma == m_system->wcaSigma() && "wca_sigma not properly set");

    spdlog::get(m_logName)->info("GSD parsing typeid");
    uint32_t types[m_system->numParticles()];
    return_bool =
        readChunk(&types, m_frame, "particles/typeid", m_system->numParticles() * 4, m_system->numParticles());
    m_system->setReturnBool(return_bool);
    checkGSDReturn();
    Eigen::VectorXi type_id = Eigen::VectorXi::Zero(m_system->numParticles());
    for (int i = 0; i < m_system->numParticles(); i++)
    {
        type_id(i) = types[i];
        spdlog::get(m_logName)->info("Particle {0} typeid : {1}", i + 1, types[i]);
    }
    m_system->setParticleTypeId(type_id);

    spdlog::get(m_logName)->info("Calculating number of bodies.");
    int M = (type_id.array() == 1).count();
    m_system->setNumBodies(M);
    spdlog::get(m_logName)->info("num_bodies : {0}", M);
    assert(M == m_system->numBodies() && "num_bodies not properly set");
}

void
GSDUtil::readParticles()
{
    // initialize vectors
    int combined_body_tensor_len = 7 * m_system->numBodies();
    int orientational_tensor_len = 4 * m_system->numParticles();
    int spatial_tensor_len       = 3 * m_system->numParticles();

    Eigen::VectorXd n7_vec = Eigen::VectorXd::Zero(combined_body_tensor_len);
    Eigen::VectorXd n4_vec = Eigen::VectorXd::Zero(orientational_tensor_len);
    Eigen::VectorXd n3_vec = Eigen::VectorXd::Zero(spatial_tensor_len);

    m_system->setOrientationsParticles(n4_vec);
    m_system->setPositionsParticles(n3_vec);
    m_system->setVelocitiesParticles(n3_vec);
    m_system->setAccelerationsParticles(n3_vec);

    m_system->setPositionsBodies(n7_vec);
    m_system->setVelocitiesBodies(n7_vec);
    m_system->setAccelerationsBodies(n7_vec);

    // quaternions
    spdlog::get(m_logName)->info("GSD parsing orientationsParticles");
    double d_quat[4 * m_system->numParticles()];
    auto return_bool = readChunk(&d_quat, m_frame, "log/particles/double_orientation", m_system->numParticles() * 4 * 8,
                                 m_system->numParticles());
    m_system->setReturnBool(return_bool);
    checkGSDReturn();
    Eigen::VectorXd quaternions = Eigen::VectorXd::Zero(4 * m_system->numParticles());
    for (int i = 0; i < m_system->numParticles(); i++)
    {
        quaternions(4 * i)     = d_quat[4 * i];
        quaternions(4 * i + 1) = d_quat[4 * i + 1];
        quaternions(4 * i + 2) = d_quat[4 * i + 2];
        quaternions(4 * i + 3) = d_quat[4 * i + 3];
        spdlog::get(m_logName)->info("Particle {0} quaternion : [{1:03.3f}, {2:03.3f}, {3:03.3f}, {4:03.3f}]", i + 1,
                                     d_quat[4 * i], d_quat[4 * i + 1], d_quat[4 * i + 2], d_quat[4 * i + 3]);
    }
    m_system->setOrientationsParticles(quaternions);

    // positions
    spdlog::get(m_logName)->info("GSD parsing position");
    double d_pos[3 * m_system->numParticles()];
    return_bool = readChunk(&d_pos, m_frame, "log/particles/double_position", m_system->numParticles() * 3 * 8,
                            m_system->numParticles());
    m_system->setReturnBool(return_bool);
    checkGSDReturn();
    Eigen::VectorXd positions = Eigen::VectorXd::Zero(3 * m_system->numParticles());
    for (int i = 0; i < m_system->numParticles(); i++)
    {
        positions(3 * i)     = d_pos[3 * i];
        positions(3 * i + 1) = d_pos[3 * i + 1];
        positions(3 * i + 2) = d_pos[3 * i + 2];
        spdlog::get(m_logName)->info("Particle {0} position : [{1:03.3f}, {2:03.3f}, {3:03.3f}]", i + 1, d_pos[3 * i],
                                     d_pos[3 * i + 1], d_pos[3 * i + 2]);
    }
    m_system->setPositionsParticles(positions);

    // velocities
    spdlog::get(m_logName)->info("GSD parsing velocity");
    double d_vel[3 * m_system->numParticles()];
    return_bool = readChunk(&d_vel, m_frame, "log/particles/double_velocity", m_system->numParticles() * 3 * 8,
                            m_system->numParticles());
    m_system->setReturnBool(return_bool);
    checkGSDReturn();
    Eigen::VectorXd velocities = Eigen::VectorXd::Zero(3 * m_system->numParticles());
    for (int i = 0; i < m_system->numParticles(); i++)
    {
        velocities(3 * i)     = d_vel[3 * i];
        velocities(3 * i + 1) = d_vel[3 * i + 1];
        velocities(3 * i + 2) = d_vel[3 * i + 2];
        spdlog::get(m_logName)->info("Particle {0} velocity : [{1:03.3f}, {2:03.3f}, {3:03.3f}]", i + 1, d_vel[3 * i],
                                     d_vel[3 * i + 1], d_vel[3 * i + 2]);
    }
    m_system->setVelocitiesParticles(velocities);

    // accelerations
    spdlog::get(m_logName)->info("GSD parsing acceleration");
    double d_acc[3 * m_system->numParticles()];
    return_bool = readChunk(&d_acc, m_frame, "log/particles/double_moment_inertia", m_system->numParticles() * 3 * 8,
                            m_system->numParticles());
    m_system->setReturnBool(return_bool);
    checkGSDReturn();
    Eigen::VectorXd accelerations = Eigen::VectorXd::Zero(3 * m_system->numParticles());
    for (int i = 0; i < m_system->numParticles(); i++)
    {
        accelerations(3 * i)     = d_acc[3 * i];
        accelerations(3 * i + 1) = d_acc[3 * i + 1];
        accelerations(3 * i + 2) = d_acc[3 * i + 2];
        spdlog::get(m_logName)->info("Particle {0} acceleration : [{1:03.3f}, {2:03.3f}, {3:03.3f}]", i + 1,
                                     d_acc[3 * i], d_acc[3 * i + 1], d_acc[3 * i + 2]);
    }
    m_system->setAccelerationsParticles(accelerations);

    spdlog::get(m_logName)->info("Calculating body kinematics");
    spdlog::get(m_logName)->critical("Time derivatives of orientational components are assumed zero.");

    Eigen::VectorXd positions_bodies     = n7_vec;
    Eigen::VectorXd velocities_bodies    = n7_vec;
    Eigen::VectorXd accelerations_bodies = n7_vec;

    int body_count{0};

    for (int i = 0; i < m_system->numParticles(); i++)
    {
        // only fill for locater particles
        if (m_system->particleTypeId()(i) == 1)
        {
            const int body_id_7{7 * body_count};
            const int particle_id_3{3 * i};

            positions_bodies.segment<3>(body_id_7).noalias() = m_system->positionsParticles().segment<3>(particle_id_3);

            positions_bodies.segment<4>(body_id_7 + 3).noalias() =
                m_system->orientationsParticles().segment<4>(particle_id_3);

            velocities_bodies.segment<3>(body_id_7).noalias() =
                m_system->velocitiesParticles().segment<3>(particle_id_3);

            accelerations_bodies.segment<3>(body_id_7).noalias() =
                m_system->accelerationsParticles().segment<3>(particle_id_3);

            body_count++;
        }
    }
    assert(body_count == m_system->numBodies() && "Incorrect number of bodies filled");

    m_system->setPositionsBodies(positions_bodies);
    m_system->setVelocitiesBodies(velocities_bodies);
    m_system->setAccelerationsBodies(accelerations_bodies);
}

void
GSDUtil::readSystemSpecifics()
{
    spdlog::get(m_logName)->info("Running writeSystemSpecifics()");

    /* REVIEW[epic=Change]: set specific parameters */

    // oscillation velocity amplitude
    spdlog::get(m_logName)->info("GSD parsing U0");
    double U0{-1.0};
    auto   return_bool = readChunk(&U0, m_frame, "log/swimmer/U0", 8);
    m_system->setReturnBool(return_bool);
    checkGSDReturn();
    m_system->setSysSpecU0(U0);
    spdlog::get(m_logName)->info("U_0 : {0}", U0);
    assert(m_system->sysSpecU0() == U0 && "U0 not properly set");

    // oscillation frequency
    spdlog::get(m_logName)->info("GSD parsing omega");
    double omega{-1.0};
    return_bool = readChunk(&omega, m_frame, "log/swimmer/omega", 8);
    m_system->setReturnBool(return_bool);
    checkGSDReturn();
    m_system->setSysSpecOmega(omega);
    spdlog::get(m_logName)->info("omega : {0}", omega);
    assert(m_system->sysSpecOmega() == omega && "omega not properly set");

    // phase shift between oscillators
    spdlog::get(m_logName)->info("GSD parsing phaseShift");
    double phase_shift{-1.0};
    return_bool = readChunk(&phase_shift, m_frame, "log/swimmer/phase_shift", 8);
    m_system->setReturnBool(return_bool);
    checkGSDReturn();
    m_system->setSysSpecPhaseShift(phase_shift);
    spdlog::get(m_logName)->info("phase_shift : {0}", phase_shift);
    assert(m_system->sysSpecPhaseShift() == phase_shift && "phase_shift not properly set");

    // average separation
    spdlog::get(m_logName)->info("GSD parsing Ravg");
    double R_avg{-1.0};
    return_bool = readChunk(&R_avg, m_frame, "log/swimmer/R_avg", 8);
    m_system->setReturnBool(return_bool);
    checkGSDReturn();
    m_system->setSysSpecRAvg(R_avg);
    spdlog::get(m_logName)->info("R_avg : {0}", R_avg);
    assert(m_system->sysSpecRAvg() == R_avg && "R_avg not properly set");
}

void
GSDUtil::writeFrame()
{
    spdlog::get(m_logName)->info("GSD writing frame");
    spdlog::get(m_logName)->info("time step: {0}", m_system->timestep());
    spdlog::get(m_logName)->info("time: {0}", m_system->t());

    writeHeader();
    writeParameters();
    writeParticles();

    spdlog::get(m_logName)->info("GSD ending frame");
    auto return_val = gsd_end_frame(m_system->handle().get());
    m_system->setReturnVal(return_val);
    checkGSDReturn();
}

void
GSDUtil::writeHeader()
{
    spdlog::get(m_logName)->info("GSD writing log/configuration/timestep");
    uint64_t step = m_system->timestep();
    auto     return_val =
        gsd_write_chunk(m_system->handle().get(), "log/configuration/step", GSD_TYPE_UINT64, 1, 1, 0, (void*)&step);
    m_system->setReturnVal(return_val);
    checkGSDReturn();

    if (gsd_get_nframes(m_system->handle().get()) == 0)
    {
        spdlog::get(m_logName)->info("GSD writing configuration/dimensions");
        uint8_t dimensions = 3;
        return_val = gsd_write_chunk(m_system->handle().get(), "configuration/dimensions", GSD_TYPE_UINT8, 1, 1, 0,
                                     (void*)&dimensions);
        m_system->setReturnVal(return_val);
        checkGSDReturn();
    }

    spdlog::get(m_logName)->info("GSD writing particles/N");
    uint32_t N = m_system->numParticles();
    return_val = gsd_write_chunk(m_system->handle().get(), "particles/N", GSD_TYPE_UINT32, 1, 1, 0, (void*)&N);
    m_system->setReturnVal(return_val);
    checkGSDReturn();
}

void
GSDUtil::writeParameters()
{
    spdlog::get(m_logName)->info("GSD writing log/integrator/dt");
    double dt = m_system->dt();
    auto   return_val =
        gsd_write_chunk(m_system->handle().get(), "log/integrator/dt", GSD_TYPE_DOUBLE, 1, 1, 0, (void*)&dt);
    m_system->setReturnVal(return_val);
    checkGSDReturn();

    spdlog::get(m_logName)->info("GSD writing log/integrator/t");
    double time = m_system->t();
    return_val  = gsd_write_chunk(m_system->handle().get(), "log/integrator/t", GSD_TYPE_DOUBLE, 1, 1, 0, (void*)&time);
    m_system->setReturnVal(return_val);
    checkGSDReturn();
}

void
GSDUtil::writeParticles()
{
    uint32_t N       = m_system->numParticles();
    uint32_t num_dim = 3 * N;
    int      return_val;
    uint64_t nframes = gsd_get_nframes(m_system->handle().get());
    spdlog::get(m_logName)->critical("vectors are assumed to have 3 spatial DoF");

    /* ANCHOR: Write kinematics using standard data structures, which are floats */
    std::vector<float> quat(uint64_t(N) * uint64_t(num_dim + N));
    quat.reserve(1); // make sure we allocate
    for (uint32_t i = 0; i < N; i++)
    {
        quat[4 * i]     = m_system->orientationsParticles()(4 * i);
        quat[4 * i + 1] = m_system->orientationsParticles()(4 * i + 1);
        quat[4 * i + 2] = m_system->orientationsParticles()(4 * i + 2);
        quat[4 * i + 3] = m_system->orientationsParticles()(4 * i + 3);
    }
    spdlog::get(m_logName)->info("GSD writing particles/orientation");
    return_val =
        gsd_write_chunk(m_system->handle().get(), "particles/orientation", GSD_TYPE_FLOAT, N, 4, 0, (void*)&quat[0]);
    m_system->setReturnVal(return_val);
    checkGSDReturn();

    std::vector<float> pos(uint64_t(N) * uint64_t(num_dim));
    pos.reserve(1); // make sure we allocate
    for (uint32_t i = 0; i < N; i++)
    {
        pos[3 * i]     = m_system->positionsParticles()(3 * i);
        pos[3 * i + 1] = m_system->positionsParticles()(3 * i + 1);
        pos[3 * i + 2] = m_system->positionsParticles()(3 * i + 2);
    }
    spdlog::get(m_logName)->info("GSD writing particles/position");
    return_val =
        gsd_write_chunk(m_system->handle().get(), "particles/position", GSD_TYPE_FLOAT, N, 3, 0, (void*)&pos[0]);
    m_system->setReturnVal(return_val);
    checkGSDReturn();

    std::vector<float> vel(uint64_t(N) * uint64_t(num_dim));
    vel.reserve(1); // make sure we allocate
    for (uint32_t i = 0; i < N; i++)
    {
        vel[3 * i]     = m_system->velocitiesParticles()(3 * i);
        vel[3 * i + 1] = m_system->velocitiesParticles()(3 * i + 1);
        vel[3 * i + 2] = m_system->velocitiesParticles()(3 * i + 2);
    }
    spdlog::get(m_logName)->info("GSD writing particles/velocity");
    return_val =
        gsd_write_chunk(m_system->handle().get(), "particles/velocity", GSD_TYPE_FLOAT, N, 3, 0, (void*)&vel[0]);
    m_system->setReturnVal(return_val);
    checkGSDReturn();

    std::vector<float> acc(uint64_t(N) * uint64_t(num_dim));
    acc.reserve(1); // make sure we allocate
    for (uint32_t i = 0; i < N; i++)
    {
        acc[3 * i]     = m_system->accelerationsParticles()(3 * i);
        acc[3 * i + 1] = m_system->accelerationsParticles()(3 * i + 1);
        acc[3 * i + 2] = m_system->accelerationsParticles()(3 * i + 2);
    }
    spdlog::get(m_logName)->info("GSD writing particles/moment_inertia");
    return_val =
        gsd_write_chunk(m_system->handle().get(), "particles/moment_inertia", GSD_TYPE_FLOAT, N, 3, 0, (void*)&acc[0]);
    m_system->setReturnVal(return_val);
    checkGSDReturn();

    /* ANCHOR: Write kinematics as doubles for higher precision */
    std::vector<double> d_quat(uint64_t(N) * uint64_t(num_dim + N));
    d_quat.reserve(1); // make sure we allocate
    for (uint32_t i = 0; i < N; i++)
    {
        d_quat[4 * i]     = m_system->orientationsParticles()(4 * i);
        d_quat[4 * i + 1] = m_system->orientationsParticles()(4 * i + 1);
        d_quat[4 * i + 2] = m_system->orientationsParticles()(4 * i + 2);
        d_quat[4 * i + 3] = m_system->orientationsParticles()(4 * i + 3);
    }
    spdlog::get(m_logName)->info("GSD writing log/particles/double_orientation");
    return_val = gsd_write_chunk(m_system->handle().get(), "log/particles/double_orientation", GSD_TYPE_DOUBLE, N, 4, 0,
                                 (void*)&d_quat[0]);
    m_system->setReturnVal(return_val);
    checkGSDReturn();

    std::vector<double> d_pos(uint64_t(N) * uint64_t(num_dim));
    d_pos.reserve(1); // make sure we allocate
    for (uint32_t i = 0; i < N; i++)
    {
        d_pos[3 * i]     = m_system->positionsParticles()(3 * i);
        d_pos[3 * i + 1] = m_system->positionsParticles()(3 * i + 1);
        d_pos[3 * i + 2] = m_system->positionsParticles()(3 * i + 2);
    }
    spdlog::get(m_logName)->info("GSD writing log/particles/double_position");
    return_val = gsd_write_chunk(m_system->handle().get(), "log/particles/double_position", GSD_TYPE_DOUBLE, N, 3, 0,
                                 (void*)&d_pos[0]);
    m_system->setReturnVal(return_val);
    checkGSDReturn();

    std::vector<double> d_vel(uint64_t(N) * uint64_t(num_dim));
    d_vel.reserve(1); // make sure we allocate
    for (uint32_t i = 0; i < N; i++)
    {
        d_vel[3 * i]     = m_system->velocitiesParticles()(3 * i);
        d_vel[3 * i + 1] = m_system->velocitiesParticles()(3 * i + 1);
        d_vel[3 * i + 2] = m_system->velocitiesParticles()(3 * i + 2);
    }
    spdlog::get(m_logName)->info("GSD wri ting log/particles/double_velocity");
    return_val = gsd_write_chunk(m_system->handle().get(), "log/particles/double_velocity", GSD_TYPE_DOUBLE, N, 3, 0,
                                 (void*)&d_vel[0]);
    m_system->setReturnVal(return_val);
    checkGSDReturn();

    std::vector<double> d_acc(uint64_t(N) * uint64_t(num_dim));
    d_acc.reserve(1); // make sure we allocate
    for (uint32_t i = 0; i < N; i++)
    {
        d_acc[3 * i]     = m_system->accelerationsParticles()(3 * i);
        d_acc[3 * i + 1] = m_system->accelerationsParticles()(3 * i + 1);
        d_acc[3 * i + 2] = m_system->accelerationsParticles()(3 * i + 2);
    }
    spdlog::get(m_logName)->info("GSD writing log/particles/double_moment_inertia");
    return_val = gsd_write_chunk(m_system->handle().get(), "log/particles/double_moment_inertia", GSD_TYPE_DOUBLE, N, 3,
                                 0, (void*)&d_acc[0]);
    m_system->setReturnVal(return_val);
    checkGSDReturn();
}