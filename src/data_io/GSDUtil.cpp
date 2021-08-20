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
    auto return_val =
        gsd_open(m_system->handle().get(), m_system->inputGSDFile().c_str(), GSD_OPEN_READWRITE);

    m_system->setReturnVal(return_val);
    m_system->check_gsd_return();

    // validate number of frames
    uint64_t nframes = gsd_get_nframes(m_system->handle().get());
    spdlog::get(m_logName)->info("{0} has {1} frames", m_system->inputGSDFile(),
                                 gsd_get_nframes(m_system->handle().get()));

    if (m_frame >= nframes)
    {
        spdlog::get(m_logName)->error(
            "data.gsd_snapshot: Cannot read frame {0} {1} only has {2} frames", m_frame,
            m_system->inputGSDFile(), gsd_get_nframes(m_system->handle().get()));
        throw std::runtime_error("Error opening GSD file");
    }

    readHeader();
    readParameters();
    readParticles();
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
        spdlog::get(m_logName)->error(
            "data.gsd_snapshot: Cannot read frame {0} {1} only has {2} frames", m_frame,
            m_system->inputGSDFile(), gsd_get_nframes(m_system->handle().get()));
        throw std::runtime_error("Error opening GSD file");
    }
}

GSDUtil::~GSDUtil()
{
    spdlog::get(m_logName)->info("GSDUtil destructor called");
}

void
GSDUtil::truncateGSD()
{
    spdlog::get(m_logName)->info("truncating GSD file: {0}", m_system->inputGSDFile());
    auto return_val = gsd_truncate(m_system->handle().get());
    m_system->setReturnVal(return_val);
    m_system->check_gsd_return();
}

/* NOTE:
 * Expected size is in units of bytes
 *     u_int8: 1
 *     float (np.single, np.float32): 4
 *     double (np.double, np.float64): 8
 */
bool
GSDUtil::readChunk(void* data, uint64_t frame, const char* name, size_t expected_size,
                   unsigned int cur_n)
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
            spdlog::get(m_logName)->error(
                "data.gsd_snapshot: Expecting {0} bytes in {1} but found {2}", expected_size, name,
                actual_size);
        }

        auto return_val = gsd_read_chunk(m_system->handle().get(), data, entry);
        m_system->setReturnVal(return_val);
        m_system->check_gsd_return();
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
    m_system->check_gsd_return();
    m_system->setTimestep(int(timestep));
    spdlog::get(m_logName)->info("time step: {0}", timestep);
    assert(int(timestep) == m_system->timestep() && "time step not properly set");

    spdlog::get(m_logName)->info("GSD parsing dimensions");
    uint8_t dim = 0;
    return_bool = readChunk(&dim, m_frame, "log/configuration/dimensions", 1);
    m_system->setReturnBool(return_bool);
    m_system->check_gsd_return();
    m_system->setNumDim(int(dim));
    spdlog::get(m_logName)->info("dim : {0}", dim);
    assert(int(dim) == m_system->numDim() && "number of dimensions not properly set");

    spdlog::get(m_logName)->info("GSD parsing number of particles");
    uint32_t N  = 0;
    return_bool = readChunk(&N, m_frame, "particles/N", 4);
    m_system->setReturnBool(return_bool);
    m_system->check_gsd_return();
    m_system->setNumParticles(int(N));
    spdlog::get(m_logName)->info("Number of particles : {0}", N);
    if (N == 0)
    {
        spdlog::get(m_logFile)->error("data.gsd_snapshot: cannot read a file with 0 particles");
        throw std::runtime_error("Error reading GSD file");
    }
    assert(int(N) == m_system->numParticles() && "number of particles not properly set");

    // Initialize tensors
    m_system->resizeTensors();
}

void
GSDUtil::readParameters()
{
    spdlog::get(m_logName)->info("GSD parsing dt");
    double dt{-1.0};
    auto   return_bool = readChunk(&dt, m_frame, "log/integrator/dt", 8);
    m_system->setReturnBool(return_bool);
    m_system->check_gsd_return();
    m_system->setDt(dt);
    spdlog::get(m_logName)->info("dt : {0}", dt);
    assert(dt == m_system->dt() && "dt not properly set");

    spdlog::get(m_logName)->info("GSD parsing t");
    double t{-1.0};
    return_bool = readChunk(&t, m_frame, "log/integrator/t", 8);
    m_system->setReturnBool(return_bool);
    m_system->check_gsd_return();
    m_system->setT(t);
    spdlog::get(m_logName)->info("t : {0}", t);
    assert(t == m_system->t() && "t not properly set");

    spdlog::get(m_logName)->info("GSD parsing tf");
    double tf{-1.0};
    return_bool = readChunk(&tf, m_frame, "log/integrator/tf", 8);
    m_system->setReturnBool(return_bool);
    m_system->check_gsd_return();
    m_system->setTf(tf);
    spdlog::get(m_logName)->info("tf : {0}", tf);
    assert(tf == m_system->tf() && "tf not properly set");

    spdlog::get(m_logName)->info("GSD parsing tau");
    double tau{-1.0};
    return_bool = readChunk(&tau, m_frame, "log/integrator/tau", 8);
    m_system->setReturnBool(return_bool);
    m_system->check_gsd_return();
    m_system->setTau(tau);
    spdlog::get(m_logName)->info("tau : {0}", tau);
    assert(tau == m_system->tau() && "tau not properly set");

    spdlog::get(m_logName)->info("GSD parsing num_steps_output");
    uint64_t num_steps_output{0};
    return_bool = readChunk(&num_steps_output, m_frame, "log/integrator/num_steps_output", 8);
    m_system->setReturnBool(return_bool);
    m_system->check_gsd_return();
    m_system->setNumStepsOutput(int(num_steps_output));
    spdlog::get(m_logName)->info("num_steps_output : {0}", num_steps_output);
    assert(int(num_steps_output) == m_system->numStepsOutput() &&
           "num_steps_output not properly set");

    spdlog::get(m_logName)->info("GSD parsing fluid_density");
    double fluid_density{-1.0};
    return_bool = readChunk(&fluid_density, m_frame, "log/material_parameters/fluid_density", 8);
    m_system->setReturnBool(return_bool);
    m_system->check_gsd_return();
    m_system->setFluidDensity(fluid_density);
    spdlog::get(m_logName)->info("fluid_density : {0}", fluid_density);
    assert(fluid_density == m_system->fluidDensity() && "fluid_density not properly set");

    spdlog::get(m_logName)->info("GSD parsing particle_density");
    double particle_density{-1.0};
    return_bool =
        readChunk(&particle_density, m_frame, "log/material_parameters/particle_density", 8);
    m_system->setReturnBool(return_bool);
    m_system->check_gsd_return();
    m_system->setParticleDensity(particle_density);
    spdlog::get(m_logName)->info("particle_density : {0}", particle_density);
    assert(particle_density == m_system->particleDensity() && "particle_density not properly set");

    spdlog::get(m_logName)->info("GSD parsing wca_epsilon");
    double wca_epsilon{-1.0};
    return_bool = readChunk(&wca_epsilon, m_frame, "log/wca/epsilon", 8);
    m_system->setReturnBool(return_bool);
    m_system->check_gsd_return();
    m_system->setWcaEpsilon(wca_epsilon);
    spdlog::get(m_logName)->info("wca_epsilon : {0}", wca_epsilon);
    assert(wca_epsilon == m_system->wcaEpsilon() && "wca_epsilon not properly set");

    spdlog::get(m_logName)->info("GSD parsing wca_sigma");
    double wca_sigma{-1.0};
    return_bool = readChunk(&wca_sigma, m_frame, "log/wca/sigma", 8);
    m_system->setReturnBool(return_bool);
    m_system->check_gsd_return();
    m_system->setWcaSigma(wca_sigma);
    spdlog::get(m_logName)->info("wca_sigma : {0}", wca_sigma);
    assert(wca_sigma == m_system->wcaSigma() && "wca_sigma not properly set");
}

void
GSDUtil::readParticles()
{
    // positions
    spdlog::get(m_logName)->info("GSD parsing position");
    double d_pos[m_system->numDim() * m_system->numParticles()];
    auto   return_bool = readChunk(&d_pos, m_frame, "log/particles/double_position",
                                 m_system->numParticles() * 3 * 8, m_system->numParticles());
    m_system->setReturnBool(return_bool);
    m_system->check_gsd_return();
    Eigen::VectorXd positions = Eigen::VectorXd::Zero(3 * m_system->numParticles());
    for (int i = 0; i < m_system->numParticles(); i++)
    {
        positions(3 * i)     = d_pos[3 * i];
        positions(3 * i + 1) = d_pos[3 * i + 1];
        positions(3 * i + 2) = d_pos[3 * i + 2];
        spdlog::get(m_logName)->info("Particle {0} position : [{1:03.3f}, {2:03.3f}, {3:03.3f}]",
                                     i + 1, d_pos[m_system->numDim() * i],
                                     d_pos[m_system->numDim() * i + 1],
                                     d_pos[m_system->numDim() * i + 2]);
    }
    m_system->setPositions(positions);

    // velocities
    spdlog::get(m_logName)->info("GSD parsing velocity");
    double d_vel[m_system->numDim() * m_system->numParticles()];
    return_bool = readChunk(&d_vel, m_frame, "log/particles/double_velocity",
                            m_system->numParticles() * 3 * 8, m_system->numParticles());
    m_system->setReturnBool(return_bool);
    m_system->check_gsd_return();
    Eigen::VectorXd velocities = Eigen::VectorXd::Zero(3 * m_system->numParticles());
    for (int i = 0; i < m_system->numParticles(); i++)
    {
        velocities(3 * i)     = d_vel[3 * i];
        velocities(3 * i + 1) = d_vel[3 * i + 1];
        velocities(3 * i + 2) = d_vel[3 * i + 2];
        spdlog::get(m_logName)->info("Particle {0} velocity : [{1:03.3f}, {2:03.3f}, {3:03.3f}]",
                                     i + 1, d_vel[m_system->numDim() * i],
                                     d_vel[m_system->numDim() * i + 1],
                                     d_vel[m_system->numDim() * i + 2]);
    }
    m_system->setVelocities(velocities);

    // accelerations
    spdlog::get(m_logName)->info("GSD parsing acceleration");
    double d_acc[m_system->numDim() * m_system->numParticles()];
    return_bool = readChunk(&d_acc, m_frame, "log/particles/double_moment_inertia",
                            m_system->numParticles() * 3 * 8, m_system->numParticles());
    m_system->setReturnBool(return_bool);
    m_system->check_gsd_return();
    Eigen::VectorXd accelerations = Eigen::VectorXd::Zero(3 * m_system->numParticles());
    for (int i = 0; i < m_system->numParticles(); i++)
    {
        accelerations(3 * i)     = d_acc[3 * i];
        accelerations(3 * i + 1) = d_acc[3 * i + 1];
        accelerations(3 * i + 2) = d_acc[3 * i + 2];
        spdlog::get(m_logName)->info(
            "Particle {0} acceleration : [{1:03.3f}, {2:03.3f}, {3:03.3f}]", i + 1,
            d_acc[m_system->numDim() * i], d_acc[m_system->numDim() * i + 1],
            d_acc[m_system->numDim() * i + 2]);
    }
    m_system->setAccelerations(accelerations);
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
    m_system->check_gsd_return();
}

void
GSDUtil::writeHeader()
{
    spdlog::get(m_logName)->info("GSD writing log/configuration/timestep");
    uint64_t step       = m_system->timestep();
    auto     return_val = gsd_write_chunk(m_system->handle().get(), "log/configuration/step",
                                      GSD_TYPE_UINT64, 1, 1, 0, (void*)&step);
    m_system->setReturnVal(return_val);
    m_system->check_gsd_return();

    if (gsd_get_nframes(m_system->handle().get()) == 0)
    {
        spdlog::get(m_logName)->info("GSD writing configuration/dimensions");
        uint8_t dimensions = m_system->numDim();
        return_val         = gsd_write_chunk(m_system->handle().get(), "configuration/dimensions",
                                     GSD_TYPE_UINT8, 1, 1, 0, (void*)&dimensions);
        m_system->setReturnVal(return_val);
        m_system->check_gsd_return();
    }

    spdlog::get(m_logName)->info("GSD writing particles/N");
    uint32_t N = m_system->numParticles();
    return_val = gsd_write_chunk(m_system->handle().get(), "particles/N", GSD_TYPE_UINT32, 1, 1, 0,
                                 (void*)&N);
    m_system->setReturnVal(return_val);
    m_system->check_gsd_return();
}

void
GSDUtil::writeParameters()
{
    spdlog::get(m_logName)->info("GSD writing log/integrator/dt");
    double dt         = m_system->dt();
    auto   return_val = gsd_write_chunk(m_system->handle().get(), "log/integrator/dt",
                                      GSD_TYPE_DOUBLE, 1, 1, 0, (void*)&dt);
    m_system->setReturnVal(return_val);
    m_system->check_gsd_return();

    spdlog::get(m_logName)->info("GSD writing log/integrator/t");
    double time = m_system->t();
    return_val  = gsd_write_chunk(m_system->handle().get(), "log/integrator/t", GSD_TYPE_DOUBLE, 1,
                                 1, 0, (void*)&time);
    m_system->setReturnVal(return_val);
    m_system->check_gsd_return();
}

void
GSDUtil::writeParticles()
{
    uint32_t N       = m_system->numParticles();
    uint32_t num_dim = m_system->numDim();
    int      return_val;
    uint64_t nframes = gsd_get_nframes(m_system->handle().get());
    spdlog::get(m_logName)->debug("vectors are assumed to have 3 spatial DoF");

    /* ANCHOR: Write kinematics using standard data structures, which are floats */
    std::vector<float> pos(uint64_t(N) * uint64_t(num_dim));
    pos.reserve(1); //! make sure we allocate
    for (uint32_t i = 0; i < N; i++)
    {
        pos[3 * i]     = m_system->positions()(3 * i);
        pos[3 * i + 1] = m_system->positions()(3 * i + 1);
        pos[3 * i + 2] = m_system->positions()(3 * i + 2);
    }
    spdlog::get(m_logName)->info("GSD writing particles/position");
    return_val = gsd_write_chunk(m_system->handle().get(), "particles/position", GSD_TYPE_FLOAT, N,
                                 3, 0, (void*)&pos[0]);
    m_system->setReturnVal(return_val);
    m_system->check_gsd_return();

    std::vector<float> vel(uint64_t(N) * uint64_t(num_dim));
    vel.reserve(1); //! make sure we allocate
    for (uint32_t i = 0; i < N; i++)
    {
        vel[3 * i]     = m_system->velocities()(3 * i);
        vel[3 * i + 1] = m_system->velocities()(3 * i + 1);
        vel[3 * i + 2] = m_system->velocities()(3 * i + 2);
    }
    spdlog::get(m_logName)->info("GSD writing particles/velocity");
    return_val = gsd_write_chunk(m_system->handle().get(), "particles/velocity", GSD_TYPE_FLOAT, N,
                                 3, 0, (void*)&vel[0]);
    m_system->setReturnVal(return_val);
    m_system->check_gsd_return();

    std::vector<float> acc(uint64_t(N) * uint64_t(num_dim));
    acc.reserve(1); //! make sure we allocate
    for (uint32_t i = 0; i < N; i++)
    {
        acc[3 * i]     = m_system->accelerations()(3 * i);
        acc[3 * i + 1] = m_system->accelerations()(3 * i + 1);
        acc[3 * i + 2] = m_system->accelerations()(3 * i + 2);
    }
    spdlog::get(m_logName)->info("GSD writing particles/moment_inertia");
    return_val = gsd_write_chunk(m_system->handle().get(), "particles/moment_inertia",
                                 GSD_TYPE_FLOAT, N, 3, 0, (void*)&acc[0]);
    m_system->setReturnVal(return_val);
    m_system->check_gsd_return();

    /* ANCHOR: Write kinematics as doubles for higher precision */
    std::vector<double> d_pos(uint64_t(N) * uint64_t(num_dim));
    d_pos.reserve(1); //! make sure we allocate
    for (uint32_t i = 0; i < N; i++)
    {
        d_pos[3 * i]     = m_system->positions()(3 * i);
        d_pos[3 * i + 1] = m_system->positions()(3 * i + 1);
        d_pos[3 * i + 2] = m_system->positions()(3 * i + 2);
    }
    spdlog::get(m_logName)->info("GSD writing log/particles/double_position");
    return_val = gsd_write_chunk(m_system->handle().get(), "log/particles/double_position",
                                 GSD_TYPE_DOUBLE, N, 3, 0, (void*)&d_pos[0]);
    m_system->setReturnVal(return_val);
    m_system->check_gsd_return();

    std::vector<double> d_vel(uint64_t(N) * uint64_t(num_dim));
    d_vel.reserve(1); //! make sure we allocate
    for (uint32_t i = 0; i < N; i++)
    {
        d_vel[3 * i]     = m_system->velocities()(3 * i);
        d_vel[3 * i + 1] = m_system->velocities()(3 * i + 1);
        d_vel[3 * i + 2] = m_system->velocities()(3 * i + 2);
    }
    spdlog::get(m_logName)->info("GSD wri ting log/particles/double_velocity");
    return_val = gsd_write_chunk(m_system->handle().get(), "log/particles/double_velocity",
                                 GSD_TYPE_DOUBLE, N, 3, 0, (void*)&d_vel[0]);
    m_system->setReturnVal(return_val);
    m_system->check_gsd_return();

    std::vector<double> d_acc(uint64_t(N) * uint64_t(num_dim));
    d_acc.reserve(1); //! make sure we allocate
    for (uint32_t i = 0; i < N; i++)
    {
        d_acc[3 * i]     = m_system->accelerations()(3 * i);
        d_acc[3 * i + 1] = m_system->accelerations()(3 * i + 1);
        d_acc[3 * i + 2] = m_system->accelerations()(3 * i + 2);
    }
    spdlog::get(m_logName)->info("GSD writing log/particles/double_moment_inertia");
    return_val = gsd_write_chunk(m_system->handle().get(), "log/particles/double_moment_inertia",
                                 GSD_TYPE_DOUBLE, N, 3, 0, (void*)&d_acc[0]);
    m_system->setReturnVal(return_val);
    m_system->check_gsd_return();
}