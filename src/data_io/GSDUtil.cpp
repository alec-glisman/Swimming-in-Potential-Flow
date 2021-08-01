//
// Created by Alec Glisman on 07/30/21
// Building off of GSDReader class in HOOMD-Blue
//

#include <GSDUtil.hpp>

GSDUtil::GSDUtil(systemData& sys)
{
    // save classes
    system = &sys;

    // Set member variables
    m_frame = 0;

    // Initialize logger
    m_logFile   = system->outputDir() + "/logs/" + m_logName + "-log.txt";
    auto logger = spdlog::basic_logger_mt(m_logName, m_logFile);
    spdlog::get(m_logName)->info("Initializing GSD reader");

    // Load GSD frame
    spdlog::get(m_logName)->info("Loading input GSD file");
    spdlog::get(m_logName)->info("GSD file path: {0}", system->inputGSDFile());
    auto return_val =
        gsd_open(system->handle().get(), system->inputGSDFile().c_str(), GSD_OPEN_READWRITE);

    system->setReturnVal(return_val);
    system->check_gsd_return();

    // validate number of frames
    uint64_t nframes = gsd_get_nframes(system->handle().get());
    spdlog::get(m_logName)->info("{0} has {1} frames", system->inputGSDFile(),
                                 gsd_get_nframes(system->handle().get()));

    if (m_frame >= nframes)
    {
        spdlog::get(m_logName)->error(
            "data.gsd_snapshot: Cannot read frame {0} {1} only has {2} frames", m_frame,
            system->inputGSDFile(), gsd_get_nframes(system->handle().get()));
        throw std::runtime_error("Error opening GSD file");
    }

    readHeader();
    readParameters();
    readParticles();
}

GSDUtil::GSDUtil(systemData& sys, uint64_t frame) : GSDUtil(sys)
{
    m_frame = frame;

    // validate number of frames
    uint64_t nframes = gsd_get_nframes(system->handle().get());
    spdlog::get(m_logName)->info("{0} has {1} frames", system->inputGSDFile(),
                                 gsd_get_nframes(system->handle().get()));

    if (m_frame >= nframes)
    {
        spdlog::get(m_logName)->error(
            "data.gsd_snapshot: Cannot read frame {0} {1} only has {2} frames", m_frame,
            system->inputGSDFile(), gsd_get_nframes(system->handle().get()));
        throw std::runtime_error("Error opening GSD file");
    }
}

GSDUtil::~GSDUtil()
{
    spdlog::get(m_logName)->info("GSDUtil destructor called");
}

/* NOTES
 * Expected size is in units of bytes
 *     u_int8: 1
 *     float (np.single, np.float32): 4
 *     double (np.double, np.float64): 8
 */
bool
GSDUtil::readChunk(void* data, uint64_t frame, const char* name, size_t expected_size,
                   unsigned int cur_n)
{
    const struct gsd_index_entry* entry = gsd_find_chunk(system->handle().get(), frame, name);

    if (entry == NULL && frame != 0)
    {
        entry = gsd_find_chunk(system->handle().get(), 0, name);
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

        auto return_val = gsd_read_chunk(system->handle().get(), data, entry);
        system->setReturnVal(return_val);
        system->check_gsd_return();
    }

    return true;
}

void
GSDUtil::readHeader()
{
    spdlog::get(m_logName)->info("GSD parsing timestep");
    uint64_t timestep    = 0;
    auto     return_bool = readChunk(&timestep, m_frame, "log/configuration/step", 8);
    system->setReturnBool(return_bool);
    system->check_gsd_return();
    system->setTimestep(int(timestep));
    spdlog::get(m_logName)->info("time step: {0}", timestep);
    assert(int(timestep) == system->timestep() && "time step not properly set");

    spdlog::get(m_logName)->info("GSD parsing dimensions");
    uint8_t dim = 0;
    return_bool = readChunk(&dim, m_frame, "log/configuration/dimensions", 1);
    system->setReturnBool(return_bool);
    system->check_gsd_return();
    system->setNumDim(int(dim));
    spdlog::get(m_logName)->info("dim : {0}", dim);
    assert(int(dim) == system->numDim() && "number of dimensions not properly set");

    spdlog::get(m_logName)->info("GSD parsing number of particles");
    uint32_t N  = 0;
    return_bool = readChunk(&N, m_frame, "particles/N", 4);
    system->setReturnBool(return_bool);
    system->check_gsd_return();
    system->setNumParticles(int(N));
    spdlog::get(m_logName)->info("Number of particles : {0}", N);
    if (N == 0)
    {
        spdlog::get(m_logFile)->error("data.gsd_snapshot: cannot read a file with 0 particles");
        throw std::runtime_error("Error reading GSD file");
    }
    assert(int(N) == system->numParticles() && "number of particles not properly set");

    // Initialize tensors
    system->resizeTensors();
}

void
GSDUtil::readParameters()
{
    spdlog::get(m_logName)->info("GSD parsing dt");
    float dt{-1.0};
    auto  return_bool = readChunk(&dt, m_frame, "log/integrator/dt", 4);
    system->setReturnBool(return_bool);
    system->check_gsd_return();
    system->setDt(double(dt));
    spdlog::get(m_logName)->info("dt : {0}", dt);
    assert(double(dt) == system->dt() && "dt not properly set");

    spdlog::get(m_logName)->info("GSD parsing t");
    float t{-1.0};
    return_bool = readChunk(&t, m_frame, "log/integrator/t", 4);
    system->setReturnBool(return_bool);
    system->check_gsd_return();
    system->setT(double(t));
    spdlog::get(m_logName)->info("t : {0}", t);
    assert(double(t) == system->t() && "t not properly set");

    spdlog::get(m_logName)->info("GSD parsing tf");
    float tf{-1.0};
    return_bool = readChunk(&tf, m_frame, "log/integrator/tf", 4);
    system->setReturnBool(return_bool);
    system->check_gsd_return();
    system->setTf(double(tf));
    spdlog::get(m_logName)->info("tf : {0}", tf);
    assert(double(tf) == system->tf() && "tf not properly set");

    spdlog::get(m_logName)->info("GSD parsing num_steps_output");
    uint64_t num_steps_output{0};
    return_bool = readChunk(&num_steps_output, m_frame, "log/integrator/num_steps_output", 8);
    system->setReturnBool(return_bool);
    system->check_gsd_return();
    system->setNumStepsOutput(int(num_steps_output));
    spdlog::get(m_logName)->info("num_steps_output : {0}", num_steps_output);
    assert(int(num_steps_output) == system->numStepsOutput() &&
           "num_steps_output not properly set");

    spdlog::get(m_logName)->info("GSD parsing fluid_density");
    float fluid_density{-1.0};
    return_bool = readChunk(&fluid_density, m_frame, "log/material_parameters/fluid_density", 4);
    system->setReturnBool(return_bool);
    system->check_gsd_return();
    system->setFluidDensity(double(fluid_density));
    spdlog::get(m_logName)->info("fluid_density : {0}", fluid_density);
    assert(double(fluid_density) == system->fluidDensity() && "fluid_density not properly set");

    spdlog::get(m_logName)->info("GSD parsing particle_density");
    float particle_density{-1.0};
    return_bool =
        readChunk(&particle_density, m_frame, "log/material_parameters/particle_density", 4);
    system->setReturnBool(return_bool);
    system->check_gsd_return();
    system->setParticleDensity(double(particle_density));
    spdlog::get(m_logName)->info("particle_density : {0}", particle_density);
    assert(double(particle_density) == system->particleDensity() &&
           "particle_density not properly set");

    spdlog::get(m_logName)->info("GSD parsing wca_epsilon");
    float wca_epsilon{-1.0};
    return_bool = readChunk(&wca_epsilon, m_frame, "log/wca/epsilon", 4);
    system->setReturnBool(return_bool);
    system->check_gsd_return();
    system->setWcaEpsilon(double(wca_epsilon));
    spdlog::get(m_logName)->info("wca_epsilon : {0}", wca_epsilon);
    assert(double(wca_epsilon) == system->wcaEpsilon() && "wca_epsilon not properly set");

    spdlog::get(m_logName)->info("GSD parsing wca_sigma");
    float wca_sigma{-1.0};
    return_bool = readChunk(&wca_sigma, m_frame, "log/wca/sigma", 4);
    system->setReturnBool(return_bool);
    system->check_gsd_return();
    system->setWcaSigma(double(wca_sigma));
    spdlog::get(m_logName)->info("wca_sigma : {0}", wca_sigma);
    assert(double(wca_sigma) == system->wcaSigma() && "wca_sigma not properly set");
}

void
GSDUtil::readParticles()
{
    // positions
    spdlog::get(m_logName)->info("GSD parsing position");
    float pos[system->numDim() * system->numParticles()];
    auto  return_bool = readChunk(&pos, m_frame, "particles/position", system->numParticles() * 12,
                                 system->numParticles());
    system->setReturnBool(return_bool);
    system->check_gsd_return();
    for (int i = 0; i < system->numParticles(); i++)
    {
        (*system->positions())(3 * i)     = pos[3 * i];
        (*system->positions())(3 * i + 1) = pos[3 * i + 1];
        (*system->positions())(3 * i + 2) = pos[3 * i + 2];
        spdlog::get(m_logName)->info("Particle {0} position : [{1:03.3f}, {2:03.3f}, {3:03.3f}]",
                                     i + 1, pos[system->numDim() * i],
                                     pos[system->numDim() * i + 1], pos[system->numDim() * i + 2]);
    }

    // velocities
    spdlog::get(m_logName)->info("GSD parsing velocity");
    float vel[system->numDim() * system->numParticles()];
    return_bool = readChunk(&vel, m_frame, "particles/velocity", system->numParticles() * 12,
                            system->numParticles());
    system->setReturnBool(return_bool);
    system->check_gsd_return();
    for (int i = 0; i < system->numParticles(); i++)
    {
        (*system->velocities())(3 * i)     = vel[3 * i];
        (*system->velocities())(3 * i + 1) = vel[3 * i + 1];
        (*system->velocities())(3 * i + 2) = vel[3 * i + 2];
        spdlog::get(m_logName)->info("Particle {0} velocity : [{1:03.3f}, {2:03.3f}, {3:03.3f}]",
                                     i + 1, vel[system->numDim() * i],
                                     vel[system->numDim() * i + 1], vel[system->numDim() * i + 2]);
    }

    // accelerations
    spdlog::get(m_logName)->info("GSD parsing acceleration");
    float acc[system->numDim() * system->numParticles()];
    return_bool = readChunk(&acc, m_frame, "particles/moment_inertia", system->numParticles() * 12,
                            system->numParticles());
    system->setReturnBool(return_bool);
    system->check_gsd_return();
    for (int i = 0; i < system->numParticles(); i++)
    {
        (*system->accelerations())(3 * i)     = acc[3 * i];
        (*system->accelerations())(3 * i + 1) = acc[3 * i + 1];
        (*system->accelerations())(3 * i + 2) = acc[3 * i + 2];
        spdlog::get(m_logName)->info(
            "Particle {0} acceleration : [{1:03.3f}, {2:03.3f}, {3:03.3f}]", i + 1,
            acc[system->numDim() * i], acc[system->numDim() * i + 1],
            acc[system->numDim() * i + 2]);
    }
}

void
GSDUtil::writeFrame()
{
    spdlog::get(m_logName)->info("GSD writing frame");
    spdlog::get(m_logName)->info("time step: {0}", system->timestep());
    spdlog::get(m_logName)->info("time: {0}", system->t());

    writeHeader();
    writeParameters();
    writeParticles();

    spdlog::get(m_logName)->info("GSD ending frame");
    auto return_val = gsd_end_frame(system->handle().get());
    system->setReturnVal(return_val);
    system->check_gsd_return();
}

void
GSDUtil::writeHeader()
{
    spdlog::get(m_logName)->info("GSD writing log/configuration/timestep");
    uint64_t step       = system->timestep();
    auto     return_val = gsd_write_chunk(system->handle().get(), "log/configuration/step",
                                      GSD_TYPE_UINT64, 1, 1, 0, (void*)&step);
    system->setReturnVal(return_val);
    system->check_gsd_return();

    if (gsd_get_nframes(system->handle().get()) == 0)
    {
        spdlog::get(m_logName)->info("GSD writing configuration/dimensions");
        uint8_t dimensions = system->numDim();
        return_val         = gsd_write_chunk(system->handle().get(), "configuration/dimensions",
                                     GSD_TYPE_UINT8, 1, 1, 0, (void*)&dimensions);
        system->setReturnVal(return_val);
        system->check_gsd_return();
    }

    spdlog::get(m_logName)->info("GSD writing particles/N");
    uint32_t N = system->numParticles();
    return_val =
        gsd_write_chunk(system->handle().get(), "particles/N", GSD_TYPE_UINT32, 1, 1, 0, (void*)&N);
    system->setReturnVal(return_val);
    system->check_gsd_return();
}

void
GSDUtil::writeParameters()
{
    spdlog::get(m_logName)->info("GSD writing log/integrator/dt");
    float_t dt      = system->dt();
    auto return_val = gsd_write_chunk(system->handle().get(), "log/integrator/dt", GSD_TYPE_FLOAT,
                                      1, 1, 0, (void*)&dt);
    system->setReturnVal(return_val);
    system->check_gsd_return();

    spdlog::get(m_logName)->info("GSD writing log/integrator/t");
    float_t time = system->t();
    return_val   = gsd_write_chunk(system->handle().get(), "log/integrator/t", GSD_TYPE_FLOAT, 1, 1,
                                 0, (void*)&time);
    system->setReturnVal(return_val);
    system->check_gsd_return();
}

void
GSDUtil::writeParticles()
{
    uint32_t N       = system->numParticles();
    uint32_t num_dim = system->numDim();
    int      return_val;
    uint64_t nframes = gsd_get_nframes(system->handle().get());
    spdlog::get(m_logName)->debug("vectors are assumed to have 3 spatial DoF");

    std::vector<float> pos(uint64_t(N) * uint64_t(num_dim));
    pos.reserve(1); //! make sure we allocate
    for (uint32_t i = 0; i < N; i++)
    {
        pos[3 * i]     = (*system->positions())(3 * i);
        pos[3 * i + 1] = (*system->positions())(3 * i + 1);
        pos[3 * i + 2] = (*system->positions())(3 * i + 2);
    }
    spdlog::get(m_logName)->info("GSD writing particles/position");
    return_val = gsd_write_chunk(system->handle().get(), "particles/position", GSD_TYPE_FLOAT, N, 3,
                                 0, (void*)&pos[0]);
    system->setReturnVal(return_val);
    system->check_gsd_return();

    std::vector<float> vel(uint64_t(N) * uint64_t(num_dim));
    vel.reserve(1); //! make sure we allocate
    for (uint32_t i = 0; i < N; i++)
    {
        vel[3 * i]     = (*system->velocities())(3 * i);
        vel[3 * i + 1] = (*system->velocities())(3 * i + 1);
        vel[3 * i + 2] = (*system->velocities())(3 * i + 2);
    }
    spdlog::get(m_logName)->info("GSD writing particles/velocity");
    return_val = gsd_write_chunk(system->handle().get(), "particles/velocity", GSD_TYPE_FLOAT, N, 3,
                                 0, (void*)&vel[0]);
    system->setReturnVal(return_val);
    system->check_gsd_return();

    std::vector<float> acc(uint64_t(N) * uint64_t(num_dim));
    acc.reserve(1); //! make sure we allocate
    for (uint32_t i = 0; i < N; i++)
    {
        acc[3 * i]     = (*system->accelerations())(3 * i);
        acc[3 * i + 1] = (*system->accelerations())(3 * i + 1);
        acc[3 * i + 2] = (*system->accelerations())(3 * i + 2);
    }
    spdlog::get(m_logName)->info("GSD writing particles/acceleration");
    return_val = gsd_write_chunk(system->handle().get(), "particles/acceleration", GSD_TYPE_FLOAT,
                                 N, 3, 0, (void*)&acc[0]);
    system->setReturnVal(return_val);
    system->check_gsd_return();
}