//
// Created by Alec Glisman on 07/30/21
// Building off of GSDReader class in HOOMD-Blue
//

#include <GSDReader.hpp>

GSDReader::GSDReader(std::shared_ptr<systemData> sys)
{
    // save classes
    system = sys;

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
    readParticles();
}

GSDReader::GSDReader(std::shared_ptr<systemData> sys, uint64_t frame) : GSDReader(sys)
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

GSDReader::~GSDReader() = default;

bool
GSDReader::readChunk(void* data, uint64_t frame, const char* name, size_t expected_size,
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
GSDReader::readHeader()
{
    uint64_t timestep    = 0;
    auto     return_bool = readChunk(&timestep, m_frame, "configuration/step", 8);
    system->setReturnBool(return_bool);
    system->check_gsd_return();
    system->setTimestep(timestep);
    spdlog::get(m_logName)->info("time step: {0}", timestep);

    uint8_t dim = 3; // FIXME: change to 0 when parsing stable
    return_bool = readChunk(&dim, m_frame, "configuration/dimensions", 1);
    system->setReturnBool(return_bool);
    system->check_gsd_return();
    system->setNumDim(dim);
    spdlog::get(m_logName)->info("dim : {0}", dim);

    unsigned int N = 0;
    return_bool    = readChunk(&N, m_frame, "particles/N", 4);
    system->setReturnBool(return_bool);
    system->check_gsd_return();
    system->setNumParticles(N);
    spdlog::get(m_logName)->info("Number of particles : {0}", N);
    if (N == 0)
    {
        spdlog::get(m_logFile)->error("data.gsd_snapshot: cannot read a file with 0 particles");
        throw std::runtime_error("Error reading GSD file");
    }

    double dt{0.0};
    return_bool = readChunk(&dt, m_frame, "integrator/dt", 32);
    system->setReturnBool(return_bool);
    system->check_gsd_return();
    system->setDt(dt);
    spdlog::get(m_logName)->info("dt : {0}", dt);

    // Initialize tensors
    system->resizeTensors();
}

void
GSDReader::readParticles()
{
    // positions
    double_t pos[system->numDim() * system->numParticles()];
    auto return_bool = readChunk(&pos, m_frame, "particles/position", system->numParticles() * 12,
                                 system->numParticles());
    system->setReturnBool(return_bool);
    system->check_gsd_return();
    // system->positions(Eigen::VectorXd{pos}); //FIXME: figure out how to convert data types
    for (int i = 0; i < system->numParticles(); i++)
    {
        spdlog::get(m_logName)->info("Particle {0} position : [{1}, {2}, {3}]", i + 1,
                                     pos[system->numDim() * i], pos[system->numDim() * i + 1],
                                     pos[system->numDim() * i + 2]);
    }

    // velocities
    double_t vel[system->numDim() * system->numParticles()];
    return_bool = readChunk(&vel, m_frame, "particles/velocity", system->numParticles() * 12,
                            system->numParticles());
    system->setReturnBool(return_bool);
    system->check_gsd_return();
    // system->velocities(Eigen::VectorXd{vel}); //FIXME: figure out how to convert data types
    for (int i = 0; i < system->numParticles(); i++)
    {
        spdlog::get(m_logName)->info("Particle {0} velocity : [{1}, {2}, {3}]", i + 1,
                                     vel[system->numDim() * i], vel[system->numDim() * i + 1],
                                     vel[system->numDim() * i + 2]);
    }

    // accelerations
    double_t acc[system->numDim() * system->numParticles()];
    return_bool = readChunk(&acc, m_frame, "particles/inertia", system->numParticles() * 12,
                            system->numParticles());
    system->setReturnBool(return_bool);
    system->check_gsd_return();
    // system->accelerations(Eigen::VectorXd{acc}); //FIXME: figure out how to convert data types
    for (int i = 0; i < system->numParticles(); i++)
    {
        spdlog::get(m_logName)->info("Particle {0} acceleration : [{1}, {2}, {3}]", i + 1,
                                     acc[system->numDim() * i], acc[system->numDim() * i + 1],
                                     acc[system->numDim() * i + 2]);
    }
}
