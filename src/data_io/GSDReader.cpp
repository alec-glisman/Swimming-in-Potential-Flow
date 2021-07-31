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

    readHeader();
    readParticles();
}

GSDReader::GSDReader(std::shared_ptr<systemData> sys, uint64_t frame) : GSDReader(sys)
{
    m_frame = frame;
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

        auto retval = gsd_read_chunk(system->handle().get(), data, entry);
        system->setReturnVal(retval);
        system->check_gsd_return();
    }

    return true;
}

void
GSDReader::readHeader()
{
    // readChunk(&m_timestep, m_frame, "configuration/step", 8);

    // uint8_t dim = 3;
    // readChunk(&dim, m_frame, "configuration/dimensions", 1);
    // m_snapshot->dimensions = dim;

    // float box[6] = {1.0f, 1.0f, 1.0f, 0.0f, 0.0f, 0.0f};
    // readChunk(&box, m_frame, "configuration/box", 6 * 4);
    // m_snapshot->global_box = BoxDim(box[0], box[1], box[2]);
    // m_snapshot->global_box.setTiltFactors(box[3], box[4], box[5]);

    // unsigned int N = 0;
    // readChunk(&N, m_frame, "particles/N", 4);
    // if (N == 0)
    // {
    //     m_exec_conf->msg->error() << "data.gsd_snapshot: "
    //                               << "cannot read a file with 0 particles" << endl;
    //     throw runtime_error("Error reading GSD file");
    // }

    // m_snapshot->particle_data.resize(N);
}

void
GSDReader::readParticles()
{
}
