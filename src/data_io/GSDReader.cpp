//
// Created by Alec Glisman on 07/30/21
// Building off of GSDReader class in HOOMD-Blue
//

#include <GSDReader.hpp>

GSDReader::GSDReader(const std::string& name, const uint64_t frame) : m_name(name), m_frame(frame)
{
    // // open the GSD file in read mode
    // spdlog::info("Welcome to spdlog!");
    // m_exec_conf->msg->notice(3) << "data.gsd_snapshot: open gsd file " << name << endl;
    // int retval = gsd_open(&m_handle, name.c_str(), GSD_OPEN_READONLY);
    // checkError(retval);

    // // validate schema
    // if (string(m_handle.header.schema) != string("hoomd"))
    // {
    //     m_exec_conf->msg->error() << "data.gsd_snapshot: "
    //                               << "Invalid schema in " << name << endl;
    //     throw runtime_error("Error opening GSD file");
    // }
    // if (m_handle.header.schema_version >= gsd_make_version(2, 0))
    // {
    //     m_exec_conf->msg->error() << "data.gsd_snapshot: "
    //                               << "Invalid schema version in " << name << endl;
    //     throw runtime_error("Error opening GSD file");
    // }

    // // set frame from the end of the file if requested
    // uint64_t nframes = gsd_get_nframes(&m_handle);
    // if (from_end && frame <= nframes)
    //     m_frame = nframes - frame;

    // // validate number of frames
    // if (m_frame >= nframes)
    // {
    //     m_exec_conf->msg->error() << "data.gsd_snapshot: "
    //                               << "Cannot read frame " << m_frame << " " << name << " only has
    //                               "
    //                               << gsd_get_nframes(&m_handle) << " frames" << endl;
    //     throw runtime_error("Error opening GSD file");

    readHeader();
    readParticles();
    readTopology();
}

GSDReader::~GSDReader()
{
    gsd_close(&m_handle);
}
