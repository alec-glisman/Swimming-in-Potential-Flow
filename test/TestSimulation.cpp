//
// Created by Alec Glisman on 07/30/21
//

/* Include all internal project dependencies */
#include <Engine.hpp>
#include <SystemData.hpp>

/* Include all external project dependencies */
#define CATCH_CONFIG_CONSOLE_WIDTH 300
#include <catch2/catch.hpp> // unit testing framework
#include <gsd.h>            // GSD File
// Logging
#include <spdlog/fmt/ostr.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/spdlog.h>
// STL
#include <fstream> // std::ifstream
#include <memory>  // for std::unique_ptr and std::shared_ptr
#include <string>  // std::string

TEST_CASE("Open GSD file", "[gsd]")
{
    // I/O Parameters
    std::string inputDataFile = "test/input/collinear_swimmer_isolated/initial_frame_dt1e-2.gsd";

    // GSD loading parameters
    std::shared_ptr<gsd_handle> handle{new gsd_handle};
    int                         return_val{-1};

    // verify inputDataFile exists
    std::ifstream inputDataStream(inputDataFile);
    REQUIRE(inputDataStream);

    // verify gsd file can be opened correctly
    REQUIRE_NOTHROW(return_val = gsd_open(handle.get(), inputDataFile.c_str(), GSD_OPEN_READONLY));
    INFO("GSD open return value: " << return_val);
    REQUIRE(return_val == 0);
}

TEST_CASE("Initialize collinear swimmer isolated simulation system",
          "[Collinear-Isolated][Engine][SystemData][RungeKutta4][PotentialHydrodynamics][gsd][GSDUtil]")
{
    // close all previous loggers
    spdlog::drop_all();

    // I/O Parameters
    std::string inputDataFile = "test/input/collinear_swimmer_isolated/initial_frame_dt1e-2.gsd";
    std::string outputDir     = "output-collinear-isolated-SystemData-init";

    // simulation classes
    std::shared_ptr<SystemData> system;
    std::shared_ptr<Engine>     eng;

    // Construct and initialize simulationSystem class
    REQUIRE_NOTHROW(system = std::make_shared<SystemData>(inputDataFile, outputDir));
    REQUIRE_NOTHROW(system->initializeData());

    // Verify data was correctly parsed from GSD to simulation
    REQUIRE(system->gSDParsed());

    // Construct simulation Engine
    REQUIRE_NOTHROW(eng = std::make_shared<Engine>(system));
}

TEST_CASE("Run collinear swimmer isolated simulation system",
          "[Collinear-Isolated][Engine][ProgressBar][SystemData][RungeKutta4][PotentialHydrodynamics]")
{
    // close all previous loggers
    spdlog::drop_all();

    // I/O Parameters
    std::string inputDataFile = "test/input/collinear_swimmer_isolated/initial_frame_dt1e-2.gsd";
    std::string outputDir     = "output-collinear-isolated-Engine-run";

    // simulation classes
    std::shared_ptr<SystemData> system;
    std::shared_ptr<Engine>     eng;

    // Construct and initialize simulation classes
    REQUIRE_NOTHROW(system = std::make_shared<SystemData>(inputDataFile, outputDir));
    REQUIRE_NOTHROW(system->initializeData());
    REQUIRE_NOTHROW(eng = std::make_shared<Engine>(system));

    // Verify simulation can run without error
    REQUIRE_NOTHROW(eng->run());
}

TEST_CASE("Initialize collinear swimmer wall simulation system",
          "[Collinear-Wall][Engine][SystemData][RungeKutta4][PotentialHydrodynamics][gsd][GSDUtil]")
{
    // close all previous loggers
    spdlog::drop_all();

    // I/O Parameters
    std::string inputDataFile = "input/collinear_swimmer_wall/initial_frame_dt1e-1_Z-height6.gsd";
    std::string outputDir     = "output-collinear-wall-SystemData-init";

    // simulation classes
    std::shared_ptr<SystemData> system;
    std::shared_ptr<Engine>     eng;

    // Construct and initialize simulationSystem class
    REQUIRE_NOTHROW(system = std::make_shared<SystemData>(inputDataFile, outputDir));
    REQUIRE_NOTHROW(system->initializeData());

    // Verify data was correctly parsed from GSD to simulation
    REQUIRE(system->gSDParsed());

    // Construct simulation Engine
    REQUIRE_NOTHROW(eng = std::make_shared<Engine>(system));
}

TEST_CASE("Run collinear swimmer wall simulation system",
          "[Collinear-Wall][Engine][ProgressBar][SystemData][RungeKutta4][PotentialHydrodynamics]")
{
    // close all previous loggers
    spdlog::drop_all();

    // I/O Parameters
    std::string inputDataFile = "input/collinear_swimmer_wall/initial_frame_dt1e-1_Z-height6.gsd";
    std::string outputDir     = "output-collinear-wall-Engine-run";

    // simulation classes
    std::shared_ptr<SystemData> system;
    std::shared_ptr<Engine>     eng;

    // Construct and initialize simulation classes
    REQUIRE_NOTHROW(system = std::make_shared<SystemData>(inputDataFile, outputDir));
    REQUIRE_NOTHROW(system->initializeData());
    REQUIRE_NOTHROW(eng = std::make_shared<Engine>(system));

    // Verify simulation can run without error
    REQUIRE_NOTHROW(eng->run());
}
