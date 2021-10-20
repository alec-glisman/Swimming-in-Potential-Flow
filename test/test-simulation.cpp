//
// Created by Alec Glisman on 07/30/21
//

/* Include all internal project dependencies */
#include <engine.hpp>
#include <systemData.hpp>

/* Include all external project dependencies */
#define CATCH_CONFIG_CONSOLE_WIDTH 300
#include <catch2/catch.hpp> // unit testing framework
#include <gsd.h>            // GSD File
// STL
#include <fstream> // std::ifstream
#include <memory>  // for std::unique_ptr and std::shared_ptr
#include <string>  // std::string

SCENARIO("Open GSD file", "[gsd]")
{
    GIVEN("Collinear wall system initial conditions GSD input.")
    {
        // I/O Parameters
        std::string inputDataFile = "input/collinear_swimmer_wall/initial_frame_dt1e-6_Z-height6.gsd";
        std::string outputDir     = "output";

        // GSD loading parameters
        std::shared_ptr<gsd_handle> handle{new gsd_handle};
        int                         return_val{-1};

        // check inputDataFile exists
        std::ifstream inputDataStream(inputDataFile);
        REQUIRE(inputDataStream);

        WHEN("GSD is opened")
        {
            REQUIRE_NOTHROW(return_val = gsd_open(handle.get(), inputDataFile.c_str(), GSD_OPEN_READONLY));

            THEN("File opened successfully")
            {
                INFO("GSD open return value: " << return_val);
                REQUIRE(return_val == 0);
            }
        }
    }
}

SCENARIO("Run simulation system",
         "[engine][progressBar][systemData][rungeKutta4][potentialHydrodynamics][gsd][GSDUtil]")
{
    GIVEN("Collinear wall system initial conditions GSD input.")
    {
        // I/O Parameters
        std::string inputDataFile = "input/collinear_swimmer_wall/initial_frame_dt1e-6_Z-height6.gsd";
        std::string outputDir     = "output";

        // simulation classes
        std::shared_ptr<systemData> system;
        std::shared_ptr<engine>     eng;

        WHEN("simulationSystem initialized")
        {
            REQUIRE_NOTHROW(system = std::make_shared<systemData>(inputDataFile, outputDir));
            REQUIRE_NOTHROW(system->initializeData());

            THEN("Initialization successful")
            {
                REQUIRE(system->gSDParsed());
            }
        }

        WHEN("engine initialized and run")
        {
            THEN("Time integration completes successfully")
            {
                REQUIRE_NOTHROW(eng = std::make_shared<engine>(system));
                REQUIRE_NOTHROW(eng->run());
            }
        }
    }
}
