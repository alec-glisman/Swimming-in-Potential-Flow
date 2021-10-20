//
// Created by Alec Glisman on 07/30/21
//

/* Include all internal project dependencies */
#include <engine.hpp>
#include <systemData.hpp>

/* Include all external project dependencies */
#include <catch2/catch.hpp> // unit testing framework
#include <gsd.h>            // GSD File
// STL
#include <memory> // for std::unique_ptr and std::shared_ptr
#include <string> // std::string

TEST_CASE("Load and parse data from GSD file", "[GSD][I/O]")
{
    /* I/O Parameters */
    std::string inputDataFile = "input/collinear_swimmer_wall/initial_frame_dt1e-6_Z-height6.gsd";
    std::string outputDir     = "output";

    SECTION("Load GSD data")
    {
        /* Initialize variables */
        std::shared_ptr<gsd_handle> handle{new gsd_handle};
        int                         return_val{-1};

        REQUIRE_NOTHROW(return_val = gsd_open(handle.get(), inputDataFile.c_str(), GSD_OPEN_READONLY));
        REQUIRE(return_val == 0);
    }

    SECTION("Initialize simulation")
    {
        /* Initialize variables */
        std::shared_ptr<systemData> system;
        std::shared_ptr<engine>     eng;

        REQUIRE_NOTHROW(system = std::make_shared<systemData>(inputDataFile, outputDir));
        REQUIRE_NOTHROW(system->initializeData());

        REQUIRE_NOTHROW(eng = std::make_shared<engine>(system));
    }
}
