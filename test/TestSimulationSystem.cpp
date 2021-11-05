
//
// Created by Alec Glisman on 10/21/21
//

/* Include all internal project dependencies */
#include <TestSystemData.hpp>

/* Include all external project dependencies */
#define CATCH_CONFIG_CONSOLE_WIDTH 300
#include <catch2/catch.hpp> // unit testing framework
// STL
#include <memory> // for std::unique_ptr and std::shared_ptr
#include <string> // std::string

TEST_CASE("Test SystemData class", "[SystemData]")
{
    // I/O Parameters
    std::string inputDataFile = "input/collinear_swimmer_wall/initial_frame_dt1e-1_Z-height6.gsd";
    std::string outputDir     = "output-SystemData-init";

    // simulation classes
    std::shared_ptr<SystemData>     system;
    std::shared_ptr<TestSystemData> testSystem;

    // Construct systenData class
    REQUIRE_NOTHROW(system = std::make_shared<SystemData>(inputDataFile, outputDir));
    // Construct TestSystemData class
    REQUIRE_NOTHROW(testSystem = std::make_shared<TestSystemData>(system));

    // Return value test
    int return_val{-1};

    SECTION("Test static helper functions")
    {
        REQUIRE_NOTHROW(return_val = testSystem->testCrossProdMat());
        REQUIRE(return_val == 0);

        REQUIRE_NOTHROW(return_val = testSystem->testEMatrix());
        REQUIRE(return_val == 0);
    }

    // REQUIRE_NOTHROW(system->initializeData());
    // Verify data was correctly parsed from GSD to simulation
    // REQUIRE(system->gSDParsed());
}