//
// Created by Alec Glisman on 07/30/21
//

/* Include all internal project dependencies */

/* Include all external project dependencies */
#include <catch2/catch.hpp> // unit testing framework
#include <gsd.h>            // GSD File
// STL
#include <memory> // for std::unique_ptr and std::shared_ptr
#include <string> // std::string

/* Forward declarations */

TEST_CASE("Load and parse data from GSD file", "[GSD][I/O]")
{
    /* Parameters
     * REVIEW[epic=Assumptions]: Assuming that the tests are called from project base dir */
    std::string inputDataFile{"input/test.gsd"};
    std::string outputDir{"temp/output"};

    /* Initialize variables */
    std::shared_ptr<gsd_handle> handle{new gsd_handle};
    int                         return_val{-1};

    SECTION("Load GSD data")
    {
        REQUIRE_NOTHROW(return_val =
                            gsd_open(handle.get(), inputDataFile.c_str(), GSD_OPEN_READONLY));
        REQUIRE(return_val == 0);
    }
}
