//
// Created by Alec Glisman on 07/30/21
//

/* Include all internal project dependencies */

/* Include all external project dependencies */
#include <catch2/catch.hpp> // unit testing framework
#include <gsd.h>            // GSD File
// STL
#include <string> // std::string

/* Forward declarations */

TEST_CASE("Load and parse data from GSD file", "[GSD][I/O]")
{
    SECTION("Load GSD data")
    {
        /* Parameters
         * REVIEW[epic=Assumptions]: Assuming that the tests are called from project base dir */
        std::string inputData{"input/three_body_collinear_x/test.gsd"};
        std::string outputDir{"temp/output"};

        /* Initialize variables */
    }
}
