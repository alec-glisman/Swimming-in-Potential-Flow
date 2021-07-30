//
// Created by Alec Glisman on 3/27/20.
//

/* Include all internal project dependencies */
#include <systemData.hpp>

/* Include all external project dependencies */
// STL
#include <iostream> // std::cout, std::endl
#include <memory>   // for std::unique_ptr and std::shared_ptr
#include <string>   // std::string
// GSD: Data IO
#include <gsd.h> // GSD File

/* Forward declarations */

int
main(const int argc, const char* argv[])
{

    /* SECTION: Parse command line input
     *      argv[0]: executable name
     *      argv[1]: input GSD filepath
     *      argv[2]: output directory to write data
     */

    // Check input is of correct length
    if (argc < 3)
    {

        std::cout << "ERROR: incorrect number of arguments!!!"
                  << "[executable] [input data] [output directory]" << std::endl;

        return EXIT_FAILURE;
    }

    // Get input files
    std::string inputDataFile, outputDir;

    inputDataFile = argv[1];
    outputDir     = argv[2];
    /* !SECTION */

    /* SECTION: Set-up and run simulation */
    // Initialize data structures
    auto system = std::make_shared<systemData>(inputDataFile, outputDir);
    /* !SECTION */

    return EXIT_SUCCESS;
}
