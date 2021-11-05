//
// Created by Alec Glisman on 3/27/20.
//

/* Include all internal project dependencies */
#include <Engine.hpp>
#include <SystemData.hpp>

/* Include all external project dependencies */
// STL
#include <memory> // for std::unique_ptr and std::shared_ptr
#include <string> // std::string

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
        throw std::runtime_error("ERROR: incorrect number of arguments!!! [executable][input data][output directory]");
    }

    // Get input files
    std::string inputDataFile, outputDir;

    inputDataFile = argv[1];
    outputDir     = argv[2];
    /* !SECTION */

    /* SECTION: Set-up and run simulation */
    // Initialize data structures
    auto system = std::make_shared<SystemData>(inputDataFile, outputDir);
    system->initializeData();
    auto eng = std::make_shared<Engine>(system);

    // Run simulations
    eng->run();

    /* !SECTION */

    return EXIT_SUCCESS;
}
