//
// Created by Alec Glisman on 3/27/20.
//

/* Include all internal project dependencies */

/* Include all external project dependencies */

// eigen3(Linear algebra)
#include <Eigen/Core>
#include <Eigen/Eigen>
#define EIGEN_USE_MKL_ALL
#include <Eigen/src/Core/util/MKL_support.h>

// spdlog (logging)
#include <spdlog/spdlog.h>

// Catch2
#include <catch2/catch.hpp>

// STL
#include <memory> // for std::unique_ptr and std::shared_ptr
#include <string> // std::string
#include <vector> // std::vector()

/* Forward declarations */

int
main(const int argc, const char* argv[])
{

    /* SECTION: Parse command line input
     *      argv[0]: executable name
     *      argv[1]: input preferences filepath
     *      argv[2]: output directory to write data to
     *      argv[3]: (optional) initial configuration data file
     */

    // Check input is of correct length
    // if (argc < 3)
    // {

    //     std::cout << "ERROR: incorrect number of arguments!!!"
    //               << "[executable] [input preferences] [output directory] [(optional) initial "
    //                  "configuration]"
    //               << std::endl;

    //     return EXIT_FAILURE;
    // }

    /* !SECTION */

    return EXIT_SUCCESS;
}
