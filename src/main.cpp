//
// Created by Alec Glisman on 3/27/20.
//

/* Include all internal project dependencies */

/* Include all external project dependencies */

// // eigen3(Linear algebra)
// #include <Eigen/Core>
// #include <Eigen/Eigen>
// #define EIGEN_USE_MKL_ALL
// #include <Eigen/src/Core/util/MKL_support.h>

// // spdlog (logging)
// #include <spdlog/spdlog.h>

// // Catch2
// #include <catch2/catch.hpp>

// // STL
#include <iostream> // std::cout, std::endl
#include <string>   // std::string, and EXIT_FAILURE/EXIT_SUCCESS
// #include <memory> // for std::unique_ptr and std::shared_ptr
// #include <vector> // std::vector()

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

    /* !SECTION */

    return EXIT_SUCCESS;
}
