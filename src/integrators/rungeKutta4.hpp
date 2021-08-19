//
// Created by Alec Glisman on 07/30/21
// Building off of GSDReader class in HOOMD-Blue
//

#ifndef BODIES_IN_POTENTIAL_FLOW_RUNGE_KUTTA_4_H
#define BODIES_IN_POTENTIAL_FLOW_RUNGE_KUTTA_4_H

#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif
/* Include all internal project dependencies */
#include <potentialHydrodynamics.hpp>
#include <systemData.hpp>

/* Include all external project dependencies */
// eigen3 (Linear algebra)
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Eigenvalues>
#define EIGEN_USE_MKL_ALL
#include <eigen3/Eigen/src/Core/util/MKL_support.h>
// Logging
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/spdlog.h>

/* Forward declarations */
class systemData;

class rungeKutta4
{
  public:
    rungeKutta4(std::shared_ptr<systemData> sys, std::shared_ptr<potentialHydrodynamics> hydro);

    ~rungeKutta4();

    void
    integrate();

  private:
    void
    initializeSpecificVars();

    void
    articulationVel(double dimensional_time);

    void
    articulationAcc(double dimensional_time);

    void
    rLoc();

    void
    momentumLinAngFree();

    void
    crossProdMat(const Eigen::Vector3d& vec, Eigen::Matrix3d& mat)
    {
        /* Function takes in vector in vector cross-product expression: c = a \times b
         * vec must be 'a' in above equation
         * Output is the matrix representation of 'a \times' operator
         */
        mat << 0, -vec(2), vec(1), vec(2), 0, -vec(0), -vec(1), vec(0), 0;
    };

    void
    accelerationUpdate(Eigen::VectorXd& acc, double dimensional_time);

    // classes
    std::shared_ptr<systemData>             m_system;
    std::shared_ptr<potentialHydrodynamics> m_potHydro;

    // logging
    std::string m_logFile;
    std::string m_logName{"rungeKutta4"};

    // swimmer parameters
    Eigen::Vector3d m_RLoc;
    Eigen::VectorXd m_velArtic;
    Eigen::VectorXd m_accArtic;

    // specific parameters
    double m_U0{-1.0};
    double m_omega{-1.0};
    double m_phase_shift{-1.0};
    double m_R_avg{-1.0};

    // time step variables
    double m_dt{-1.0};
    double m_c1_2_dt{-1.0};
    double m_c1_6_dt{-1.0};

    // constants
    const double m_c1_2{0.50};
    const double m_c1_6{1.0 / 6.0};
    const double m_c2{2.0};
};

#endif // BODIES_IN_POTENTIAL_FLOW_RUNGE_KUTTA_4_H