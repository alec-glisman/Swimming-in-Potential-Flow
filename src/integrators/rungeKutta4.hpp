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

    void
    GSDOutput();

  private:
    /* REVIEW[epic=Change]: Change initializeSpecificVars() for different systems and load data
     * into m_systemParam */
    void
    initializeSpecificVars();

    /* REVIEW[epic=Change]: Change constraint linear system (A, b) for each system. May need to
     * change dimensions of Eigen matrixes and vectors */
    void
    initializeConstraintLinearSystem();

    void
    updateConstraintLinearSystem(double dimensional_time);

    void
    accelerationUpdate(Eigen::VectorXd& acc, double dimensional_time);

    void
    udwadiaKalaba(Eigen::VectorXd& acc, double dimensional_time);

    /* REVIEW[epic=Change]: Change rbmconn if particles do not all move identically with respect to
     * a single locater point */
    /* NOTE: before calling this method, make sure to call articulationVel(), articulationAcc(), and
     * rLoc() */
    void
    momentumLinAngFree(Eigen::VectorXd& acc, double dimensional_time);

    /* REVIEW[epic=Change]: Change assignment of m_velArtic for different systems */
    void
    articulationVel(double dimensional_time);

    /* REVIEW[epic=Change]: Change assignment of m_accArtic for different systems */
    void
    articulationAcc(double dimensional_time);

    /* REVIEW[epic=Change]: Change assignment of m_RLoc for different systems */
    void
    rLoc();

    /* Function takes in vector in vector cross-product expression: c = a \times b
     * vec must be 'a' in above equation
     * Output is the matrix representation of 'a \times' operator */
    void
    crossProdMat(const Eigen::Vector3d& vec, Eigen::Matrix3d& mat)
    {
        mat << 0, -vec(2), vec(1), vec(2), 0, -vec(0), -vec(1), vec(0), 0;
    };

    /* REVIEW[epic=Change,order=0]: change parameters stored for different systems */
    // system specific data
    struct systemParameters
    {
        double          U0{-1.0};
        double          omega{-1.0};
        double          phaseShift{-1.0};
        double          RAvg{-1.0};
        Eigen::VectorXd U_swim;
        Eigen::VectorXd A_swim;
    };

    // classes
    std::shared_ptr<systemData>             m_system;
    std::shared_ptr<potentialHydrodynamics> m_potHydro;

    // structs
    systemParameters m_systemParam;

    // logging
    std::string       m_logFile;
    const std::string m_logName{"rungeKutta4"};

    // (linear and angular) "momentum" and "force" balance parameters
    Eigen::Vector3d m_RLoc;
    Eigen::VectorXd m_velArtic;
    Eigen::VectorXd m_accArtic;

    // constraint parameters
    Eigen::MatrixXd m_A;
    Eigen::VectorXd m_b;

    // "identity" tensors
    const Eigen::Matrix3d m_I = Eigen::Matrix3d::Identity(3, 3);
    Eigen::Matrix3d       m_I_tilde;

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