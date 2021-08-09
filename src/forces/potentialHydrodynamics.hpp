//
// Created by Alec Glisman on 07/30/21
// Building off of GSDReader class in HOOMD-Blue
//

#ifndef BODIES_IN_POTENTIAL_FLOW_POTENTIAL_HYDRODYNAMICS_H
#define BODIES_IN_POTENTIAL_FLOW_POTENTIAL_HYDRODYNAMICS_H

#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif
/* Include all internal project dependencies */
#include <systemData.hpp>

/* Include all external project dependencies */
// eigen3 (Linear algebra)
#include <Eigen/Core>
#include <Eigen/Eigen>
#define EIGEN_USE_MKL_ALL
#include <Eigen/src/Core/util/MKL_support.h>
// Logging
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/spdlog.h>

/* Forward declarations */
class systemData;

class potentialHydrodynamics
{
  public:
    potentialHydrodynamics(systemData& sys);

    ~potentialHydrodynamics();

    void
    update();

  private:
    // classes
    systemData* system;

    // logging
    std::string m_logFile;
    std::string m_logName{"potentialHydrodynamics"};

    // For-loop variables
    int num_inter{-1}; // Number of pairwise interactions to count

    // tensor variables
    int len_tensor{-1}; // length of tensor quantities

    // Distance between particle pairs
    Eigen::VectorXd alphaVec; // particle number; dim = (r) x (1)
    Eigen::VectorXd betaVec;  // particle number; dim = (r) x (1)
    Eigen::VectorXd r_mag_ab; // [1]; dim = (r) x (1)
    Eigen::MatrixXd r_ab;     // [1]; dim = (3)  x (r)

    // Identity Matrices
    Eigen::MatrixXd       I3N;                                  // dim = (3N) x (3N)
    const Eigen::Matrix3d I3 = Eigen::Matrix3d::Identity(3, 3); // dim = 3 x 3

    // hydrodynamic quantities
    Eigen::MatrixXd M_added;
    Eigen::MatrixXd M_intrinsic;
    Eigen::MatrixXd M_total;
    Eigen::MatrixXd grad_M_added;

    Eigen::VectorXd F_hydro;
    Eigen::VectorXd F_hydroNoInertia;

    // constants
    const double unitSphereVol{4.0 / 3.0 * M_PI};
    const double c1_2{0.50};
    const double c2_3{2.0 / 3.0};
    const double c3_2{1.50};
    const double c15_2{7.50};

    void
    calcParticleDistances();

    void
    calcHydroTensors();

    void
    calcHydroForces();

    void
    calcAddedMass();

    void
    calcAddedMassGrad();

    // Helper function to calculate 3D --> 2D horizontal stack tensor reduction in dimension
    static int
    flattenedCol(int realColParticle, int depthParticle, int spatialDim, int dim3N);
};

#endif // BODIES_IN_POTENTIAL_FLOW_POTENTIAL_HYDRODYNAMICS_H