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
// eigen3(Linear algebra)
#define EIGEN_NO_AUTOMATIC_RESIZING
#define EIGEN_USE_MKL_ALL
#define EIGEN_USE_THREADS
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Eigen>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include <eigen3/unsupported/Eigen/CXX11/ThreadPool>
#include <helper_eigenTensorConversion.hpp>
// Logging
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/spdlog.h>

/* Forward declarations */
class systemData;

class potentialHydrodynamics
{
  public:
    potentialHydrodynamics(std::shared_ptr<systemData> sys);

    ~potentialHydrodynamics();

    void
    update();

    void
    updateForcesOnly();

  private:
    // classes
    std::shared_ptr<systemData> m_system;

    // logging
    std::string       m_logFile;
    const std::string m_logName{"potentialHydrodynamics"};

    // For-loop variables
    int m_num_inter{-1}; // Number of pairwise interactions to count

    // tensor variables
    int m_3N{-1}; // length of tensor quantities

    // Distance between particle pairs
    Eigen::VectorXd m_alphaVec; // particle number; dim = (r) x (1)
    Eigen::VectorXd m_betaVec;  // particle number; dim = (r) x (1)
    Eigen::VectorXd m_r_mag_ab; // [1]; dim = (r) x (1)
    Eigen::MatrixXd m_r_ab;     // [1]; dim = (3)  x (r)

    // Identity Matrices
    Eigen::MatrixXd       m_I3N;                                  // dim = (3N) x (3N)
    Eigen::MatrixXd       m_c1_2_I3N;                             // dim = (3N) x (3N)
    const Eigen::Matrix3d m_I3 = Eigen::Matrix3d::Identity(3, 3); // dim = 3 x 3

    // hydrodynamic quantities
    Eigen::MatrixXd          m_M_added;
    Eigen::MatrixXd          m_M_intrinsic;
    Eigen::MatrixXd          m_M_total;
    Eigen::Tensor<double, 3> m_grad_M_added;

    Eigen::VectorXd m_F_hydro;
    Eigen::VectorXd m_F_hydroNoInertia;

    Eigen::VectorXd m_t1_Inertia;
    Eigen::VectorXd m_t2_VelGrad;
    Eigen::VectorXd m_t3_PosGrad;

    // constants
    const double m_unitSphereVol{4.0 / 3.0 * M_PI};
    const double m_c1_2{0.50};
    const double m_c2_3{2.0 / 3.0};
    const double m_c3_2{1.50};
    const double m_c15_2{7.50};

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

  public:
    const Eigen::VectorXd&
    fHydro() const
    {
        return m_F_hydro;
    }

    const Eigen::VectorXd&
    fHydroNoInertia() const
    {
        return m_F_hydroNoInertia;
    }

    const Eigen::MatrixXd&
    mTotal() const
    {
        return m_M_total;
    }

    const Eigen::MatrixXd&
    mIntrinsic() const
    {
        return m_M_intrinsic;
    }

    const Eigen::MatrixXd&
    mAdded() const
    {
        return m_M_added;
    }

    const Eigen::VectorXd&
    t2VelGrad() const
    {
        return m_t2_VelGrad;
    }

    const Eigen::Tensor<double, 3>&
    gradMAdded() const
    {
        return m_grad_M_added;
    }
};

#endif // BODIES_IN_POTENTIAL_FLOW_POTENTIAL_HYDRODYNAMICS_H