//
// Created by Alec Glisman on 10/20/21
//

#include <testSystemData.hpp>

int
testSystemData::testCrossProdMat()
{
    int num_failed_tests{0};

    std::uniform_real_distribution<double> unif(-1000, 1000);
    std::default_random_engine             re;
    double                                 a_random_double = unif(re);

    // test 1: unit vector in x direction
    Eigen::Vector3d e_x;
    e_x << 1.0, 0.0, 0.0;
    Eigen::Matrix3d mat_e_x;
    mat_e_x << 0.0, -e_x(2), e_x(1), e_x(2), 0.0, -e_x(0), -e_x(1), e_x(0), 0.0;

    Eigen::Matrix3d test_mat_e_x;
    m_system->crossProdMat(e_x, test_mat_e_x);
    num_failed_tests += !(test_mat_e_x.isApprox(mat_e_x));

    // test 2: unit vector in y direction
    Eigen::Vector3d e_y;
    e_y << 0.0, 1.0, 0.0;
    Eigen::Matrix3d mat_e_y;
    mat_e_y << 0.0, -e_y(2), e_y(1), e_y(2), 0.0, -e_y(0), -e_y(1), e_y(0), 0.0;

    Eigen::Matrix3d test_mat_e_y;
    m_system->crossProdMat(e_y, test_mat_e_y);
    num_failed_tests += !(test_mat_e_y.isApprox(mat_e_y));

    // test 3: unit vector in z direction
    Eigen::Vector3d e_z;
    e_z << 0.0, 0.0, 1.0;
    Eigen::Matrix3d mat_e_z;
    mat_e_z << 0.0, -e_z(2), e_z(1), e_z(2), 0.0, -e_z(0), -e_z(1), e_z(0), 0.0;

    Eigen::Matrix3d test_mat_e_z;
    m_system->crossProdMat(e_z, test_mat_e_z);
    num_failed_tests += !(test_mat_e_z.isApprox(mat_e_z));

    // test 4: zero vector
    Eigen::Vector3d x_0;
    x_0 << 0.0, 0.0, 0.0;
    Eigen::Matrix3d mat_x_0;
    mat_x_0 << 0.0, -x_0(2), x_0(1), x_0(2), 0.0, -x_0(0), -x_0(1), x_0(0), 0.0;

    Eigen::Matrix3d test_mat_x_0;
    m_system->crossProdMat(x_0, test_mat_x_0);
    num_failed_tests += !(test_mat_x_0.isApprox(mat_x_0));

    // test 4: linear combination of Cartesian basis vectors (triplicate)
    Eigen::Vector3d x_1 = unif(re) * e_x + unif(re) * e_y + unif(re) * e_z;
    Eigen::Matrix3d mat_x_1;
    mat_x_1 << 0.0, -x_1(2), x_1(1), x_1(2), 0.0, -x_1(0), -x_1(1), x_1(0), 0.0;

    Eigen::Matrix3d test_mat_x_1;
    m_system->crossProdMat(x_1, test_mat_x_1);
    num_failed_tests += !(test_mat_x_1.isApprox(mat_x_1));

    x_1 = unif(re) * e_x + unif(re) * e_y + unif(re) * e_z;
    mat_x_1 << 0.0, -x_1(2), x_1(1), x_1(2), 0.0, -x_1(0), -x_1(1), x_1(0), 0.0;
    m_system->crossProdMat(x_1, test_mat_x_1);
    num_failed_tests += !(test_mat_x_1.isApprox(mat_x_1));

    x_1 = unif(re) * e_x + unif(re) * e_y + unif(re) * e_z;
    mat_x_1 << 0.0, -x_1(2), x_1(1), x_1(2), 0.0, -x_1(0), -x_1(1), x_1(0), 0.0;
    m_system->crossProdMat(x_1, test_mat_x_1);
    num_failed_tests += !(test_mat_x_1.isApprox(mat_x_1));

    return num_failed_tests;
}

int
testSystemData::testEMatrix()
{
    int num_failed_tests{0};

    std::uniform_real_distribution<double> unif(-1000, 1000);
    std::default_random_engine             re;
    double                                 a_random_double = unif(re);

    // test 1: quaternion with first component unity
    Eigen::Vector4d qw;
    qw << 1.0, 0.0, 0.0, 0.0;
    Eigen::Matrix<double, 3, 4> E_qw;
    E_qw << -qw(1), qw(0), -qw(3), qw(2), -qw(2), qw(3), qw(0), -qw(1), -qw(3), -qw(2), qw(1), qw(0);

    Eigen::Matrix<double, 3, 4> test_E_qw;
    m_system->eMatrix(qw, test_E_qw);
    num_failed_tests += !(test_E_qw.isApprox(E_qw));

    // test 2: quaternion with first second unity
    Eigen::Vector4d qx;
    qx << 0.0, 1.0, 0.0, 0.0;
    Eigen::Matrix<double, 3, 4> E_qx;
    E_qx << -qx(1), qx(0), -qx(3), qx(2), -qx(2), qx(3), qx(0), -qx(1), -qx(3), -qx(2), qx(1), qx(0);

    Eigen::Matrix<double, 3, 4> test_E_qx;
    m_system->eMatrix(qx, test_E_qx);
    num_failed_tests += !(test_E_qx.isApprox(E_qx));

    // test 3: quaternion with third component unity
    Eigen::Vector4d qy;
    qy << 0.0, 0.0, 1.0, 0.0;
    Eigen::Matrix<double, 3, 4> E_qy;
    E_qy << -qy(1), qy(0), -qy(3), qy(2), -qy(2), qy(3), qy(0), -qy(1), -qy(3), -qy(2), qy(1), qy(0);

    Eigen::Matrix<double, 3, 4> test_E_qy;
    m_system->eMatrix(qy, test_E_qy);
    num_failed_tests += !(test_E_qy.isApprox(E_qy));

    // test 4: quaternion with fourth component unity
    Eigen::Vector4d qz;
    qz << 0.0, 0.0, 0.0, 1.0;
    Eigen::Matrix<double, 3, 4> E_qz;
    E_qz << -qz(1), qz(0), -qz(3), qz(2), -qz(2), qz(3), qz(0), -qz(1), -qz(3), -qz(2), qz(1), qz(0);

    Eigen::Matrix<double, 3, 4> test_E_qz;
    m_system->eMatrix(qz, test_E_qz);
    num_failed_tests += !(test_E_qz.isApprox(E_qz));

    // test 5: quaternion with all zero components
    Eigen::Vector4d q0;
    q0 << 0.0, 0.0, 0.0, 0.0;
    Eigen::Matrix<double, 3, 4> E_q0;
    E_q0 << -q0(1), q0(0), -q0(3), q0(2), -q0(2), q0(3), q0(0), -q0(1), -q0(3), -q0(2), q0(1), q0(0);

    Eigen::Matrix<double, 3, 4> test_E_q0;
    m_system->eMatrix(q0, test_E_q0);
    num_failed_tests += !(test_E_q0.isApprox(E_q0));

    // test 6: linear combinations of basis vectors (triplicate)
    Eigen::Vector4d             q = unif(re) * q + unif(re) * qx + unif(re) * qy + unif(re) * qz;
    Eigen::Matrix<double, 3, 4> E_q;
    E_q << -q(1), q(0), -q(3), q(2), -q(2), q(3), q(0), -q(1), -q(3), -q(2), q(1), q(0);

    Eigen::Matrix<double, 3, 4> test_E_q;
    m_system->eMatrix(q, test_E_q);
    num_failed_tests += !(test_E_q.isApprox(E_q));

    q = unif(re) * q + unif(re) * qx + unif(re) * qy + unif(re) * qz;
    E_q << -q(1), q(0), -q(3), q(2), -q(2), q(3), q(0), -q(1), -q(3), -q(2), q(1), q(0);
    m_system->eMatrix(q, test_E_q);
    num_failed_tests += !(test_E_q.isApprox(E_q));

    q = unif(re) * q + unif(re) * qx + unif(re) * qy + unif(re) * qz;
    E_q << -q(1), q(0), -q(3), q(2), -q(2), q(3), q(0), -q(1), -q(3), -q(2), q(1), q(0);
    m_system->eMatrix(q, test_E_q);
    num_failed_tests += !(test_E_q.isApprox(E_q));

    return num_failed_tests;
}

int
testSystemData::testRigidBodyMotionTensors()
{
    const int num_tests{4};
    int       num_failed_tests{0};

    return num_failed_tests;
}