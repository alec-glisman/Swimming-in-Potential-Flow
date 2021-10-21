//
// Created by Alec Glisman on 10/20/21
//

#include <testSystemData.hpp>

int
testSystemData::testCrossProdMat()
{
    int num_failed_tests{0};

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

    // test 4: linear combination of Cartesian basis vectors
    const Eigen::Vector3d x_1 = 0.25 * e_x - 30 * e_y + 2.0 * e_z;
    Eigen::Matrix3d       mat_x_1;
    mat_x_1 << 0.0, -x_1(2), x_1(1), x_1(2), 0.0, -x_1(0), -x_1(1), x_1(0), 0.0;

    Eigen::Matrix3d test_mat_x_1;
    m_system->crossProdMat(x_1, test_mat_x_1);
    num_failed_tests += !(test_mat_x_1.isApprox(mat_x_1));

    return num_failed_tests;
}

int
testSystemData::testEMatrix()
{
    int num_failed_tests{0};

    return num_failed_tests;
}
