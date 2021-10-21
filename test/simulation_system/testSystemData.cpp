//
// Created by Alec Glisman on 10/20/21
//

#include <testSystemData.hpp>

testSystemData::testSystemData(std::shared_ptr<systemData> sys)
{
    m_system = sys;
}

testSystemData::~testSystemData()
{
}

int
testSystemData::testCrossProdMat()
{
    int num_failed_tests{0};

    // testing inputs
    Eigen::Vector3d e_x;
    e_x << 1.0, 0.0, 0.0;
    Eigen::Matrix3d mat_e_x;
    mat_e_x << 0.0, -e_x(2), e_x(1), e_x(2), 0.0, -e_x(0), -e_x(1), e_x(0), 0.0;

    Eigen::Vector3d e_y;
    e_y << 0.0, 1.0, 0.0;

    Eigen::Vector3d e_z;
    e_z << 0.0, 0.0, 1.0;

    Eigen::Vector3d x_0;
    x_0 << 0.0, 0.0, 0.0;

    const Eigen::Vector3d x_1 = 0.25 * e_x - 30 * e_y + 2.0 * e_z;

    // test outputs
    Eigen::Matrix3d test_mat_e_x;
    m_system->crossProdMat(e_x, test_mat_e_x);
    num_failed_tests += !(test_mat_e_x.isApprox(mat_e_x));

    return num_failed_tests;
}
