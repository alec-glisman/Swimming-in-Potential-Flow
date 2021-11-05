//
// Created by Alec Glisman on 10/20/21
//

#ifndef BODIES_IN_POTENTIAL_FLOW_TEST_SystemData_HPP
#define BODIES_IN_POTENTIAL_FLOW_TEST_SystemData_HPP

#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif

/* Include all internal project dependencies */
#include <SystemData.hpp>

/* Include all external project dependencies */
#include <random> // std::uniform_real_distribution, std::default_random_engine

/**
 * @class TestSystemData
 *
 * @brief Friend class to test `SystemData` class
 *
 */
class TestSystemData
{
  public:
    // classes
    /// shared pointer reference to SystemData class
    std::shared_ptr<SystemData> m_system;

    /**
     * @brief Construct a new test System Data object
     *
     * @param sys SystemData class to test
     */
    explicit TestSystemData(std::shared_ptr<SystemData> sys) : m_system(sys){};

    /**
     * @brief Destroy the test System Data object
     *
     */
    ~TestSystemData() = default;

    /**
     * @brief Test `SystemData::crossProdMat()`
     *
     * @return int Number of failed tests
     */
    int
    testCrossProdMat();

    /**
     * @brief Test `SystemData::eMatrix()`
     *
     * @return int Number of failed tests
     */
    int
    testEMatrix();

    /**
     * @brief Test `SystemData::rigidBodyMotionTensors()`
     *
     * @return int Number of failed tests
     */
    int
    testRigidBodyMotionTensors();

  private:
};

#endif // BODIES_IN_POTENTIAL_FLOW_TEST_SystemData_HPP
