//
// Created by Alec Glisman on 10/20/21
//

#ifndef BODIES_IN_POTENTIAL_FLOW_TEST_SYSTEMDATA_HPP
#define BODIES_IN_POTENTIAL_FLOW_TEST_SYSTEMDATA_HPP

#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif

/* Include all internal project dependencies */
#include <systemData.hpp>

/**
 * @class testSystemData
 *
 * @brief Friend class to test `systemData` class
 *
 */
class testSystemData
{
  public:
    // classes
    /// shared pointer reference to systemData class
    std::shared_ptr<systemData> m_system;

    /**
     * @brief Construct a new test System Data object
     *
     * @param sys systemData class to test
     */
    explicit testSystemData(std::shared_ptr<systemData> sys) : m_system(sys){};

    /**
     * @brief Destroy the test System Data object
     *
     */
    ~testSystemData() = default;

    /**
     * @brief Test `systemData::crossProdMat()`
     *
     * @return int Number of failed tests
     */
    int
    testCrossProdMat();

    /**
     * @brief Test `systemData::eMatrix()`
     *
     * @return int Number of failed tests
     */
    int
    testEMatrix();

  private:
};

#endif // BODIES_IN_POTENTIAL_FLOW_TEST_SYSTEMDATA_HPP
