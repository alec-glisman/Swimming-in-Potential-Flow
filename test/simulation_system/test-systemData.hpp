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

class testSystemData : public systemData
{
};

#endif // BODIES_IN_POTENTIAL_FLOW_TEST_SYSTEMDATA_HPP
