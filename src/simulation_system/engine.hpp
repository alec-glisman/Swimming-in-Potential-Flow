//
// Created by Alec Glisman on 07/31/21
//

#ifndef BODIES_IN_POTENTIAL_FLOW_ENGINE_H
#define BODIES_IN_POTENTIAL_FLOW_ENGINE_H

#ifdef NVCC
#error This header cannot be compiled by nvcc
#endif

/* Include all internal project dependencies */
#include <GSDUtil.hpp>
#include <gsd.h>
#include <systemData.hpp>

/* Include all external project dependencies */
// Logging
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/spdlog.h>
// STL
#include <memory>    // for std::unique_ptr and std::shared_ptr
#include <stdexcept> // std::errors
#include <string>    // std::string

/* Forward declarations */
class systemData;
class GSDUtil;

class engine
{
  public:
  private:
};

#endif // BODIES_IN_POTENTIAL_FLOW_ENGINE_H