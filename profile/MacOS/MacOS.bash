#!/usr/bin/env bash

# NOTE: Should be called from project base directory

# variables
base_dir="$(pwd)"
output_dir="${base_dir}/profile/MacOS"

# Apple instruments arguments
instruments_template="Time Profiler"
trace_output="${output_dir}/bodies-in-potential-flow-profile.trace"
executable="${base_dir}/build_profile/src/bodies-in-potential-flow"
args="${base_dir}/input/starter_gsd/collinear_swimmer_wall.gsd ${output_dir}/output"

# configure project
cmake -G "Unix Makefiles" -B "build_profile" -DCMAKE_BUILD_TYPE=Profile -DENABLE_TESTING=True -DCMAKE_EXPORT_COMPILE_COMMANDS=True

# build project
cmake --build "build_profile" -j

# profile project
instruments -t "${instruments_template}" -D ${trace_output} ${executable} ${args}
