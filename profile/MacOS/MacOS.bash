#!/usr/bin/env bash

# NOTE: Should be called from project base directory

# variables
base_dir="$(pwd)"
output_dir="${base_dir}/profile/MacOS"

# Apple instruments arguments
instruments_template="Time Profiler"
trace_output="${output_dir}/bodies-in-potential-flow-profile.trace"
executable="${base_dir}/build/profile/src/bodies-in-potential-flow"
args="${base_dir}/test/input/collinear_swimmer_isolated/initial_frame_dt1e-2.gsd ${output_dir}/output"

# configure project
cmake -G "Ninja" -B "build/profile" -DCMAKE_BUILD_TYPE=Profile -DENABLE_TESTING=False -DENABLE_COVERAGE=False

# build project
cmake --build "build/profile" -j

# profile project
instruments -t "${instruments_template}" -D ${trace_output} ${executable} ${args}
