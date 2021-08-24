#!/usr/bin/env bash

# NOTE: Should be called from project base directory

# Commonly used perf commands
# perf stat
#   This command provides overall statistics for common performance events, including instructions executed and clock cycles consumed. Options allow for selection of events other than the default measurement events.
# perf record
#   This command records performance data into a file, 'perf.data', which can be later analyzed using the 'perf report' command.
# perf report
#   This command reads and displays the performance data from the 'perf.data' file created by 'perf record'.
# perf list
#   This command lists the events available on a particular machine. These events will vary based on performance monitoring hardware and software configuration of the system.
# perf top
#   This command performs a similar function to the 'top' utility. It generates and displays a performance counter profile in realtime.
# perf trace
#   This command performs a similar function to the 'strace' tool. It monitors the system calls used by a specified thread or process and all signals received by that application.
# perf help
#   This command displays a complete list of 'perf' commands.


# variables
base_dir="$(pwd)"
output_dir="${base_dir}/profile/Linux"
flame_graph_base_dir="${base_dir}/profile/FlameGraph"

# perf instruments arguments
trace_name="bodies-in-potential-flow-profile.perf"
trace_output="${output_dir}/${trace_name}"
executable="${base_dir}/build_profile/src/bodies_in_potential_flow"
args="${base_dir}/input/starter_gsd/collinear_swimmer_wall.gsd ${output_dir}/output"

# configure project
cmake -G "Unix Makefiles" -B "build_profile" -DCMAKE_BUILD_TYPE=Profile -DENABLE_TESTING=True -DCMAKE_EXPORT_COMPILE_COMMANDS=True

# build project
cmake --build "build_profile" -j

# profile project
perf record --all-cpus -g --output="${trace_output}" "./${executable}" ${args}
cp "${trace_output}" "${flame_graph_base_dir}/${trace_name}"

perf script | "./${flame_graph_base_dir}/stackcollapse-perf.pl" > out.perf_folded
"./${flame_graph_base_dir}/flamegraph.pl" out.perf-folded > perf-kernel.svg