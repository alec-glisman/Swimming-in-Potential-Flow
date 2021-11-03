#!/usr/bin/env bash

# NOTE: Should be called from project base directory

# Commonly used perf commands
# perf stat
#   provides overall statistics for common performance events, including instructions executed and clock cycles consumed. Options allow for selection of events other than the default measurement events.
# perf record
#   records performance data into a file, 'perf.data', which can be later analyzed using the 'perf report' command.
# perf report
#   reads and displays the performance data from the 'perf.data' file created by 'perf record'.
# perf list
#   lists the events available on a particular machine. These events will vary based on performance monitoring hardware and software configuration of the system.
# perf top
#   performs a similar function to the 'top' utility. It generates and displays a performance counter profile in realtime.
# perf trace
#   performs a similar function to the 'strace' tool. It monitors the system calls used by a specified thread or process and all signals received by that application.
# perf help
#   displays a complete list of 'perf' commands.
# perf script
#   reads perf.data (created by perf record) and display trace output


# variables
base_dir="$(pwd)"
output_dir="${base_dir}/profile/Linux"
flame_graph_base_dir="${base_dir}/profile/FlameGraph"

# perf arguments
executable="${base_dir}/build_profile/src/./bodies-in-potential-flow"
args="${base_dir}/input/starter_gsd/collinear_swimmer_wall.gsd ${output_dir}/output"

# configure project
cmake -G "Unix Makefiles" -B "build_profile" -DCMAKE_BUILD_TYPE=Profile -DENABLE_TESTING=True -DCMAKE_EXPORT_COMPILE_COMMANDS=True

# build project
cmake --build "build_profile" -j

# profile project
sudo perf record --all-cpus -g "${executable}" ${args}  # default output is 'perf.data'

# visualize results
sudo chmod a+wrx "perf.data"
mv "perf.data" "${flame_graph_base_dir}/perf.data"
cd "${flame_graph_base_dir}"
perf script | ./stackcollapse-perf.pl | ./flamegraph.pl > perf.svg

# cleanup
mv "perf.svg" "${output_dir}/perf.svg"
mv "perf.data" "${output_dir}/perf.data"
cd "${base_dir}"
