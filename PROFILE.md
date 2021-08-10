# Profiling Code

## MacOS

Use the instruments app ([documentation](https://help.apple.com/instruments/mac/current/#/devb14ffaa5)). Hard to pass command line arguments that are required. This is made easier by calling instruments via command line rather than GUI.

**Sample:**
Assuming function is called from base directory of project.

```[shell]
instruments_template="Time Profiler"
trace_output="~/Desktop/bodies-in-potential-flow-profile.trace"
executable="$(pwd)/build_profile/src/bodies_in_potential_flow"
args="$(pwd)/temp/test.gsd $(pwd)/temp/output"

instruments -t "${instruments_template}" -D "${trace_output}" "${executable}" "${args}"
```
