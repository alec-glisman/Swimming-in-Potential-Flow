# Profiling Code

## MacOS

Use the instruments app ([documentation](https://help.apple.com/instruments/mac/current/#/devb14ffaa5)). Hard to pass command line arguments that are required. This is made easier by calling instruments via command line rather than GUI.

**Sample:**

```[shell]
instruments_template="Time Profiler"
trace_output="~/Desktop/test.trace"
executable="$(pwd)/build_debug/src/bodies_in_potential_flow"
args="$(pwd)/temp/test.gsd" "$(pwd)/temp/output"

instruments -t ${instruments_template} -D ${trace_output} ${executable} ${args}
```
