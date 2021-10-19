# Directory: tests

C++ code containing unit tests for simulation.

`CMakeLists.txt` links C++ files in this directory to the dependencies in the `src` directory.  
`data.gsd` is a sample GSD file that is loaded by the unit tests.

`testsMain.cpp` simply defines a main file to Catch2.
Commented out code gives a quick example of possible commands.  
`testsSimulationBuild.cpp` contains unit tests verifying the GSD can be loaded into the simulation and the simulation can initialize free of errors.
