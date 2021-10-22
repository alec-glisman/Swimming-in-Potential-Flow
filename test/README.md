# Directory: test

C++ code containing unit test for simulation.

`CMakeLists.txt` links C++ files in this directory to the dependencies in the `src` directory.  
`data.gsd` is a sample GSD file that is loaded by the unit test.

`testMain.cpp` simply defines a main file to Catch2.
Commented out code gives a quick example of possible commands.  
`testSimulationBuild.cpp` contains unit test verifying the GSD can be loaded into the simulation and the simulation can initialize free of errors.
