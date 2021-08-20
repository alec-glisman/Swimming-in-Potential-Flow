# Directory: src

C++ code for the simulation.  

## Subdirectory: data_io

## Subdirectory: forces

## Subdirectory: integrators

## Subdirectory: simulation_system

## Files

`CMakeLists.txt` links all files and external libraries into an executable named `bodies_in_potential_flow`.  
`main.cpp` main file of the project that constructs and runs the simulation.
Requires two command line inputs: (1) input GSD filepath and (2) output directory to write data, respectively.
The output directory must already exist, as the C++ code will not create it.
