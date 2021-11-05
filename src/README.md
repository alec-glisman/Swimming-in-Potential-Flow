# Directory: src

C++ code for the simulation.

---

## Subdirectory: cuda_helpers

Currently unused. Contains files that help to verify CUDA functions run and exit successfully.

- These files come from the [NVIDIA cuda-samples repository](https://github.com/NVIDIA/cuda-samples)
- Located in `Common` directory

---

## Subdirectory: data_io

C and C++ Files for importing and exporting data from simulation.

### Class: gsd

[HOOMD GSD](https://gsd.readthedocs.io/en/stable/python-module-gsd.hoomd.html) library for direct GSD I/O

### Class: GSDUtil

Wrapper for gsd class that loads data relevant to simulation.

---

## Subdirectory: forces

### Class: potentialHydrodynamics

Calculates hydrodynamic tensors (mass matrices and gradients) and forces for spheres in potential flow.
We only look at the leading-order dipole-dipole interactions $\mathcal{O}(r^{-3})$.
Errors are of $\mathcal{O}(r^{-6})$.

---

## Subdirectory: integrators

### Class: rungeKutta4

Numerical integration of kinematic data over one time step using the public `integrate()` method with a Runge-Kutta 4th order method.
Modify this class for alternate system configurations and constraints.

The private method `accelerationUpdate()` updates the acceleration using the constrained Lagrangian mechanics framework developed by Udwadia and Kalaba in 1996 ([paper](https://royalsocietypublishing.org/doi/pdf/10.1098/rspa.1992.0158?casa_token=FB12tYItZ5sAAAAA:y6w2zNwEaNEnvyNY1DXZVS5f67E4tZ52a0sja6w9TSHFJFbDpKvt9wdgPIuHbHZWCZOOpjb3l8LyPQ), [Wikipedia](https://en.wikipedia.org/wiki/Udwadia%E2%80%93Kalaba_formulation)).

Systems invariant to rigid spatial translation and rotation can also use the private `momentumLinAngFree()` method to calculate the rigid body motion portion of the motion. This assumes the articulation velocity and acceleration is known as well as the locater point. Before calling this method, call `articulationVel()`, `articulationAcc()` and `rLoc()` to update the relevant member variables.

---

## Subdirectory: simulation_system

### Class: engine

The engine class assembles the simulation system and runs the time integration with the public `run()` method.
The constructor also constructs the integrators and forces needed for the dynamics.

### Class: progressBar

`progressBar.hpp` modified version of [prakhar1989/progress-cpp](https://github.com/prakhar1989/progress-cpp.git) that displays simulation progress to terminal during execution.

### Class: SystemData

The SystemData class contains all relevant data for the general simulation and can be accessed through relevant getter and setter functions.
The constructor also constructs the GSD parser class and loads data into itself.

---

## Files

`CMakeLists.txt` links all files and external libraries into an executable named `bodies-in-potential-flow`.  
`main.cpp` main file of the project that constructs and runs the simulation.
Requires two command line inputs: (1) input GSD filepath and (2) output directory to write data, respectively.
The output directory must already exist, as the C++ code will not create it.
