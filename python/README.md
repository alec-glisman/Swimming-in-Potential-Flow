# Directory: python

Python scripts to generate GSD files to input to simulation as well as analyze GSD files output from simulation.

`analysis`: Directory contains Python scripts that analyze GSD files (output from C++ simulation).  
`initial_configurations`: Directory contains Python scripts that create GSD files (input to C++ simulation).  

`GSDUtil.py`: Python class to create/load GSD files.  

- Initializer will either create a GSD and output to input path, or load a GSD from a given path based on the value of `create_gsd`
- `setLogParameters()` will output integrator, material, and other parameters.
- `setParticleParameters()` will output configuration and particle type parameters.
- `setKinematics()` will output kinematic data in both float and double data types. 
**Note that the acceleration components are output to moment_inertia as GSD does not store acceleration by default.**
