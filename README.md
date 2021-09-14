# Bodies in Potential Flow

Study dynamics of inertial active matter in a potential fluid (irrotational, incompressible).  

Author: Alec Glisman

## Software tested

* ZSH: $\geq$ v5.8
* Intel oneAPI: $\geq$ v2021.1.1
* CUDA (not currently used): $\geq$ v11.3
* CMake: $\geq$ v3.16.3
* GCC: $\geq$ v11.1.0
* Python (Anaconda): $\geq$ v3.9.5
* Perl: $\geq$ v5.30.0

## Docker integration

The project can be built and run inside a Docker container.
The `docker-compose.yml` file will specify all build options and `Dockerfile` builds the image using a Ubuntu 20.04 base image.
Data I/O is handled automatically with an attached volume.

Commands:

```[shell]
docker-compose --env-file [docker_environment_file] up --detach  # build
docker-compose down         # tear down
```

## Project structure: links to relevant readme files

`.vscode`: [Files relevant for developing the project in VSCode.](.vscode/)
`data_tag`: Data relevant for validating results for a selection of git tags.  
`include`: External dependencies required for the project.
Further information can be found at the end of this readme.  
`input`: Data files that are used in Perl scripts to modify parameters of interest across a range of simulations.  
`profile`: [Scripts to profile code on various platforms and find performance improvement areas.](profile/README.md)  
`python`: [Python scripts to generate GSD files to input to simulation as well as analyze GSD files output from simulation.](python/README.md)  
`requirements`: [Files and scripts relevant for loading C++, Perl, and Python dependencies.](requirements/README.md)  
`scripts`: Perl scripts to run many simulations simultaneously and analyze the results.  
`src`: [C++ code for the simulation.](src/README.md)  
`tests`: [C++ code containing unit tests for simulation.](tests/README.md)

## Data I/O

All data is input and output from simulation using the [HOOMD GSD format](https://gsd.readthedocs.io/en/stable/index.html).
The [schema](https://gsd.readthedocs.io/en/stable/python-module-gsd.fl.html) is well-documented.

I will be making a few modifications and make extensive use of the `log` section of the schema to save data I desire.
Many variables are stored as floats, but I output a number of variables as doubles in the logs section for accuracy in further computations.

Frame 0 is created using a Python script and passed into the C++ simulation.
The simulation then updates parameters (such as kinematics), and outputs the "true" initial frame as Frame 1.
For this reason, there could be issues when loading data from frame 0 and any data that is not an input parameter should not be used for further work.

## Modifications for other systems

The simulation system can be easily adapted for other configurations and constraints.
The code that must be changed inside the C++ framework is tagged with comments of the form `// REVIEW[epic=Change,order=5]: Change constraint linear system for each system`.
The relevant class is the `rungeKutta4` integration class.
Of course, separate Python scripts for GSD initialization and numerical analysis must also be generated.

## External dependencies

Further information found in `requirements` directory readme.

### Must install separately

* [Intel MKL](https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl.html#gs.7owc4e) (oneAPI v2021.1.1): Optimized mathematical instructions
* [Boost](https://www.boost.org/) (v1.71.0): All-purpose STL extension
* [Eigen3](https://gitlab.com/libeigen/eigen) (3.4.0): Linear algebra
* [spdlog](https://github.com/gabime/spdlog) (v1.9.1): Logging
* [Catch2](https://github.com/catchorg/Catch2) (v2.13.6): Unit testing

### Integrated into project

* [gsd](https://github.com/glotzerlab/gsd) (v2.4.2): Simulation data output
