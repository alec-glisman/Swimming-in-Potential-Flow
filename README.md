# Bodies in Potential Flow

Study dynamics of inertial active matter in a potential fluid (irrotational, incompressible).  

Author: Alec Glisman

## Data I/O

All data is input and output from simulation using the [HOOMD GSD format](https://gsd.readthedocs.io/en/stable/index.html).
The [schema](https://gsd.readthedocs.io/en/stable/python-module-gsd.fl.html) is well-documented.

I will be making a few modifications and make extensive use of the `log` section of the schema to save data I desire.
Many variables are stored as floats, but I output a number of variables as doubles in the logs section for accuracy in further computations.

Frame 0 is created using a Python script and passed into the C++ simulation.
The simulation then updates parameters (such as kinematics), and outputs the "true" initial frame as Frame 1.
For this reason, there could be issues when loading data from frame 0 and any data that is not an input parameter should not be used for further work.

## External Dependencies

### Must install separately

* [Intel MKL](https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/onemkl.html#gs.7owc4e) (oneAPI v2021.1.1): Optimized mathematical instructions
* [Boost](https://www.boost.org/) (v1.76.0): All-purpose STL extension
* [Eigen3](https://gitlab.com/libeigen/eigen) (master 66499f0f): Linear algebra
* [spdlog](https://github.com/gabime/spdlog) (v1.9.1): Logging
* [Catch2](https://github.com/catchorg/Catch2) (v2.13.6): Unit testing

```[shell]
# Install using homebrew (works on MacOS and Linux)
brew install boost spdlog catch2
brew install eigen --HEAD  # Needed for Eigen::seqN()

# Intel OneAPI must be downloaded online
# @SOURCE: https://software.intel.com/content/www/us/en/develop/tools/oneapi/base-toolkit/download.html
```

Note that there are more dependencies for python scripts via a conda environment and perl scripts via cpan. The relevant files to generate these environments are located in the `requirements` directory.

```[shell]
conda env export --from-history  > requirements/Python/environment.yml  # Export conda environment
conda env create -f requirements/Python/environment.yml # Create conda environment

# Install Perl modules via cpanm
cd requirements/Perl
cpanm --installdeps .
```

### Integrated into project

* [gsd](https://github.com/glotzerlab/gsd) (v2.4.2): Simulation data output

## Software tested

* CMake: v3.21.1
* GCC: v11.2.0
* Python (Anaconda): v3.9.5
* Perl: v5.30.0
* ZSH: v5.8
