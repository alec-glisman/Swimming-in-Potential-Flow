# Project Development Progress

## Get project running

- [x] Decide on input data schema
  - May be a good idea to make GSD input in python and then load in C++.
- [x] Write data parsing scheme
  - [x] Kinematic input
  - [x] Parameter input
- [x] Add logging of startup behavior
- [x] Write data output (preferably in GSD file)
- [x] Figure out why input assertions seem to pass, but data not loaded into class??
- [x] Add logging of startup behavior

- [x] Write leading order hydrodynamics calculation

  - [x] Check sign of added mass matrix and gradient

- [x] Write RK4 integration function

  - [x] Add function for t=0 step to correctly set velocity and acceleration
  - [x] Calculate (linear) velocity of CoM
  - [x] Calculate (linear) acceleration of CoM
  - [x] Calculate angular CoM components
  - [x] Update accelerations based on constrained lagrangian approach in paper

- [x] Convert raw pointers into smart pointers
- [x] Figure out why profile build is not linking as expected

- [x] Perl script to run many simulations
- [x] Python analysis script to analyze data and generate plots

## Updated hydrodynamics

- Add quaternions to simulation

  - [x] Update GSD creation to also input quaternions
  - [x] Load quaternions on C++ startup
  - [x] Output quaternions on frame-write

- Move data and calculations to `simulationData` class

  - [x] Load particle-typeID and group-ID on simulation startup (use for body and locater-point information)
  - [x] Make number of bodies a parameter
  - [x] Add DoF vector (body center and orientation) for all kinematic quantities (used in integration and EoM).
  - [x] Move kinematic constraints (velocities and accelerations) from integration to data class
    - [x] Add setter/getter functions (pass as constant)
  - [x] Update kinematic constraints in integration functions at each step (even intermediate steps)
  - [x] Add new data to `update()` function inside the class

- General-use functions for data class

  - [x] Calculate $\bm{E}$ from an input 4-vector.
  - [x] Compute general tensors
    - [x] $\bm \epsilon$
    - [x] $\tilde{\bm{\kappa}}$
    - [x] $\bm{C}^{(i)}$
    - [x] $\bm{\beta}$
    - [x] $\tilde{\bm{\epsilon}}$
  - [x] Calculate rigid body motion tensors
    - [x] $\bm{\Sigma}$
    - [x] $\bm{A}$
    - [x] $\nabla_{\xi} \bm{A}$
  - [x] Move 3rd axis indexing function to data from hydrodynamics
  - [x] Move crossProduct matrix function from hydrodynamics
  - [x] `systemData` add function to calculate all inter-particle distances from relevant locater points
  - [x] `systemData` save body number data for each particle
  - [x] Save initial position data in systemData class for quaternion rotation (normalized).
  - [x] Update integration kinematic calculations with only body components
    - [x] Write function to convert body velocities to particle velocities
    - [x] Write function to convert body accelerations to particle velocities
    - [x] Write function to convert body positions to particle positions

- General-use functions for hydrodynamics

  - [x] Update $\nabla_{R} \bm{M}$ calculation with $\bm{\beta}$
  - [x] Make $\nabla_{R} \bm{M}$ 3rd order tensor
  - [x] Compute the $\bm{N}^{(i)}$ 3rd order tensors and $\tilde{\mathbf{M}}$ 2nd order tensors
  - [x] Update how hydrodynamic forces are calculated

- Integration
  - [x] Figure out how to calculate the initial quaternions
  - [x] Check that udwadia method is being used
  - [x] Make sure particle positions are updated first, then update systemData, then update potentialHydrodynamics
  - [x] Update linear constraint system such that the only constraint is the quaternion unit norm for each body (systemData).
  - [x] Change 2nd order integration to use body components.
  - [x] Rewrite kinematic constraints using quaternions
  - [x] Figure out how to calculate particle positions from body positions simply and with no assumptions.

## Python modules

- [x] Add doxygen comments to all 5 python files
- [x] Add examples to call scripts in the doxygen notation
- [ ] Create/modify files to input to testing scripts in input directory and give GSD files better names. Also maybe use larger dt
- [ ] Change image_sys var
  - [ ] C++: change parsing to type int and then convert to bool after parsing
  - [ ] Python: change type to integer as well

## Code coverage

- [ ] See if `make test` can also run the tests
- [ ] Find out how to analyze code coverage
- [ ] Start code coverage with `simulationData` class
- [ ] Testing class should inherit from parent data class. May need to change private to protected to be able to access them correctly.

## Tensor module optimization

- [x] move thread-pool and device to engine class
- [x] Make sure all calls to TensorCast and MatrixCast have a device to compute on
- [x] Make sure all tensor calculations have a device to compute on
- [x] Optimize all contraction operations with a device
  - [ ] Try using the default device to check performance
- [ ] Add GPU offloading (#define EIGEN_USE_GPU)
  - <https://www.tensorflow.org/guide/create_op#compile_the_op_using_your_system_compiler_tensorflow_binary_installation>
  - See if this can be done with type doubles. I think I read that this only works with floats?
- [ ] Remove all setters/getters in systemData that are unused

## Bugs to Fix

- Figure out why python configuration script is not outputting correct data
  - [x] Quaternions
  - [x] Particle ID type

## Verifications

- [ ] Move assertions to testing classes
- [ ] Check that all quaternions obey the unit-norm requirement
- [ ] Check that the collinear swimmer (isolated) recovers the expected results
  - [ ] Verify collinear constraint upheld for the wall system.
  - [ ] Verify that the collinear swimmer wall system recovers the isolated result when swimmers are very far from wall.

## Additional features

- [x] Make kinematic variables in systemData class private
  - [x] Add getter and setter functions
- [ ] Write WCA potential
- [ ] Convert Perl scripts to Python
- [ ] Separate GSDUtil and systemData
  - [ ] Move GSDUtil creation to engine class
  - [ ] Move systemData creation to engine class
    - [ ] Main constructor now only requires the IO path strings
    - [ ] Make 2nd engine constructor that takes in a pre-made systemData class
  - [ ] Move GSD return values from systemData to GSDUtil
