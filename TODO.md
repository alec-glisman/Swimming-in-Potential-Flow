# To-Do

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
  - [ ] Compute general tensors
    - [x] $\bm \epsilon$
    - [x] $\tilde{\bm{\kappa}}$
    - [ ] $\bm{C}^{(i)}$
    - [ ] $\bm{\beta}$
    - [ ] $\tilde{\bm{\epsilon}}$
  - [ ] Calculate rigid body motion tensors
    - [x] $\bm{\Sigma}$
    - [x] $\bm{A}$
    - [ ] $\nabla_{\xi} \bm{A}$
  - [x] Move 3rd axis indexing function to data from hydrodynamics
  - [x] Move crossProduct matrix function from hydrodynamics

- General-use functions for hydrodynamics

  - [ ] Update $\nabla_{R} \bm{M}$ calculation with $\bm{\beta}$
  - [ ] Make $\nabla_{R} \bm{M}$ 3rd order tensor
  - [ ] Compute the $\bm{N}^{(i)}$ 3rd order tensors and $\tilde{\mathbf{M}}$ 2nd order tensors
  - [ ] Update how hydrodynamic forces are calculated

- Integration
  - [x] Check that udwadia method is being used
  - [ ] Update linear constraint system such that the only constraint is the quaternion unit norm for each body.
  - [ ] move hyper-parameter booleans into `systemData` with clear comment anchors saying how to change (and when)
  - [ ] Update acceleration calculations (only use body DoF and also output particle DoF)

## Tensor module optimization

- [ ] Set a general device to call from systemData class
- [ ] Make sure all calls to TensorCast have a device to compute on
- [ ] Add GPU offloading (#define EIGEN_USE_GPU)
  - <https://www.tensorflow.org/guide/create_op#compile_the_op_using_your_system_compiler_tensorflow_binary_installation>

## Bugs to Fix

- Figure out why python configuration script is not outputting correct data
  - [x] Quaternions
  - [x] Particle ID type

## Verifications

- [ ] Check that all quaternions obey the unit-norm requirement
- [ ] Check that the collinear swimmer (isolated) recovers the expected results
  - [ ] Verify collinear constraint upheld for the wall system.
  - [ ] Verify that the collinear swimmer wall system recovers the isolated result when swimmers are very far from wall.

## Additional features

- [x] Make kinematic variables in systemData class private
  - [x] Add getter and setter functions
- [ ] Write WCA potential
