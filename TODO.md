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
  - [ ] Update GSD creation to also input quaternions
  - [ ] Load quaternions on C++ startup
  - [ ] Output quaternions on frame-write

- Create new hydrodynamics class
  - [ ] Update `CMakeLists.txt` to get new file version in build.
  - [ ] Copy information from previous files

- Move data and calculations to `simulationData` class
  - [ ] Load particle-typeID and group-ID on simulation startup (use for body and locater-point information)
  - [ ] Make number of bodies a parameter
  - [ ] Add DoF vector (body center and orientation) for all kinematic quantities (used in integration and EoM).
  - [ ] Move kinematic constraints (velocities and accelerations) from integration to data class
    - [ ] Add setter/getter functions (pass as constant)
  - [ ] Move mass matrices (and gradient) homes into this class.
  - [ ] Update kinematic constraints in integration functions at each step (even intermediate steps)
  - [ ] Add new data to `update()` function inside the class

- General-use functions for data class
  - [ ] Compute $\bm \epsilon$, $\tilde{\bm{\epsilon}}$, and $\tilde{\bm{\kappa}}$ in the constructor.
  - [ ] Calculate $\bm{E}$ from an input 4-vector.
  - [ ] Calculate rigid body motion tensors ($\bm{\Sigma}$, $\bm{A}$, $\nabla_{\xi} \bm{A}$)
  - [ ] Calculate $\bm{C}^{(i)}$
  - [ ] Calculate $\bm{\beta}$
  - [ ] Move 3rd axis indexing function to data from hydrodynamics
  - [ ] Write function to calculate contraction between 3rd order tensor with 2nd order tensor (both left and right)

- General-use functions for hydrodynamics
  - [ ] Update $\nabla_{R} \bm{M}$ calculation with $\bm{\beta}$
  - [ ] Compute the $\bm{N}^{(i)}$ 3rd order tensors and $\tilde{\mathbf{M}}$ 2nd order tensors
  - [ ] Update how hydrodynamic forces are calculated

- Integration
  - [ ] Check that udwadia method is being used
  - [ ] Update linear constraint system such that the only constraint is the quaternion unit norm for each body.
  - [ ] move hyper-parameter booleans into `systemData` with clear comment anchors saying how to change (and when)
  - [ ] Update acceleration calculations

- Verifications
  - [ ] Check that all quaternions obey the unit-norm requirement
  - [ ] Check that the collinear swimmer (isolated) recovers the expected results
  - [ ] Verify collinear constraint upheld for the wall system.
  - [ ] Verify that the collinear swimmer wall system recovers the isolated result when swimmers are very far from wall.

## Additional features

- [x] Make kinematic variables in systemData class private
  - [x] Add getter and setter functions
- [ ] **Try to speed-up code execution by explicitly declaring size of all matrices at compile time in Eigen types.**
- [ ] Write WCA potential
- [ ] Calculate constraint forces
- [ ] Write RKF45 integration function
