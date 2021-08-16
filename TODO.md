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
  - [ ] Check sign of added mass matrix and gradient

- [x] Write RK4 integration function
  - [ ] Add function for t=0 step to correctly set velocity and acceleration
  - [ ] Calculate (linear) velocity of CoM
  - [ ] Calculate (linear) acceleration of CoM
  - [ ] Calculate angular CoM components

- [x] Convert raw pointers into smart pointers
- [x] Figure out why profile build is not linking as expected

- [x] Perl script to run many simulations
- [x] Python analysis script to analyze data and generate plots

## Additional features

- [ ] Write WCA potential
- [ ] Calculate constraint forces
- [ ] Write RKF45 integration function
